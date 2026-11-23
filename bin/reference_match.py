#!/usr/bin/env python3
"""
REFERENCE_MATCH: find the current run's clonotypes in external reference clonotype sets
(earlier cohorts or trial iterations, published tumour-reactive or antigen-specific clones).

Reference input (--reference): either
  * a reference samplesheet CSV with a `file` column, one row per reference table:
        study,file[,cdr3_col,v_col,count_col,freq_col,label_cols,patient_col]
    `study` names the reference set (default: the file name); the optional *_col columns
    name that table's columns when they are not auto-detected; or
  * a single reference table (CSV or TSV), used as one study named after the file.

Reference tables are AIRR rearrangement tables or any clone list. Columns are auto-detected
(first match wins, case-insensitive):
    CDR3 amino acid   junction_aa, cdr3_aa, cdr3aa, cdr3b, cdr3, aminoAcid, TCR_ Clonotype
    V gene (optional) v_call, v_gene, vGeneName, v
    count (optional)  duplicate_count, templates, count, cloneCount
    freq (optional)   frequency, freq, productive_frequency, cloneFraction
    patient (opt.)    patient, subject_id, subject, donor
    labels (optional) antigen, epitope, label, putative_target, Putative Target, species
Only the CDR3 column is required.

Matching (--match-key): exact match on `junction_aa,v_call` (default) or `junction_aa`.
V genes are compared at gene level (allele dropped; Adaptive names such as TCRBV05-01 are
mapped to IMGT TRBV5-1). A reference without a V column is matched on CDR3 alone, and the
`match_key` column of every output says which key was used.
# TODO: similarity-based matching (TCRdist3 radius or GLIPH2 motif membership). It needs a
# cross-dataset distance step (reference x sample tcrdist3 rectangular matrix, or a GLIPH2
# run on the pooled sets) that the pipeline does not have yet; only exact matching is done.

Outputs:
  reference_matches.tsv            every (sample, reference clone) match, with the clone's
                                   count and frequency in both studies
  per_sample/<sample>_reference_matches.tsv   the same, split per sample
  reference_overlap_summary.tsv    per sample and per patient x study: clone counts, shared
                                   clones, fraction of clones / repertoire shared, Jaccard and
                                   Morisita-Horn (the latter only when the reference has counts
                                   or frequencies)
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

CANDIDATES = {
    "cdr3": ["junction_aa", "cdr3_aa", "cdr3aa", "cdr3b", "cdr3_b", "cdr3", "aminoacid",
             "amino_acid", "tcr_ clonotype", "tcr_clonotype", "clonotype"],
    "v": ["v_call", "v_gene", "vgenename", "v_gene_name", "trbv", "v"],
    "count": ["duplicate_count", "templates", "count", "counts", "clonecount", "clone_count"],
    "freq": ["frequency", "freq", "productive_frequency", "clonefraction", "clone_fraction"],
    "patient": ["patient", "patient_id", "subject_id", "subject", "donor"],
    "labels": ["antigen", "epitope", "antigen.epitope", "label", "putative_target", "putative target",
               "species", "antigen.species", "category"],
}
PATIENT_CANDIDATES = ("patient", "patient_id", "subject_id", "subject", "donor")

MATCH_COLS = ["sample", "patient", "study", "match_key", "junction_aa", "v_gene",
              "count", "freq", "ref_count", "ref_freq", "ref_patient", "ref_labels"]
SUMMARY_COLS = ["level", "id", "patient", "study", "match_key", "clones", "reference_clones",
                "shared_clones", "frac_clones_shared", "frac_repertoire_shared",
                "frac_reference_clones_found", "frac_reference_repertoire_found",
                "jaccard", "morisita_horn"]


def v_gene(v) -> str:
    """Gene-level IMGT V call: first of an ambiguous list, allele dropped, Adaptive names mapped."""
    if pd.isna(v) or str(v).strip() == "":
        return ""
    s = str(v).split(",")[0].split("*")[0].strip()
    m = re.match(r"^TCR([ABDG])V0*(\d+)(?:-0*(\d+))?", s, re.I)
    if m:  # Adaptive: TCRBV05-01 -> TRBV5-1
        s = f"TR{m.group(1).upper()}V{m.group(2)}" + (f"-{m.group(3)}" if m.group(3) else "")
    return s


def clean_cdr3(x) -> str:
    return re.sub(r"[^A-Z]", "", str(x).upper()) if pd.notna(x) else ""


def find_col(df: pd.DataFrame, kind: str, override: str = ""):
    if override:
        if override not in df.columns:
            raise SystemExit(f"[REFERENCE_MATCH] column '{override}' not found; available: {list(df.columns)}")
        return override
    lower = {c.strip().lower(): c for c in df.columns}
    return next((lower[c] for c in CANDIDATES[kind] if c in lower), None)


def read_table(path: str) -> pd.DataFrame:
    sep = "\t" if path.endswith((".tsv", ".txt", ".tsv.gz")) else ","
    return pd.read_csv(path, sep=sep, dtype=str, low_memory=False)


def load_reference(row: dict, base_dir: str) -> pd.DataFrame:
    # As given, relative to the reference samplesheet, or (Nextflow-staged) by basename next to it.
    tries = [row["file"], os.path.join(base_dir, row["file"]),
             os.path.join(base_dir, os.path.basename(row["file"]))]
    path = next((p for p in tries if os.path.exists(p)), None)
    if path is None:
        raise SystemExit(f"[REFERENCE_MATCH] reference file not found: {row['file']}")
    raw = read_table(path)
    cdr3 = find_col(raw, "cdr3", row.get("cdr3_col", ""))
    if cdr3 is None:
        raise SystemExit(f"[REFERENCE_MATCH] {path}: no CDR3 column found (set cdr3_col)")
    v = find_col(raw, "v", row.get("v_col", ""))
    cnt = find_col(raw, "count", row.get("count_col", ""))
    frq = find_col(raw, "freq", row.get("freq_col", ""))
    pat = find_col(raw, "patient", row.get("patient_col", ""))
    if row.get("label_cols"):
        labels = [c.strip() for c in row["label_cols"].split(";") if c.strip()]
        missing = [c for c in labels if c not in raw.columns]
        if missing:
            raise SystemExit(f"[REFERENCE_MATCH] {path}: label column(s) {missing} not found")
    else:
        lab = find_col(raw, "labels")
        labels = [lab] if lab else []

    ref = pd.DataFrame({
        "junction_aa": raw[cdr3].map(clean_cdr3),
        "v_gene": raw[v].map(v_gene) if v else "",
        "ref_count": pd.to_numeric(raw[cnt], errors="coerce") if cnt else np.nan,
        "ref_freq": pd.to_numeric(raw[frq], errors="coerce") if frq else np.nan,
        "ref_patient": raw[pat].astype(str) if pat else "",
    })
    if labels:
        ref["ref_labels"] = raw[labels].fillna("").astype(str).apply(
            lambda r: " | ".join(re.sub(r"\s+", " ", x).strip() for x in r if x.strip()), axis=1)
    else:
        ref["ref_labels"] = ""
    ref = ref[ref["junction_aa"] != ""]
    ref["study"] = row.get("study") or os.path.splitext(os.path.basename(path))[0]
    ref["has_v"] = bool(v) and (ref["v_gene"] != "").any()
    if ref["ref_freq"].isna().all() and ref["ref_count"].notna().any():
        ref["ref_freq"] = ref["ref_count"] / ref["ref_count"].sum()
    print(f"[REFERENCE_MATCH] {ref['study'].iloc[0] if len(ref) else path}: {len(ref)} clones "
          f"(cdr3={cdr3}, v={v}, count={cnt}, freq={frq}, patient={pat}, labels={labels})")
    return ref


def morisita_horn(a: np.ndarray, b: np.ndarray) -> float:
    if np.nansum(a) == 0 or np.nansum(b) == 0:
        return float("nan")
    pa, pb = a / a.sum(), b / b.sum()
    den = (pa ** 2).sum() + (pb ** 2).sum()
    return float(2 * (pa * pb).sum() / den) if den > 0 else float("nan")


def summarise(cur: pd.DataFrame, ref: pd.DataFrame, key: list, level: str, ident: str,
              patient: str, study: str, key_name: str) -> dict:
    """cur: one sample's (or patient's) clones keyed by `key`, with count/freq."""
    c = cur.groupby(key, as_index=False)["count"].sum()
    c["freq"] = c["count"] / c["count"].sum()
    # min_count=1: a reference without abundances keeps NaN rather than summing to 0.
    r = ref.groupby(key, as_index=False).agg(ref_count=("ref_count", lambda s: s.sum(min_count=1)),
                                             ref_freq=("ref_freq", lambda s: s.sum(min_count=1)))
    m = c.merge(r, on=key, how="inner")
    n_c, n_r, n_s = len(c), len(r), len(m)
    has_abund = r["ref_freq"].notna().any()
    if has_abund:
        allk = c.merge(r, on=key, how="outer")
        mh = morisita_horn(allk["count"].fillna(0).to_numpy(float), allk["ref_freq"].fillna(0).to_numpy(float))
        ref_rep = m["ref_freq"].sum() / r["ref_freq"].sum() if r["ref_freq"].sum() else np.nan
    else:
        mh, ref_rep = np.nan, np.nan
    return {"level": level, "id": ident, "patient": patient, "study": study, "match_key": key_name,
            "clones": n_c, "reference_clones": n_r, "shared_clones": n_s,
            "frac_clones_shared": n_s / n_c if n_c else np.nan,
            "frac_repertoire_shared": float(m["freq"].sum()),
            "frac_reference_clones_found": n_s / n_r if n_r else np.nan,
            "frac_reference_repertoire_found": ref_rep,
            "jaccard": n_s / (n_c + n_r - n_s) if (n_c + n_r - n_s) else np.nan,
            "morisita_horn": mh}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--concat", required=True, help="concatenated_cdr3_sorted.tsv of this run")
    ap.add_argument("--samplesheet", required=True, help="this run's samplesheet CSV")
    ap.add_argument("--reference", required=True, help="reference samplesheet CSV or a single reference table")
    ap.add_argument("--subject-col", default="patient")
    ap.add_argument("--match-key", default="junction_aa,v_call", choices=["junction_aa,v_call", "junction_aa"])
    args = ap.parse_args()

    # ── References ──
    first = read_table(args.reference)
    base_dir = os.path.dirname(os.path.abspath(args.reference))
    rows = (first.fillna("").to_dict("records") if "file" in first.columns
            else [{"file": args.reference, "study": ""}])
    refs = [load_reference(r, base_dir) for r in rows]

    # ── This run ──
    meta = pd.read_csv(args.samplesheet, dtype=str).fillna("")
    subj = next((c for c in (args.subject_col,) + PATIENT_CANDIDATES if c in meta.columns), None)
    cur = pd.read_csv(args.concat, sep="\t", usecols=["sample", "junction_aa", "v_call", "duplicate_count"],
                      dtype={"sample": str, "junction_aa": str, "v_call": str})
    cur["junction_aa"] = cur["junction_aa"].map(clean_cdr3)
    cur = cur[cur["junction_aa"] != ""]
    cur["v_gene"] = cur["v_call"].map(v_gene)
    cur = cur.rename(columns={"duplicate_count": "count"})
    cur["patient"] = cur["sample"].map(meta.set_index("sample")[subj]) if subj else ""
    cur["patient"] = cur["patient"].fillna("")
    cur = cur.groupby(["sample", "patient", "junction_aa", "v_gene"], as_index=False)["count"].sum()
    cur["freq"] = cur["count"] / cur.groupby("sample")["count"].transform("sum")

    want_v = args.match_key == "junction_aa,v_call"
    matches, summary = [], []
    for ref in refs:
        if ref.empty:
            continue
        study = ref["study"].iloc[0]
        use_v = want_v and bool(ref["has_v"].iloc[0])
        if want_v and not use_v:
            print(f"[REFERENCE_MATCH] {study}: no V gene column - matching on CDR3 only", file=sys.stderr)
        key = ["junction_aa", "v_gene"] if use_v else ["junction_aa"]
        key_name = "junction_aa,v_call" if use_v else "junction_aa"

        # Reference clones collapsed on the key (labels / patients of duplicates joined).
        rk = ref.groupby(key, as_index=False).agg(
            ref_count=("ref_count", lambda s: s.sum(min_count=1)),
            ref_freq=("ref_freq", lambda s: s.sum(min_count=1)),
            ref_patient=("ref_patient", lambda s: ";".join(sorted({x for x in s if x}))),
            ref_labels=("ref_labels", lambda s: " ; ".join(dict.fromkeys(x for x in s if x))))
        # On a CDR3-only key, v_gene is the current run's V gene.
        m = cur.merge(rk, on=key, how="inner")
        m["study"], m["match_key"] = study, key_name
        matches.append(m)

        for (sample, patient), sdf in cur.groupby(["sample", "patient"]):
            summary.append(summarise(sdf, ref, key, "sample", sample, patient, study, key_name))
        for patient, pdf in cur.groupby("patient"):
            if patient:
                summary.append(summarise(pdf, ref, key, "patient", patient, patient, study, key_name))

    mt = (pd.concat(matches, ignore_index=True) if matches else pd.DataFrame()).reindex(columns=MATCH_COLS)
    mt = mt.sort_values(["study", "sample", "freq"], ascending=[True, True, False])
    mt.to_csv("reference_matches.tsv", sep="\t", index=False)
    os.makedirs("per_sample", exist_ok=True)
    for sample in cur["sample"].unique():
        mt[mt["sample"] == sample].to_csv(f"per_sample/{sample}_reference_matches.tsv", sep="\t", index=False)
    pd.DataFrame(summary).reindex(columns=SUMMARY_COLS).to_csv("reference_overlap_summary.tsv", sep="\t", index=False)
    print(f"[REFERENCE_MATCH] {len(mt)} match row(s) across {mt['sample'].nunique()} sample(s); "
          f"{len(refs)} reference set(s)")


if __name__ == "__main__":
    sys.exit(main())
