#!/usr/bin/env python3
"""
RESULTS_EXPORT: the run's results as a few tidy tables, written whether or not the report
notebooks render (--run_reports false).

Every figure in the reports is a view of one of these tables, so the set is kept small and
non-redundant: one row per sample, one per sample x clone, one per clone test, one per
sample pair and, with single-cell data, one per cell. Clone-level annotations (generation
probability, VDJdb match, cluster labels, clonal dynamics call, publicity) are columns of
the clone table rather than separate files.

    samples.tsv              sample metadata + repertoire metrics
    clones.tsv.gz            sample x TCRbeta clonotype: counts, frequency, annotations
    clone_tests.tsv.gz       Fisher tests between timepoints (expansion and contraction)
    sample_overlap.tsv       pairwise sample overlap (Jaccard, Sorensen, Morisita-Horn)
    cells.tsv.gz             (single-cell) one row per cell: state, clone, UMAP, per-cell
                             CoNGA / TCRi / cluster calls
    marker_expression.tsv    (GEX) mean expression of marker genes per state and clone group
    README.md                column dictionary and the table behind each report figure

Inputs are all optional; a table whose source did not run on this route is not written.
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd


def log(msg):
    print(f"[results_export] {msg}", file=sys.stderr)


def usable(p):
    return p and os.path.exists(p) and os.path.getsize(p) > 0 and os.path.basename(p) != "NO_FILE"


def read(p, **kw):
    if not usable(p):
        return None
    sep = "," if p.endswith(".csv") else "\t"
    try:
        d = pd.read_csv(p, sep=sep, low_memory=False, **kw)
    except Exception as exc:  # a malformed optional input must not sink the export
        log(f"WARN could not read {p}: {exc}")
        return None
    return d if len(d) else None


def find(d, name):
    """First file called `name` under directory d (inputs are staged by module)."""
    if not d:
        return None
    hits = sorted(glob.glob(os.path.join(d, "**", name), recursive=True))
    return hits[0] if hits else None


def first_col(df, cands):
    return next((c for c in cands if c in df.columns), None)


def v_gene(s):
    return s.astype("string").str.split("*").str[0]


# ── clones ───────────────────────────────────────────────────────────────────
def build_clones(a):
    clono = read(find(a.indir, "clonotypes.csv"))
    if clono is None:
        clono = read(find(a.indir, "concatenated_cdr3_sorted.tsv"))
        if clono is None:
            return None, None
        log("clonotypes.csv absent (no compare stage): clones from concatenated_cdr3_sorted.tsv")
        sheet = read(a.samplesheet) if a.samplesheet else None
        if sheet is not None and "sample" in sheet.columns:
            meta = [c for c in ("sample", "alias", "patient", "timepoint", "origin") if c in sheet.columns]
            clono = clono.merge(sheet[meta].drop_duplicates("sample"), on="sample", how="left")
        tot = clono.groupby("sample")["duplicate_count"].transform("sum")
        clono["frequency"] = clono["duplicate_count"] / tot

    pcol = first_col(clono, [a.patient_col, "patient", "patient_id", "subject_id"])
    keep = {"sample": "sample", pcol: "patient", "origin": "origin", "alias": "alias",
            "timepoint": "timepoint", "timepoint_rank": "timepoint_rank",
            "junction_aa": "junction_aa", "v_call": "v_call", "j_call": "j_call",
            "duplicate_count": "count", "frequency": "frequency"}
    keep = {k: v for k, v in keep.items() if k and k in clono.columns}
    c = clono[list(keep)].rename(columns=keep)
    if "patient" not in c:
        c["patient"] = "all"
    if "origin" not in c:
        c["origin"] = "all"
    c["v_gene"] = v_gene(c["v_call"]) if "v_call" in c else pd.NA
    # Sample metadata (one row per sample) before clone annotations widen the table.
    smeta = c.drop_duplicates("sample")[[x for x in ("sample", "alias", "patient", "timepoint",
                                                     "timepoint_rank", "origin") if x in c.columns]]

    # Generation probability (OLGA), per CDR3b.
    pg = read(find(a.indir, "olga_pgen.tsv"))
    if pg is not None and "junction_aa" in pg.columns:
        pg = pg.drop_duplicates("junction_aa")[[x for x in ("junction_aa", "pgen", "log10_pgen") if x in pg.columns]]
        c = c.merge(pg, on="junction_aa", how="left")

    # VDJdb: best-scoring match per CDR3b over all samples (as the reports pick it).
    vfiles = sorted(glob.glob(os.path.join(a.indir, "**", "*.vdjmatch.txt"), recursive=True))
    if vfiles:
        parts = []
        for f in vfiles:
            d = read(f)
            if d is not None and "cdr3aa" in d.columns:
                parts.append(d[[x for x in ("cdr3aa", "antigen.species", "antigen.gene", "antigen.epitope",
                                            "mhc.class", "vdjdb.score") if x in d.columns]])
        if parts:
            v = pd.concat(parts, ignore_index=True)
            if "vdjdb.score" in v:
                v = v.sort_values("vdjdb.score", ascending=False)
            v = v.drop_duplicates("cdr3aa").rename(columns={
                "cdr3aa": "junction_aa", "antigen.species": "vdjdb_species", "antigen.gene": "vdjdb_gene",
                "antigen.epitope": "vdjdb_epitope", "mhc.class": "vdjdb_mhc_class", "vdjdb.score": "vdjdb_score"})
            c = c.merge(v, on="junction_aa", how="left")

    # Clustering labels, per patient x CDR3b, from the export the consensus step annotated.
    ex = read(a.cluster_export)
    if ex is not None:
        cd = first_col(ex, ["CDR3b", "cdr3b", "junction_aa", "cdr3_b_aa"])
        pc = first_col(ex, [a.patient_col, "patient", "patient_id", "subject_id"])
        lab = [x for x in ("giana_cluster", "gliph_cluster", "tcrdist_cluster", "consensus_cluster",
                           "consensus_methods") if x in ex.columns]
        if cd and pc and lab:
            m = ex[[pc, cd] + lab].rename(columns={pc: "patient", cd: "junction_aa"})
            m = m.dropna(subset=["junction_aa"])
            m["patient"] = m["patient"].astype(str)
            # One label set per clone: the most frequent non-empty one.
            m = (m.replace("", np.nan).groupby(["patient", "junction_aa"], as_index=False)
                   .agg({x: lambda s: s.dropna().mode().iat[0] if s.notna().any() else pd.NA for x in lab}))
            c["patient"] = c["patient"].astype(str)
            c = c.merge(m, on=["patient", "junction_aa"], how="left")

    # Clonal dynamics call, per patient x compartment x clonotype.
    calls = []
    for f, lab in (("expanded_clones_summary.csv", "expanding"), ("contracting_clones_summary.csv", "contracting")):
        d = read(find(a.indir, f))
        if d is None:
            continue
        pc = first_col(d, [a.patient_col, "patient", "patient_id", "subject_id"])
        k = [x for x in ("origin", "junction_aa", "v_call", "j_call") if x in d.columns]
        d = d[[pc] + k].rename(columns={pc: "patient"}).drop_duplicates()
        d["patient"] = d["patient"].astype(str)
        d["dynamics_call"] = lab
        calls.append(d)
    if calls:
        d = pd.concat(calls, ignore_index=True)
        k = [x for x in ("patient", "origin", "junction_aa", "v_call", "j_call") if x in d.columns and x in c.columns]
        d = d.groupby(k, as_index=False)["dynamics_call"].agg(
            lambda s: "both" if s.nunique() > 1 else s.iat[0])
        c["patient"] = c["patient"].astype(str)
        c = c.merge(d, on=k, how="left")
        c["dynamics_call"] = c["dynamics_call"].fillna("not called")

    # Publicity: patients and samples sharing the CDR3b.
    c["n_patients_with_cdr3"] = c.groupby("junction_aa")["patient"].transform("nunique")
    c["n_samples_with_cdr3"] = c.groupby("junction_aa")["sample"].transform("nunique")

    sort = [x for x in ("patient", "origin", "timepoint_rank", "sample") if x in c.columns]
    c = c.sort_values(sort + ["count"], ascending=[True] * len(sort) + [False])
    return c, smeta


# ── samples ──────────────────────────────────────────────────────────────────
def build_samples(a, smeta):
    st = read(find(a.indir, "sample_stats.csv"))
    if smeta is None and st is None:
        return None
    if st is None:
        return smeta
    if smeta is None:
        return st
    st = st.drop(columns=[x for x in smeta.columns if x in st.columns and x != "sample"])
    return smeta.merge(st, on="sample", how="outer")


# ── clone tests ──────────────────────────────────────────────────────────────
def build_tests(a):
    parts = []
    for f, direction in (("detailed_comparisons_expansion.csv", "expansion"),
                         ("detailed_comparisons_contraction.csv", "contraction")):
        d = read(find(a.indir, f))
        if d is not None:
            d.insert(0, "test", direction)
            parts.append(d)
    return pd.concat(parts, ignore_index=True) if parts else None


# ── cells (single-cell) ──────────────────────────────────────────────────────
def build_cells(a):
    cells = read(a.cells)
    if cells is None:
        return None
    cells["cell_id"] = cells["cell_id"].astype(str)
    for p, cols in ((a.cluster_export, ["giana_cluster", "gliph_cluster", "tcrdist_cluster",
                                        "consensus_cluster", "consensus_methods"]),
                    (a.conga_cells, ["conga_score", "conga_gex_cluster", "conga_tcr_cluster", "is_invariant"]),
                    (a.tcri_cells, ["tcri_score", "tcri_group", "tcri_predicted_phenotype"])):
        d = read(p)
        if d is None or "cell_id" not in d.columns:
            continue
        cols = [x for x in cols if x in d.columns and x not in cells.columns]
        if not cols:
            continue
        d = d[["cell_id"] + cols].drop_duplicates("cell_id")
        d["cell_id"] = d["cell_id"].astype(str)
        cells = cells.merge(d, on="cell_id", how="left")
    # Samplesheet metadata (alias, compartment, timepoint) the GEX object may lack.
    sheet = read(a.samplesheet) if a.samplesheet else None
    if sheet is not None and "sample" in sheet.columns:
        add = [x for x in ("alias", "origin", "timepoint") if x in sheet.columns and x not in cells.columns]
        if add:
            cells = cells.merge(sheet[["sample"] + add].drop_duplicates("sample"), on="sample", how="left")
    return cells


# ── README ───────────────────────────────────────────────────────────────────
README = """# TCRtoolkit results tables

Written by RESULTS_EXPORT on every run, with or without the report notebooks
(`--run_reports false`). Each report figure is a view of one of these tables; the reports'
own per-figure TSVs are not written by default (`--write_figure_tables true` restores them).

| Table | One row per | Present |
|---|---|---|
{inventory}

## Columns

**samples.tsv**: `sample`, `alias`, `patient`, `timepoint`, `timepoint_rank`, `origin`
(compartment), then the repertoire metrics of `sample_stats.csv`: `num_clones`, `num_TCRs`,
`clonality`, `shannon_entropy`, `simpson_index(_corrected)`, `inverse_simpson`,
`gini_coefficient`, `hill_q0/q1/q2`, productive / convergent counts, `clone_unit`.

**clones.tsv.gz**: one TCRbeta clonotype (`junction_aa` + `v_call` + `j_call`) in one sample:
`count` (cells, templates or reads), `frequency` (of the sample), sample metadata, `v_gene`,
`pgen` / `log10_pgen` (OLGA), `vdjdb_species` / `vdjdb_gene` / `vdjdb_epitope` /
`vdjdb_mhc_class` / `vdjdb_score` (best VDJdb match), `giana_cluster`, `gliph_cluster`,
`tcrdist_cluster`, `consensus_cluster`, `consensus_methods` (patient-scoped cluster labels),
`dynamics_call` (expanding / contracting / both / not called, between any two timepoints of
the patient and compartment), `n_patients_with_cdr3`, `n_samples_with_cdr3` (publicity).

**clone_tests.tsv.gz**: one Fisher's exact test of one clonotype between two timepoints
(`t_pre`, `t_post`) of a patient and compartment: `test` (expansion or contraction, each
one-sided), counts and frequencies before and after, `fold_change`, `p_value`, `q_value`
(BH within the comparison).

**sample_overlap.tsv**: one pair of samples: `n_shared`, `jaccard`, `sorensen`,
`morisita_horn`, on the run's clone unit.

**cells.tsv.gz** (single-cell): `cell_id`, `sample`, `patient`, `timepoint`, `alias`, `origin`,
`annot` (cell state), `clone_id` (paired alpha/beta), `cdr3a`, `cdr3b`, V/J genes, `has_tcr`,
`paired_tcr`, `U_1` / `U_2` (UMAP), per-cell cluster labels, `conga_score` (< 1 = hit),
`conga_gex_cluster`, `conga_tcr_cluster`, `tcri_score`, `tcri_group`,
`tcri_predicted_phenotype`.

**marker_expression.tsv** (GEX): `gene`, `annot`, `group` (`expanded`: clone of >= 2 cells
in its sample; `singleton`; `no_tcr`), `n_cells`, `mean_logexpr` (mean log-normalised
expression), `frac_expressing`.

## Which table each report figure comes from

| Figure | Table | How |
|---|---|---|
| Clonality, diversity, Hill numbers per sample | samples.tsv | one value per sample |
| Clone-size bins, rank-abundance | clones.tsv.gz | `frequency` per sample |
| Repertoire overlap heatmaps | sample_overlap.tsv | pivot `sample_a` x `sample_b` |
| Persistent / unique clones (UpSet) | clones.tsv.gz | presence of `junction_aa` per patient's samples |
| Clonal dynamics status bars, volcano | clone_tests.tsv.gz | `q_value < 0.05` by `test` |
| Top clones over time (Master Summary T1) | clones.tsv.gz | per patient and compartment, top `frequency` peaks; key = `junction_aa` + `v_gene` |
| TCRdist3 / GIANA / GLIPH2 / consensus clusters | clones.tsv.gz | group by `*_cluster` |
| VDJdb specificity | clones.tsv.gz | `vdjdb_species` |
| Publicity, generation probability | clones.tsv.gz | `n_patients_with_cdr3`, `log10_pgen` |
| Cell-state maps, composition, expansion by state | cells.tsv.gz | `U_1`, `U_2`, `annot`, clone size per `sample` x `clone_id` |
| Top clones' cell states over time, on the map (Master Summary) | cells.tsv.gz | key = `cdr3b` + `trbv`, share of the sample's TCR+ cells |
| Cell states by clonal dynamics | cells.tsv.gz + clones.tsv.gz | join `patient`, `cdr3b` = `junction_aa`, use `dynamics_call` |
| CoNGA / TCRi views | cells.tsv.gz | `conga_score`, `tcri_score`, `tcri_predicted_phenotype` |
| Marker genes by state; expanded vs singleton | marker_expression.tsv | pivot `gene` x `annot`, difference of `group` means |
"""


def main():
    ap = argparse.ArgumentParser(description="Write the run's core results tables.")
    ap.add_argument("--indir", default="in", help="directory holding the staged module tables")
    ap.add_argument("--cluster_export", default=None, help="export annotated with cluster labels")
    ap.add_argument("--cells", default=None, help="per-cell T-cell export (single-cell GEX route)")
    ap.add_argument("--conga_cells", default=None)
    ap.add_argument("--tcri_cells", default=None)
    ap.add_argument("--markers", default=None, help="marker_expression.tsv (SC_MARKER_EXPRESSION)")
    ap.add_argument("--samplesheet", default=None)
    ap.add_argument("--patient_col", default="patient")
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    written = []

    def write(df, name, unit):
        if df is None or not len(df):
            log(f"{name}: no input on this route, not written")
            return
        df.to_csv(os.path.join(a.outdir, name), sep="\t", index=False)
        written.append((name, unit, f"{len(df):,} rows"))
        log(f"wrote {name} ({len(df):,} rows x {df.shape[1]} columns)")

    clones, smeta = build_clones(a)
    write(build_samples(a, smeta), "samples.tsv", "sample")
    write(clones, "clones.tsv.gz", "sample x TCRbeta clonotype")
    write(build_tests(a), "clone_tests.tsv.gz", "clonotype x timepoint pair x test")
    write(read(find(a.indir, "overlap_samples.tsv")), "sample_overlap.tsv", "sample pair")
    write(build_cells(a), "cells.tsv.gz", "cell")
    write(read(a.markers), "marker_expression.tsv", "gene x cell state x clone group")

    inv = "\n".join(f"| `{n}` | {u} | {r} |" for n, u, r in written) or "| (none) | | |"
    with open(os.path.join(a.outdir, "README.md"), "w") as fh:
        fh.write(README.format(inventory=inv))


if __name__ == "__main__":
    main()
