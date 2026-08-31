#!/usr/bin/env python3
"""
CLUSTER_ROLLUP

Summarises the clonotype-clustering modules (GIANA / GLIPH2 / TCRdist3) into the
rollup tables the Master Summary expects. Those three modules write raw per-patient
/ per-sample output but no rollup, so MASTER_SUMMARY's giana_summary_file /
gliph2_summary_file / tcrdist3_summary_file had nothing to read.

Outputs (written into --outdir, matching the .qmd's expected filenames):

    giana_summary_rollup.tsv        metric / value
    gliph2_summary_rollup.tsv       metric / value
    tcrdist3_summary_rollup.tsv     metric / value

    method_presence_summary.tsv     method / frac_cells      (fallback: see below)
    method_cluster_counts.tsv       method / n_clusters      (fallback: see below)

    annotation_giana_summary.tsv    annot / frac_giana       (only when informative)
    annotation_gliph2_summary.tsv   annot / frac_gliph
    annotation_tcrdist3_summary.tsv annot / frac_tcrdist

The two method_* tables are also produced by CONSENSUS_CLUSTERING, which sees CoNGA
and consensus labels too and is therefore the better source. The workflow prefers the
CONSENSUS copies when consensus ran and falls back to these otherwise (VDJ-only).

A file is simply not written when its inputs are absent, so the Master Summary's
0-byte/NO_FILE handling treats it as missing rather than empty.
"""

import argparse
import os
import sys

import pandas as pd

# CDR3 amino-acid column, in resolution order, across the various export flavours.
CDR3_COLS = ["cdr3b", "CDR3b", "junction_aa", "cdr3_b_aa", "cdr3"]
ANNOT_COLS = ["annot", "Annotation", "celltype", "CellType", "seurat_clusters"]


def usable(path):
    """A staged NO_FILE placeholder is 0 bytes; treat it (and absent paths) as missing."""
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def tag_from_name(path, suffix):
    """Patient/sample identity is carried in the filename prefix, e.g. HRS371754_giana.txt."""
    base = os.path.basename(path)
    return base[: -len(suffix)] if base.endswith(suffix) else os.path.splitext(base)[0]


def fmt(v):
    """Keep the metric/value column as text so integer counts don't render as '263.0'."""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def write_rollup(outdir, name, rows):
    df = pd.DataFrame([(m, fmt(v)) for m, v in rows], columns=["metric", "value"])
    df.to_csv(os.path.join(outdir, name), sep="\t", index=False)
    print(f"[cluster_rollup] wrote {name} ({len(df)} metrics)")


def read_giana(paths):
    """GIANA writes two '##' comment lines before the real header."""
    frames = []
    for p in paths:
        if not usable(p):
            continue
        try:
            df = pd.read_csv(p, sep="\t", comment="#", low_memory=False)
        except Exception as exc:  # a truncated/empty per-patient file must not sink the run
            print(f"[cluster_rollup] WARN: could not read {p}: {exc}", file=sys.stderr)
            continue
        if df.empty:
            continue
        df["__group"] = tag_from_name(p, "_giana.txt")
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else None


def read_gliph2(paths):
    frames = []
    for p in paths:
        if not usable(p):
            continue
        try:
            df = pd.read_csv(p, sep="\t", low_memory=False)
        except Exception as exc:
            print(f"[cluster_rollup] WARN: could not read {p}: {exc}", file=sys.stderr)
            continue
        if df.empty:
            continue
        df["__group"] = tag_from_name(p, "_cluster_member_details.txt")
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else None


def read_tcrdist(paths):
    # SAMPLE emits clone_df.csv alongside the .hdf5/.csv distance matrices; only the
    # clone tables are parseable here, so ignore everything else rather than warning on it.
    frames = []
    for p in paths:
        if not usable(p) or not p.endswith("_clone_df.csv"):
            continue
        try:
            df = pd.read_csv(p, low_memory=False)
        except Exception as exc:
            print(f"[cluster_rollup] WARN: could not read {p}: {exc}", file=sys.stderr)
            continue
        if df.empty:
            continue
        df["__group"] = tag_from_name(p, "_clone_df.csv")
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else None


def tcrdist_clusters(matrix_paths, radius):
    """
    Threshold each sample's pairwise distance matrix at `radius` and count connected
    components of size >= 2 - the same definition enrich_seurat.R uses.

    Only the Full-SC route ran that clustering (it lives in enrich_seurat.R), so on the
    VDJ-only route TCRdist3 produced distance matrices that were never thresholded and the
    method reported 0 clusters. Doing it here makes TCRdist3 comparable to GIANA/GLIPH2 on
    both routes.

    Returns (n_clusters_total, {sample: n_clusters}, n_clustered_clones).
    """
    try:
        import numpy as np
        import scipy.sparse as sp
        from scipy.sparse.csgraph import connected_components
    except ImportError as exc:
        print(f"[cluster_rollup] WARN: scipy/numpy unavailable, skipping tcrdist clustering: {exc}",
              file=sys.stderr)
        return 0, {}, 0

    per_sample = {}
    total_clustered = 0

    for p in matrix_paths:
        if not usable(p):
            continue
        name = os.path.basename(p)
        sample = name.split("_distance_matrix")[0]
        ext = os.path.splitext(p)[1].lower()
        try:
            if ext == ".csv":
                full = np.loadtxt(p, delimiter=",")
                adj = sp.csr_matrix((full <= radius) & (full > 0))
            elif ext == ".hdf5":
                import h5py
                with h5py.File(p, "r") as f:
                    m = sp.csr_matrix(
                        (f["data"][:], f["indices"][:], f["indptr"][:]),
                        shape=tuple(f["shape"][:]),
                    )
                # Stored sparse: an absent entry means "not within the sparsity cutoff",
                # so only explicit non-zero entries can be within radius.
                coo = m.tocoo()
                keep = (coo.data <= radius) & (coo.data > 0)
                adj = sp.coo_matrix(
                    (np.ones(keep.sum()), (coo.row[keep], coo.col[keep])), shape=m.shape
                ).tocsr()
            else:
                continue
        except Exception as exc:
            print(f"[cluster_rollup] WARN: could not read {p}: {exc}", file=sys.stderr)
            continue

        n_comp, labels = connected_components(adj, directed=False)
        sizes = pd.Series(labels).value_counts()
        n_real = int((sizes >= 2).sum())
        per_sample[sample] = n_real
        total_clustered += int(sizes[sizes >= 2].sum())

    return sum(per_sample.values()), per_sample, total_clustered


def first_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def resolve_cdr3(df):
    """Return a Series of CDR3b strings, falling back to the '<CDR3>_<TRBV>' clone_id form."""
    col = first_col(df, CDR3_COLS)
    if col is not None:
        return df[col].astype("string").str.strip()
    for key in ("clone_id", "CTaa", "clonotype"):
        if key in df.columns:
            return df[key].astype("string").str.split("_", n=1).str[0].str.strip()
    return None


def cluster_key(df, cluster_col):
    """
    GIANA and GLIPH2 both number/label clusters *within* a patient, so the same id
    reappears across patients for entirely unrelated clusters (45 of 77 collide in the
    8-sample test set). Counting on the raw id silently merges them - scope it by the
    patient the file came from.
    """
    if cluster_col not in df.columns:
        return None
    return df["__group"].astype("string") + "::" + df[cluster_col].astype("string")


def cluster_series(df, cdr3, cluster_col):
    """Map CDR3 -> set of CDR3s that the method actually assigned to a cluster."""
    if cluster_col not in df.columns or cdr3 is None:
        return set()
    keep = df[cluster_col].notna() & (df[cluster_col].astype("string").str.strip() != "")
    return set(cdr3[keep].dropna().unique())


def annotation_table(export, annot_col, cdr3_export, clustered, frac_name):
    """Per-annotation fraction of cells whose clonotype was clustered by this method."""
    tmp = pd.DataFrame({"annot": export[annot_col].astype("string"), "cdr3": cdr3_export})
    tmp = tmp.dropna(subset=["annot"])
    if tmp.empty:
        return None
    tmp["hit"] = tmp["cdr3"].isin(clustered)
    out = (
        tmp.groupby("annot", dropna=True)["hit"]
        .mean()
        .reset_index()
        .rename(columns={"hit": frac_name})
        .sort_values(frac_name, ascending=False)
    )
    return out


def main():
    ap = argparse.ArgumentParser(description="Roll up GIANA / GLIPH2 / TCRdist3 outputs")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--giana", nargs="*", default=[])
    ap.add_argument("--gliph2", nargs="*", default=[])
    ap.add_argument("--tcrdist", nargs="*", default=[])
    ap.add_argument("--tcrdist-matrix", nargs="*", default=[],
                    help="per-sample distance matrices (.hdf5/.csv) to threshold at --tcrdist-radius")
    ap.add_argument("--tcrdist-radius", type=float, default=24.0)
    ap.add_argument("--export-cells", default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    giana = read_giana(args.giana)
    gliph = read_gliph2(args.gliph2)
    tcrd = read_tcrdist(args.tcrdist)

    export = None
    cdr3_export = None
    annot_col = None
    if usable(args.export_cells):
        try:
            export = pd.read_csv(args.export_cells, sep="\t", low_memory=False)
            cdr3_export = resolve_cdr3(export)
            annot_col = first_col(export, ANNOT_COLS)
        except Exception as exc:
            print(f"[cluster_rollup] WARN: could not read export_cells: {exc}", file=sys.stderr)
            export = None

    clustered = {}
    counts = {}

    # ── GIANA ────────────────────────────────────────────────────────────────
    if giana is not None:
        cdr3 = resolve_cdr3(giana)
        gkey = cluster_key(giana, "cluster")
        n_clusters = gkey.nunique() if gkey is not None else 0
        sizes = gkey.value_counts() if gkey is not None else pd.Series(dtype=int)
        write_rollup(args.outdir, "giana_summary_rollup.tsv", [
            ("Patients with clusters", giana["__group"].nunique()),
            ("Clonotypes clustered", int(len(giana))),
            ("Unique CDR3b clustered", int(cdr3.nunique()) if cdr3 is not None else 0),
            ("Clusters detected", int(n_clusters)),
            ("Largest cluster size", int(sizes.max()) if len(sizes) else 0),
            ("Median cluster size", float(sizes.median()) if len(sizes) else 0.0),
            ("Samples represented", int(giana["sample"].nunique()) if "sample" in giana.columns else 0),
        ])
        clustered["GIANA"] = cluster_series(giana, cdr3, "cluster")
        counts["GIANA"] = int(n_clusters)

    # ── GLIPH2 ───────────────────────────────────────────────────────────────
    if gliph is not None:
        cdr3 = resolve_cdr3(gliph)
        pkey = cluster_key(gliph, "tag")
        n_clusters = pkey.nunique() if pkey is not None else 0
        sizes = pkey.value_counts() if pkey is not None else pd.Series(dtype=int)
        write_rollup(args.outdir, "gliph2_summary_rollup.tsv", [
            ("Patients with clusters", gliph["__group"].nunique()),
            ("Clonotype-cluster memberships", int(len(gliph))),
            ("Unique CDR3b clustered", int(cdr3.nunique()) if cdr3 is not None else 0),
            ("Motif clusters detected", int(n_clusters)),
            ("Largest cluster size", int(sizes.max()) if len(sizes) else 0),
            ("Median cluster size", float(sizes.median()) if len(sizes) else 0.0),
            ("Samples represented", int(gliph["sample"].nunique()) if "sample" in gliph.columns else 0),
        ])
        clustered["GLIPH2"] = cluster_series(gliph, cdr3, "tag")
        counts["GLIPH2"] = int(n_clusters)

    # ── TCRdist3 ─────────────────────────────────────────────────────────────
    # clone_df carries no cluster labels (those come from thresholding the distance
    # matrix in enrich_seurat.R), so cluster counts are taken from the per-cell export
    # when CLUSTER_TO_SC has already written a tcrdist_cluster column.
    if tcrd is not None:
        cdr3 = resolve_cdr3(tcrd)

        # Prefer the per-cell tcrdist_cluster column when CLUSTER_TO_SC already wrote it
        # (Full SC); otherwise threshold the distance matrices here so the VDJ-only route
        # reports real cluster counts instead of 0.
        n_tcrdist_clusters = 0
        per_sample = {}
        n_clustered_clones = 0
        if export is not None and "tcrdist_cluster" in export.columns:
            n_tcrdist_clusters = int(export["tcrdist_cluster"].nunique(dropna=True))
        elif args.tcrdist_matrix:
            n_tcrdist_clusters, per_sample, n_clustered_clones = tcrdist_clusters(
                args.tcrdist_matrix, args.tcrdist_radius
            )

        rows = [
            ("Samples analysed", tcrd["__group"].nunique()),
            ("Clones in distance matrices", int(len(tcrd))),
            ("Unique CDR3b", int(cdr3.nunique()) if cdr3 is not None else 0),
            ("Median clones per sample", float(tcrd.groupby("__group").size().median())),
            ("Max clones per sample", int(tcrd.groupby("__group").size().max())),
            (f"Clusters detected (radius {args.tcrdist_radius:g})", n_tcrdist_clusters),
        ]
        if n_clustered_clones:
            rows.append(("Clones in a cluster", n_clustered_clones))
        if per_sample:
            rows.append(("Median clusters per sample",
                         float(pd.Series(list(per_sample.values())).median())))
        write_rollup(args.outdir, "tcrdist3_summary_rollup.tsv", rows)

        if cdr3 is not None:
            clustered["TCRdist3"] = set(cdr3.dropna().unique())
        if n_tcrdist_clusters:
            counts["TCRdist3"] = n_tcrdist_clusters

    # ── method_* fallbacks (CONSENSUS supersedes these when it runs) ──────────
    if counts:
        pd.DataFrame(
            sorted(counts.items()), columns=["method", "n_clusters"]
        ).to_csv(os.path.join(args.outdir, "method_cluster_counts.tsv"), sep="\t", index=False)
        print(f"[cluster_rollup] wrote method_cluster_counts.tsv ({len(counts)} methods)")

    if clustered and export is not None and cdr3_export is not None:
        n_cells = len(export)
        rows = [
            (m, float(cdr3_export.isin(s).sum()) / n_cells if n_cells else 0.0)
            for m, s in sorted(clustered.items())
        ]
        pd.DataFrame(rows, columns=["method", "frac_cells"]).to_csv(
            os.path.join(args.outdir, "method_presence_summary.tsv"), sep="\t", index=False
        )
        print(f"[cluster_rollup] wrote method_presence_summary.tsv ({len(rows)} methods)")

    # ── per-annotation coverage ──────────────────────────────────────────────
    # Only meaningful when the export actually carries >1 distinct annotation;
    # in VDJ-only mode every cell is "Unannotated" and the panel would be a single bar.
    if export is not None and annot_col is not None and cdr3_export is not None:
        if export[annot_col].nunique(dropna=True) > 1:
            for method, frac_name, fname in (
                ("GIANA", "frac_giana", "annotation_giana_summary.tsv"),
                ("GLIPH2", "frac_gliph", "annotation_gliph2_summary.tsv"),
                ("TCRdist3", "frac_tcrdist", "annotation_tcrdist3_summary.tsv"),
            ):
                if method not in clustered:
                    continue
                tbl = annotation_table(export, annot_col, cdr3_export, clustered[method], frac_name)
                if tbl is not None and not tbl.empty:
                    tbl.to_csv(os.path.join(args.outdir, fname), sep="\t", index=False)
                    print(f"[cluster_rollup] wrote {fname} ({len(tbl)} annotations)")
        else:
            print("[cluster_rollup] single annotation value in export; skipping annotation panels")

    print("[cluster_rollup] done")


if __name__ == "__main__":
    main()
