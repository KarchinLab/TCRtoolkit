#!/usr/bin/env python3
"""Rebuild anchors_human_trb.tsv.gz, the reference sample SELECTION_TEST calibrates against.

Anchors are productive CDR3beta sequences drawn from the OLGA default human TRB model, with
their P_gen. They depend only on the model, so they ship with the pipeline instead of being
recomputed per run (about two core-hours). Usage: make_anchors.py --out FILE [--n N --seed S]
"""
import argparse
import multiprocessing as mp
import os

import numpy as np
import olga
import olga.generation_probability as gp
import olga.load_model as lm
import olga.sequence_generation as sg
import pandas as pd

MODEL = os.path.join(os.path.dirname(olga.__file__), "default_models", "human_T_beta") + os.sep


def load():
    genomic = lm.GenomicDataVDJ()
    genomic.load_igor_genomic_data(MODEL + "model_params.txt", MODEL + "V_gene_CDR3_anchors.csv",
                                   MODEL + "J_gene_CDR3_anchors.csv")
    model = lm.GenerativeModelVDJ()
    model.load_and_process_igor_model(MODEL + "model_marginals.txt")
    return genomic, model


def work(job):
    n, seed = job
    genomic, model = load()
    np.random.seed(seed)
    gen = sg.SequenceGenerationVDJ(model, genomic)
    pgen = gp.GenerationProbabilityVDJ(model, genomic)
    seqs = [gen.gen_rnd_prod_CDR3()[1] for _ in range(n)]
    return [(s, pgen.compute_aa_CDR3_pgen(s, print_warnings=False)) for s in seqs]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=300000)
    ap.add_argument("--seed", type=int, default=20261104)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    jobs = [(a.n // a.workers + 1, a.seed + i) for i in range(a.workers)]
    with mp.Pool(a.workers) as pool:
        rows = [r for part in pool.map(work, jobs) for r in part]
    pd.DataFrame(rows, columns=["cdr3", "pgen"]).to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
