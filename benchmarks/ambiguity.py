#!/usr/bin/env python3
"""Ambiguity groups benchmark (Phase 14, F57, D39).

Runs NanoVI's RUN_VI and WRITE_OUTPUT scripts (bin/estimate_abundance.py --groups,
bin/write_abundance_table.py --groups-output; NANOVI_REPO at a NanoVI commit with
them, which the report names) on the score files Phase 13 stored for the exact
database (the Phase 7a merged build, cross-genus rule keep, D38; minimap2 2.31 +
samtools 1.24, -N 3): the 92 simulations (phase13/runs/<sample>/lik/exact.npz) and
the 13 real sup mocks (phase13/real/runs/<sample>/lik/exact.npz). Each sample is
run at every --group-min-shared in TAUS with VB (and, with --credible-intervals,
the interval columns); the 10k simulations also with EM at the default. The
identical-16S relatives come from the merged build's identical_16s.tsv.

Scores, per sample, of the species table and of the group table against the truth
summed over the same groups (evaluate.compare on the simulations' realised shares,
evaluate.compare_real on the real mocks' known shares):
  strict   a group holds its reported members only
  lenient  a group also holds its identical_16s_relatives (unreported species that
           share an identical 16S with a member); a relative named by two groups
           joins them
plus the abundance and reads in groups of two or more rows, the species the species
table had right (within RIGHT_TOLERANCE of the truth) that sit in such a group, the
groups that are right themselves, groups across genera, and at 10k reads the
coverage of the group intervals (design abundance renormalised over the reported
species, as intervals.py scores the species).

Writes <phase14>/runs/<sample>/<arm>/ (est.json, species.tsv, groups.tsv) and
<phase14>/results/ambiguity.tsv (one row per sample, arm and table),
ambiguity_groups.tsv (one row per group of two or more rows) and
ambiguity_coverage.tsv.
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from common import BENCH, BIN, NANOVI, PHASE13, PHASE14, write_tsv
from evaluate import compare, compare_real, load_real_truth, load_truth

PYTHON = os.environ.get("BENCH_PYTHON", sys.executable)
DB = BENCH / "db" / "exact"
IDENTICAL = PHASE13 / "db" / "merged_keep" / "identical_16s.tsv"
TAUS = [0.2, 0.25, 0.3, 0.4]
DEFAULT_TAU = 0.25
RIGHT_TOLERANCE = 0.25
PLACEHOLDERS = {"Below detection", "Unassigned"}


def samples():
    sims = sorted(p.parent.parent.name for p in (PHASE13 / "runs").glob("*_s*/lik/exact.npz"))
    real = sorted(p.parent.parent.name for p in (PHASE13 / "real" / "runs").glob("*/lik/exact.npz"))
    return [(s, "sim") for s in sims] + [(s, "real") for s in real]


def lik_dir(sample, kind):
    return (PHASE13 / "runs" / sample if kind == "sim" else PHASE13 / "real" / "runs" / sample) / "lik"


def arms(sample):
    out = [(f"vb_tau{t}", "vb", t) for t in TAUS]
    if "_10k_" in sample:
        out.append((f"em_tau{DEFAULT_TAU}", "em", DEFAULT_TAU))
    return out


def run_arm(sample, kind, arm, estimator, tau):
    out = PHASE14 / "runs" / sample / arm
    if (out / "groups.tsv").exists():
        return
    out.mkdir(parents=True, exist_ok=True)
    lik = lik_dir(sample, kind)
    taxonomy = DB / "taxonomy.tsv"
    subprocess.run([PYTHON, str(BIN / "estimate_abundance.py"), "--scores", str(lik / "exact.npz"),
                    "--output", str(out / "tmp_est.json"), "--estimator", estimator, "--taxonomy", str(taxonomy),
                    "--groups", "--group-min-shared", str(tau)],
                   check=True, stdout=subprocess.DEVNULL, stderr=open(out / "est.log", "w"))
    (out / "tmp_est.json").rename(out / "est.json")
    subprocess.run([PYTHON, str(BIN / "write_abundance_table.py"), "--abundance", str(out / "est.json"),
                    "--read-counts", str(lik / "exact.counts.json"), "--taxonomy", str(taxonomy),
                    "--output-prefix", str(out / "species"), "--counts", "--identical-16s", str(IDENTICAL),
                    "--groups-output", str(out / "tmp_groups")]
                   + (["--credible-intervals"] if estimator == "vb" else []),
                   check=True, stdout=subprocess.DEVNULL, stderr=open(out / "write.log", "w"))
    (out / "tmp_groups.tsv").rename(out / "groups.tsv")


def _job(args):
    run_arm(*args)
    return args[0], args[2]


class Union:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def join(self, a, b):
        self.parent[self.find(a)] = self.find(b)


def group_keys(groups, lenient):
    """species name -> group key: the group of a reported member, and (lenient) of its relatives."""
    u = Union()
    for g in groups.itertuples():
        for m in g.members.split("; "):
            u.join(m, g.group)
        if lenient:
            for r in filter(None, str(g.identical_16s_relatives).split("; ")):
                u.join(r, g.group)
    return {k: u.find(k) for k in u.parent if not k.startswith("G")}


def collapse(abund, keys):
    out = {}
    for name, v in abund.items():
        k = keys.get(name, name)
        out[k] = out.get(k, 0.0) + v
    return out


def collapse_truth(truth, keys):
    """Real-mock truth: known shares summed; a group of presence-only species stays NaN."""
    out = {}
    for name, v in truth.items():
        k = keys.get(name, name)
        if v == v:
            current = out.get(k)
            out[k] = v if current is None or current != current else current + v
        else:
            out.setdefault(k, float("nan"))
    return out


def right_species(est, truth, kind):
    """Truth species the species table had within RIGHT_TOLERANCE of their share."""
    if kind == "real":
        known = {k: v for k, v in truth.items() if v == v}
        s = sum(est.get(k, 0.0) for k in known)
        tk = sum(known.values())
        return {k for k, v in known.items() if s and abs(est.get(k, 0.0) / s - v / tk) <= RIGHT_TOLERANCE * v / tk}
    return {k for k, v in truth.items() if v > 0 and abs(est.get(k, 0.0) - v) <= RIGHT_TOLERANCE * v}


def score_arm(sample, kind, arm):
    out = PHASE14 / "runs" / sample / arm
    read = dict(sep="\t", keep_default_na=False)
    species = pd.read_csv(out / "species.tsv", **read)
    groups = pd.read_csv(out / "groups.tsv", **read)
    species = species[~species["species"].isin(PLACEHOLDERS)]
    groups = groups[~groups["group"].isin(PLACEHOLDERS)].copy()
    groups["n_members"] = groups["n_members"].astype(int)
    est = dict(zip(species["species"], species["abundance"].astype(float)))
    if kind == "sim":
        t = load_truth(BENCH / "sims" / sample)
        truth = dict(zip(t["species"], t["realised"]))
    else:
        truth = load_real_truth(sample)
    n_reads = json.loads((lik_dir(sample, kind) / "exact.counts.json").read_text())["assigned_count"]
    assigned = species["estimated_counts"].astype(int).sum()
    multi = groups[groups["n_members"] > 1]
    in_groups = {m for ms in multi["members"] for m in ms.split("; ")}
    right = right_species(est, truth, kind)
    base = {"sample": sample, "kind": kind, "arm": arm,
            "groups": len(multi), "abundance_in_groups": multi["abundance"].astype(float).sum(),
            "reads_in_groups": multi["estimated_counts"].astype(int).sum() / assigned if assigned else 0.0,
            "species_right": len(right), "species_right_hidden": len(right & in_groups),
            "groups_cross_genus": int(sum(len({m.split()[0] for m in ms.split("; ")}) > 1 for ms in multi["members"]))}
    rows, group_rows = [], []
    for table in ("species", "strict", "lenient"):
        keys = {} if table == "species" else group_keys(groups, lenient=table == "lenient")
        e = collapse(est, keys)
        if kind == "sim":
            tr = collapse(truth, keys)
            m = compare(e, tr, "species")
        else:
            tr = collapse_truth(truth, keys)
            m = compare_real(e, tr, n_reads, "species")
        rows.append({**base, "table": table, **{k.replace("species_", ""): v for k, v in m.items()}})
        if table != "species":
            for g in multi.itertuples():
                k = keys[g.members.split("; ")[0]]
                tg = tr.get(k, 0.0)
                group_rows.append({"sample": sample, "kind": kind, "arm": arm, "scoring": table, "group": g.group,
                                   "members": g.members, "n_members": g.n_members, "abundance": float(g.abundance),
                                   "truth": tg, "lowest_common_rank": g.lowest_common_rank,
                                   "max_shared": float(g.max_shared),
                                   "relatives": g.identical_16s_relatives,
                                   "right": tg == tg and tg > 0 and abs(float(g.abundance) - tg) <= RIGHT_TOLERANCE * tg,
                                   "members_right": len(right & set(g.members.split("; ")))})
    return rows, group_rows


def coverage(sample, arm):
    """At 10k reads: does each interval contain the design abundance renormalised over the reported species?"""
    out = PHASE14 / "runs" / sample / arm
    read = dict(sep="\t", keep_default_na=False)
    species = pd.read_csv(out / "species.tsv", **read)
    groups = pd.read_csv(out / "groups.tsv", **read)
    species = species[~species["species"].isin(PLACEHOLDERS)]
    groups = groups[~groups["group"].isin(PLACEHOLDERS)]
    t = load_truth(BENCH / "sims" / sample)
    design = dict(zip(t["species"], t["design_abundance"]))
    total = sum(design.get(s, 0.0) for s in species["species"])
    rows = []
    for s in species.itertuples():
        if s.species in design:
            truth = design[s.species] / total
            rows.append({"sample": sample, "arm": arm, "level": "species", "name": s.species, "truth": truth,
                         "covers": float(s.ci_lower) <= truth <= float(s.ci_upper)})
    for g in groups.itertuples():
        members = g.members.split("; ")
        if int(g.n_members) > 1 and all(m in design for m in members):
            truth = sum(design[m] for m in members) / total
            rows.append({"sample": sample, "arm": arm, "level": "group", "name": g.members, "truth": truth,
                         "covers": float(g.ci_lower) <= truth <= float(g.ci_upper)})
            rows += [{**r, "level": "group_member"} for r in rows if r["level"] == "species" and r["name"] in members]
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--score-only", action="store_true")
    parser.add_argument("--only", help="run the samples whose name contains this (a smoke test; no results written)")
    args = parser.parse_args()
    print(f"NanoVI: {NANOVI}", file=sys.stderr)
    todo = [(s, k, a, e, t) for s, k in samples() for a, e, t in arms(s) if not args.only or args.only in s]
    if not args.score_only:
        with ProcessPoolExecutor(args.jobs) as pool:
            for i, (s, a) in enumerate(pool.map(_job, todo), 1):
                print(f"[{i}/{len(todo)}] {s} {a}", file=sys.stderr, flush=True)
    rows, group_rows, cov = [], [], []
    for s, k, a, e, t in todo:
        r, g = score_arm(s, k, a)
        rows += r
        group_rows += g
        if k == "sim" and "_10k_" in s and e == "vb":
            cov += coverage(s, a)
    if args.only:
        d = pd.DataFrame(rows)
        print(d.groupby(["arm", "table"])[[c for c in ("l1", "l1_known", "groups") if c in d]].mean())
        return
    write_tsv(pd.DataFrame(rows), PHASE14 / "results" / "ambiguity.tsv")
    write_tsv(pd.DataFrame(group_rows), PHASE14 / "results" / "ambiguity_groups.tsv")
    write_tsv(pd.DataFrame(cov), PHASE14 / "results" / "ambiguity_coverage.tsv")


if __name__ == "__main__":
    main()
