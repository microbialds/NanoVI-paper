#!/usr/bin/env python3
"""Run NanoVI's estimate_abundance.py with ambiguity-group links limited to one family.

A test of a change to D39 that NanoVI does not have: two table rows link only when
their species are in the same family of the database's taxonomy (the --taxonomy
argument, taken as estimate_abundance.py takes it). The rows' f is computed as in
NanoVI (bin/ambiguity_groups.py), and the f between rows of different families is set
to 0 before the rows are linked; families are transitive, so every group stays inside
one family. A row's family is that of its first tax ID (a row's tax IDs share a species
name, F46). The arguments are estimate_abundance.py's, passed through unchanged.
"""
import runpy
import sys

import numpy as np
import pandas as pd

from common import BIN

import ambiguity_groups as ag  # noqa: E402  (bin/, put on sys.path by common)


def _argument(name):
    return sys.argv[sys.argv.index(name) + 1]


def main():
    family = pd.read_csv(_argument("--taxonomy"), sep="\t", dtype=str).set_index("tax_id")["family"].to_dict()
    nanovi_groups = ag.ambiguity_groups
    nanovi_pairwise = ag.pairwise_shared

    def family_groups(idx, phi, reported_ids, row_keys=None, min_shared=ag.GROUP_MIN_SHARED):
        row_of_column, n_rows = ag.table_rows(reported_ids, row_keys)
        row_family = [None] * n_rows
        for sid, r in zip(reported_ids, row_of_column):
            if row_family[r] is None:
                row_family[r] = family[str(sid)]
        row_family = np.array(row_family, dtype=object)

        def pairwise_within_family(*args, **kwargs):
            expected, shared = nanovi_pairwise(*args, **kwargs)
            shared[row_family[:, None] != row_family[None, :]] = 0.0
            return expected, shared

        ag.pairwise_shared = pairwise_within_family
        try:
            return nanovi_groups(idx, phi, reported_ids, row_keys, min_shared)
        finally:
            ag.pairwise_shared = nanovi_pairwise

    ag.ambiguity_groups = family_groups
    script = str(BIN / "estimate_abundance.py")
    sys.argv = [script] + sys.argv[1:]
    runpy.run_path(script, run_name="__main__")


if __name__ == "__main__":
    main()
