# Ambiguity groups (Phase 14, F57, November 2026)

When two species' 16S genes are identical or nearly so, their reads fit both, and the
estimator splits those reads by abundance: one member of the pair takes them. NanoVI's
species table then reports each species alone, with a split the reads do not support.
`--ambiguity_groups` (decision D39) writes a second per-sample table,
`<sample>_groups.tsv`, in which such species are reported together with their sum. This
benchmark asks how much closer to the truth that table is, what it costs (correct species
merged into a group, groups that make no biological sense) and which threshold to use.
The species table never changes. The option stays off by default.

The design and the options considered are in `docs/design/ambiguity_groups.md`. Every
table cited here is in `docs/benchmarks/data/`, and the figure is in
`docs/benchmarks/figures/` with its source data under the same name in `data/`.
Per-sample outputs are under `~/pipelines/NanoVI_benchmarks/phase14` (D11) and are not in
the repository.

## 1. The rule

NanoVI branch `roadmap/p14-ambiguity-groups`, commit `9137b7b` (`bin/ambiguity_groups.py`).
RUN_VI (`estimate_abundance.py --groups`) works on the rows of the species table (tax IDs
merged as the table merges them, F46), using the responsibilities φ of the final fit (VB,
or EM's posterior):

  f_a(b) = Σ_r φ_ra · [read r's share on b ≥ 0.01] / E_a,   E_a = Σ_r φ_ra

f_a(b) is the share of a's expected reads that sit on reads b also takes at least 1% of.
Rows a and b are linked when max(f_a(b), f_b(a)) ≥ τ. In sections 3-6 there is no other
condition. D39 was amended after the benchmark (2026-11-05, section 7): rows link only
within one family of the database's taxonomy, and a row with a blank family links to
nothing. That is the rule NanoVI ships (commit `3b24374`). Groups are the connected
components of the links, so every species row is in exactly one group, a row that links
to nothing being a group of one. τ = 0.25 is fixed in the pipeline;
`--group-min-shared` exists in the script for this benchmark only. WRITE_OUTPUT
(`write_abundance_table.py --groups-output`) writes each group with:
- its members and their summed abundance and counts;
- its lowest common rank;
- `max_shared` (the largest f between two members);
- its own `shared_read_fraction` (reads split with rows outside the group);
- with `--credible_intervals`, an interval from the summed α;
- `identical_16s_relatives`: the species that share an identical 16S with a member in the
  database (`identical_16s_species`, D24) and are not reported in the sample.

## 2. How it ran

`benchmarks/ambiguity.py` (NANOVI_REPO at `9137b7b`) runs the two scripts on the score
files Phase 13 stored for the exact database, so no read was realigned:
- **Database:** the Phase 7a merged GTDB r226 build with `--cross-genus-rule keep` (D38),
  minimap2 2.31 + samtools 1.24, `-N 3`, fastplong 0.7.0.
- **Simulations (92):** even and log-staggered, R9.4.1-like and R10.4.1-like, 1k and 10k
  reads at seeds 1-10, 50k at seeds 1-3: `phase13/runs/<sample>/lik/exact.npz`.
- **Real mocks (13 sup samples, 10,000 reads each):** `phase13/real/runs/<sample>/lik/exact.npz`.
- **Arms:** τ ∈ {0.2, 0.25, 0.3, 0.4} with VB and `--credible-intervals`, and EM at
  τ = 0.25 on the 40 samples at 10k.
- **Identical-16S relatives:** from the merged build's `identical_16s.tsv` (4,335 species).

460 runs, no failures. After the code review the branch computes a group's shared reads on the sparse
reads × rows matrix instead of through `posterior_read_summary`, and a group of one row
copies its species row's values. Rerun with the reviewed code on three samples
(emu_gutmock_SRR14326315, zymo_MSPlus_sup_rep3, even_r941_10k_s01), the group tables have
the same groups, members, abundances and counts. `shared_read_fraction` and the
intervals differ by at most 1.1e-16. No number below depends on it. `analyse.py ambiguity` writes the tables (`ambiguity_*.tsv`) and
`fig_ambiguity.png`.

Scoring (`evaluate.compare` on the simulations' realised read shares, `compare_real` on
the real mocks' known shares) puts the truth over the same groups as the estimate. There
are two ways to do that:

- **strict**: a group holds its reported members only. A truth species that is not
  reported stays on its own.
- **lenient**: a group also holds its `identical_16s_relatives`, so their truth counts
  towards it (a relative named by two groups joins them). This is the reading the table
  invites: "*S. aureus*, which may include *S. argenteus*".

A species is **right** in the species table when its estimate is within 25% of its truth
(of its known share renormalised over the known species, on the real mocks). A right
species **hidden** by a group is one that ends up in a group of two or more rows.

Groups of two or more rows are counted under "groups" below. "Reads in groups" is their
`estimated_counts` over the assigned reads.

## 3. Simulations

Means over the 92 simulations (`ambiguity_sims_overall.tsv`; paired counts from
`ambiguity_sims_paired.tsv`, group table against the species table of the same sample):

| τ | species L1 | strict L1 | lenient L1 | groups | reads in groups | right species hidden (of 12.1) | strict L1 lower / higher than species |
|---|---|---|---|---|---|---|---|
| 0.20 | 0.132 | 0.068 | 0.0095 | 2.36 | 0.34 | 1.83 | 88 / 0 |
| **0.25** | 0.132 | **0.075** | **0.0096** | 2.08 | 0.29 | 1.54 | 83 / 0 |
| 0.30 | 0.132 | 0.091 | 0.011 | 1.61 | 0.21 | 1.07 | 75 / 0 |
| 0.40 | 0.132 | 0.116 | 0.018 | 0.61 | 0.06 | 0.25 | 38 / 4 |

The lenient score is lower than the species score in all 92 samples at every τ. With
no in-sample groups at all, the relatives alone take it to 0.021 (preliminary,
`docs/design/ambiguity_groups.md`). The strict score is never higher than the species
score at τ ≤ 0.3.

By set, at τ = 0.25 (`ambiguity_sims_summary.tsv`):

| Set | species L1 | strict L1 | lenient L1 | groups | reads in groups | right species | hidden |
|---|---|---|---|---|---|---|---|
| even R10.4.1 1k / 10k / 50k | 0.183 / 0.185 / 0.188 | 0.111 / 0.109 / 0.108 | 0.005 / 0.005 / 0.005 | 2.9 / 2.9 / 3.0 | 0.28 / 0.29 / 0.30 | 16 | 3.3 / 3.8 / 4.0 |
| even R9.4.1 1k / 10k / 50k | 0.306 / 0.275 / 0.278 | 0.193 / 0.129 / 0.123 | 0.020 / 0.008 / 0.007 | 3.1 / 2.9 / 3.0 | 0.25 / 0.28 / 0.30 | 12 | 0.5 / 0 / 0 |
| log R10.4.1 1k / 10k / 50k | 0.023 / 0.030 / 0.030 | 0.016 / 0.009 / 0.010 | 0.002 / 0.005 / 0.005 | 0.8 / 1.1 / 1.0 | 0.17 / 0.41 / 0.40 | 10-11 | 0.8 / 1.1 / 1.0 |
| log R9.4.1 1k / 10k / 50k | 0.030 / 0.026 / 0.028 | 0.027 / 0.019 / 0.017 | 0.019 / 0.015 / 0.014 | 1.1 / 1.6 / 2.0 | 0.23 / 0.36 / 0.49 | 9-11 | 1.0 / 1.6 / 2.0 |

The strict score still carries about 0.10 on the even communities, because
*S. argenteus* (realised 0.05) is pruned to 0 in every even sample. It cannot be a member
of any group, so *S. aureus*'s group of one has 0.10 and *S. argenteus* 0. It appears in
*S. aureus*'s `identical_16s_relatives`, which is what the lenient score credits.

**Which groups form** (τ = 0.25, `ambiguity_groups_formed.tsv`; 191 groups of two or
more rows over the 92 samples):

| Group | Samples | Sets | Group within 25% of its truth | Right members hidden per group |
|---|---|---|---|---|
| *E. coli* + *E. fergusonii* | 73 | every set | 73/73 | 0.42 |
| *K. pneumoniae* + *K. variicola* | 22 | log only | 22/22 | 1.0 |
| *B. subtilis* + *B. spizizenii* (± *B. stercoris*, other *Bacillus*; 9 member sets) | 45 | even | 45/45 | 0.98 |
| *S. pneumoniae* + *S. mitis* (± *S. pseudopneumoniae*, *S. chosunense*, …; 5 member sets) | 34 | even | 34/34 | 1.06 |
| *A. muciniphila* + *Lepagella* sp009775355 | 2 | even R9.4.1 1k | 2/2 | 1.0 |
| others (13 member sets) | 15 | 14 of them at 1k | 8/15 | 0.47 |

- **The pairs the species table gets wrong** (*E. coli* / *E. fergusonii* in every set,
  *Bacillus* and *Streptococcus* in even communities) become groups whose sum is within
  25% of the truth in every sample.
- **Separable pairs:** *K. pneumoniae* / *K. variicola* split correctly in even
  communities (0.485 / 0.499 in Phase 7a) and are never grouped there. In log
  communities (*K. variicola* 5-9 ranks below) they are grouped in 22 of 40 samples,
  each time hiding a member that was right.
- **Where the hidden species come from.** On even R10.4.1 the *Bacillus* and
  *Streptococcus* pairs are often within 25% for both members, yet still grouped: 3.3-4
  of the 16 right species are hidden. The cost falls where the species table was already
  close.
- **Groups across genera:** 0.03 per sample, 3 in all (*Akkermansia* / *Lepagella* twice, *Collinsella* / *Lactobacillus* once).
- **The other groups** (15) are mostly a community species with a false positive at 1k reads; 8 of them are within 25% of their truth, and the 4 made only of false positives have truth 0.

Groups also absorb false-positive species: a false positive that shares its reads with a
community species joins its group. Strict precision rises from 0.889 to 0.908 and false
positives fall from 2.66 to 1.90 per sample. Strict recall falls from 0.859 to 0.846 (84
of 92 samples lower). That fall is a counting effect: grouping a pair reduces the number
of present truth keys while *S. argenteus* stays missing.

**EM** (`ambiguity_em_vb.tsv`, 40 samples at 10k): the same groups and L1 within 0.0002
of VB in every set (strict even R10.4.1 0.1091 vs 0.1092).

**Intervals** (`ambiguity_coverage.tsv`, VB, 10k, τ = 0.25, truth the design abundance
renormalised over the reported species as in `intervals.py`):

| Level | Coverage | n |
|---|---|---|
| all reported species | 0.61 | 623 |
| members of groups of two or more rows (their species intervals) | 0.26 | 128 |
| those groups (summed-α intervals) | 0.69 | 64 |

The 0.61 for all reported species is not comparable with Phase 7a's 81% at 10k
(`methods_2026.md` §8): that figure left out the pair members (42%), and this one pools
both, on minimap2 2.31 alignments. The group intervals cover the group's truth far more often than the members' intervals
cover theirs, as the Phase 7a toy case predicted (the pair's sum is identifiable, its
split is not). They are still sampling intervals and still undercover.

## 4. Real mock communities

Means over the 13 real sup mocks (`ambiguity_real_summary.tsv`, `ambiguity_real_paired.tsv`;
L1 over the species with a known share):

| τ | species L1 | strict L1 | lenient L1 | strict precision | FP species (species → strict) | groups | reads in groups | right hidden (of 4.0) | groups across genera |
|---|---|---|---|---|---|---|---|---|---|
| 0.20 | 0.352 | 0.325 | 0.332 | 0.822 | 5.15 → 2.46 | 1.62 | 0.13 | 0.54 | 0.62 |
| **0.25** | 0.352 | **0.325** | 0.332 | 0.799 | 5.15 → 2.85 | 1.38 | 0.09 | 0.46 | 0.54 |
| 0.30 | 0.352 | 0.337 | 0.331 | 0.771 | 5.15 → 3.38 | 1.38 | 0.08 | 0.46 | 0.54 |
| 0.40 | 0.352 | 0.342 | 0.343 | 0.740 | 5.15 → 4.08 | 0.77 | 0.03 | 0.23 | 0.38 |

At τ = 0.25 the strict L1 is lower than the species L1 on 6 samples and higher on 3, all
of them by 0.001 (MSPlus reps 3-5). Precision rises on 12 of 13, recall is unchanged.
Per sample (`ambiguity_real_per_sample.tsv`):

| Sample | species L1 | strict L1 | groups | reads in groups |
|---|---|---|---|---|
| emu_gutmock_SRR14326315 | 0.278 | 0.164 | 4 | 0.36 |
| emu_zymo_SRR10391201 | 0.428 | 0.246 | 1 | 0.20 |
| preprint_zymo_SRR23636353 | 0.353 | 0.305 | 1 | 0.23 |
| zymo_GutMock_sup_rep1-4 | 0.592-0.623 | 0.587-0.619 | 0-1 | ≤ 0.01 |
| zymo_MSPlus_sup_rep1-6 | 0.168-0.193 | 0.168-0.194 | 1-2 | 0.01-0.11 |

The gains come from the three R9.4.1 mocks:
- **Zymo's *Bacillus*** (GTDB *B. spizizenii*) is grouped with *B. subtilis*,
  *B. stercoris*, *B. halotolerans* and *B. velezensis* in emu_zymo (group 0.201 against
  0.174) and in the preprint sample (0.234 against 0.174).
- **The Emu gut mock** gets *Dorea*, *Paraclostridium* and *Enterococcus* groups.

On the sup mocks groups are rare and small.

Groups on the real mocks (`ambiguity_groups_formed.tsv`, kind `real`) cross genera far
more often than in the simulations: 7 of the 18 groups do. Examples:
- *Enterococcus faecalis* + *Staphylococcus hominis* (class level, MSPlus reps 3-5,
  group 0.092 against 0.099);
- *Akkermansia muciniphila* + *Faecalibacterium taiwanense* (GutMock);
- *Romboutsia hominis* + *Parabacteroides merdae* + *Faecalicatena contorta* (Emu gut mock,
  three phyla, 0.135 against 0.122).

Every link inside these groups was checked (`ambiguity.py --cross-links`,
`ambiguity_real_cross_links.tsv`: 10 links at f ≥ 0.25). No link joins two mock species:
- 9 links join a mock species to a false positive with 10-23 expected reads (one more
  joins two false positives).
- Every species in them except *Faecalicatena contorta* carries a sequence that
  `cross_genus_identical.tsv` lists, the D23 sequences the `keep` database retains.
- In the Emu gut mock, *Romboutsia hominis* (a mock species) links to two false
  positives: *Parabacteroides merdae* (f = 0.27) and *F. contorta* (0.35).

The group's sum is close to the truth, but its label names a species that is not there,
at a rank as high as superkingdom. Two groups contain only false-positive
species: *B. subtilis* + *B. velezensis* in MSPlus, 0.014-0.021, truth 0.

The lenient score is not better than the strict one on the real mocks (0.332 against
0.325). On MSPlus it is worse than the species score (for example rep1 0.203 against
0.184), because the mocks' species are listed as relatives of one another's groups and
the joined groups are scored together.

## 5. e2e fixture

`nextflow run main.nf -profile test,docker --ambiguity_groups true`
(`phase14/runs/e2e_br_groups`, and `e2e_br_groups_em` with `--estimator em`):

- **The pair:** *E. coli* and *S. flexneri* (truth 0.25 and 0.05; 51 of their reads tie)
  form one group at 0.300 (VB 0.30006, EM 0.30000) with 180 reads.
- **Why they link:** f is 0.771 for *S. flexneri* (16.8 of its 21.8 expected reads) and
  0.351 for *E. coli*. The group's lowest common rank is the family Enterobacteriaceae and
  its `shared_read_fraction` is 0.
- **The other four species** are groups of one with the species table's values.
  *S. pneumoniae* / *S. mitis*, which the fixture splits correctly (0.150 / 0.050), share
  0.101 of *S. pneumoniae*'s reads and stay apart.

With the option off, the five tables and `software_versions.yml` are byte-identical to
dev @ 2f6c031: with the defaults, with `--estimator em`, and with
`--credible_intervals --binomial_names`. With it on, the species table and the
consolidated tables are byte-identical too.

## 6. Cost

RUN_VI on the Phase 6 synthetic set (200,000 scored reads over 2,000 species, 173
reported; VB, one core, two runs each): 3.86 / 3.76 s and 153 MB on dev, 3.91 / 3.73 s
and 155 MB on the branch without the option, 3.97 / 3.81 s and 162 MB with it. The JSON
is otherwise identical to dev's. No group formed there, since the synthetic reads' decoys
score lower.

## 7. Links limited to one family

A test of a change to D39, made after sections 3-6 and since adopted (D39 amended
2026-11-05): two rows link only when their species are in the same family of the
database's taxonomy. `ambiguity.py --family` runs it as one
more arm, `family_tau0.25` (VB, τ = 0.25, `--credible-intervals`), on all 105 samples:
105 runs, no failures. The run goes through `benchmarks/family_groups.py`, which runs
NanoVI's `estimate_abundance.py` at `9137b7b` and sets f to 0 between rows of different
families before the rows are linked. Families are transitive, so every group stays inside
one family. `analyse.py ambiguity_family` compares it with `vb_tau0.25`
(`ambiguity_family_summary.tsv`, `_paired.tsv`, `_changed.tsv`, `_coverage.tsv`).

**What changes.** The species tables are byte-identical to the default arm's in all 105
samples, and the group tables in 95. The 10 that change are the 10 with a group above the
family:
- 7 real mocks: the 7 cross-genus groups of section 4 (class or superkingdom level);
- 3 simulations (even R9.4.1 1k, seeds 2, 5, 9): *Akkermansia muciniphila* +
  *Lepagella* sp009775355 (twice, 0.043 against 0.044 and 0.055 against 0.057) and
  *Collinsella* sp938045065 + *Lactobacillus delbrueckii* (0.003, truth 0).

In this database every cross-genus group also crosses a family, so a genus rule would remove
the same groups. Other taxonomies differ: GTDB places *Shigella* in *Escherichia*, but the
e2e fixture's taxonomy keeps *Shigella* as its own genus. There, a genus rule would split
*E. coli* / *S. flexneri*; the family rule keeps them (both Enterobacteriaceae). The
fixture's RUN_VI JSON is byte-identical with and without the family rule.

| | default (D39) | family rule |
|---|---|---|
| Real mocks, strict L1 | 0.3251 | 0.3257 (4 lower, 3 higher) |
| Real mocks, lenient L1 | 0.3319 | 0.3334 |
| Real mocks, false positives per sample (strict) | 2.85 | 3.54 (species table 5.15) |
| Real mocks, groups above family | 7 | 0 |
| Simulations, strict L1 | 0.07516 | 0.07526 (2 higher) |
| Simulations, false positives per sample (strict) | 1.90 | 1.93 |
| Group interval coverage at 10k | 0.69 (64 groups) | 0.69 (the same groups) |

**Per sample** (strict L1):
- **Worse with the family rule:** *Akkermansia*'s groups. GutMock reps 1, 3 and 4 rise by
  0.0040, 0.0040 and 0.0053. There the false positives' reads belong to *Akkermansia* (the
  group's 0.0075-0.0093 against its truth 0.0097, the species row alone 0.005-0.007), so
  the group was nearer the truth. The two simulated *Akkermansia* + *Lepagella* groups
  rise by 0.0051 and 0.0039 for the same reason.
- **Better with the family rule:** MSPlus reps 3-5 fall by 0.0010-0.0013, and the Emu gut
  mock falls by 0.0024 (its *Romboutsia* group was 0.135 against 0.122).
- **False positives:** 9 come back as their own rows over the 7 real mocks. They keep
  their high `shared_read_fraction` in both tables.

**The shipped rule.** NanoVI commit `3b24374` implements the rule in `bin/ambiguity_groups.py`.
Its families come from `taxonomy.tsv` with the ranks above them, so a family name used in two
orders is two families. A blank or missing family links to nothing. `ambiguity.py` (all 460
runs, `NANOVI_PHASE14=phase14/confirm_3b24374`) reproduces `family_tau0.25` at
`vb_tau0.25`: in all 105 samples the species tables are byte-identical, and so are the group
tables' groups, members, abundances, counts, ranks and `max_shared`. `shared_read_fraction`
and the intervals differ by at most 2.2e-16, as they do in samples the family rule does not
change (the review fix of section 2). The e2e fixture's tables with the option off are
byte-identical to dev @ 2f6c031.

**Reading.** On these data the family rule is neutral on accuracy: means within 0.0006
(real) and 0.0001 (simulations), in both directions per sample. It costs 0.7
false-positive rows per real mock. What it removes are the groups whose label is wrong in
kind: in each of them a database sequence (D23) makes a species from another phylum or
class look like a 16S look-alike.

## 8. Summary

- **What the group table does.** It reports what the reads support where the species
  table splits by abundance. On the simulations at τ = 0.25, L1 falls from 0.132 to 0.075
  (strict, 83/92 samples lower, none higher) and to 0.0096 when the database relatives
  count (92/92). The *E. coli* / *E. fergusonii*, *Bacillus* and *Streptococcus* groups
  are within 25% of their truth in all 152 samples where they form. On the real mocks the strict L1
  falls from 0.352 to 0.325, almost all of it on the three R9.4.1 mocks.
- **The costs:**
  - 1.5 of the 12.1 right species per simulated sample are merged into a group, more
    where both members were already close (even R10.4.1: 3.3-4 of 16).
  - Separable pairs are grouped when one member is much rarer (*Klebsiella* in 22/40
    log samples).
  - On the real mocks, groups cross genera through false-positive species in 7 of 18
    cases.
- **The threshold.** τ = 0.2 scores slightly lower (simulations 0.068, real 0.325) but
  hides more right species (1.8) and forms more groups across genera. τ = 0.3 gives up
  most of the real-mock gain (0.337). 0.25 (D39) is supported.
- **Intervals:** the group intervals cover 0.69 of group truths at 10k, against 0.26 for
  their members' intervals.
- **The cross-genus groups** all come from a false positive linked to a mock species,
  almost always through the cross-genus identical sequences the `keep` database retains
  (D38).
- **Links limited to one family** (section 7; D39 amended to it) remove all 10 groups
  above the family. Mean L1 changes by at most 0.0006, and the real mocks gain 0.7
  false-positive rows per sample. *E. coli* / *S. flexneri* stays grouped on the e2e
  fixture. With the rule as shipped, the strict L1 is 0.075 on the simulations and 0.326
  on the real mocks.
- **Not addressed here:** the species table still reports the split, and
  `shared_read_fraction` still flags it.
