# Ambiguity groups (design note, F57)

Status: proposal, waiting for decision D39 (Phase 14). Nothing here is implemented yet. The numbers below are a preliminary look at the Phase 13 runs, made to choose the options; the Phase 14 benchmark (`docs/benchmarks/ambiguity_2026.md`) redoes them with the shipped code. Preliminary scripts and tables: `~/pipelines/NanoVI_benchmarks/phase14/prelim/` (`pair_sharing.py`, `group_rules.py`, `group_rules2.py`, `group_rules3.py`, `ties.py`), run on NanoVI dev @ 2f6c031 `bin/` against the Phase 13 score files (`phase13/runs/<sample>/lik/exact.npz`: the Phase 7a exact database, merged GTDB r226 build, `--cross-genus-rule keep` (D38), minimap2 2.31, `-N 3`), refitted with the default VB.

## Problem

When species' 16S genes are identical or nearly so, their reads carry little or no evidence for one species over the other, and the estimator splits such reads by abundance: the more abundant species takes them. The species table then reports each species alone, with an abundance that the data do not support. Phase 7a's examples (D24, `methods_2026.md` §5): *S. argenteus* is estimated at 0 beside *S. aureus* in every even simulation (realised 0.05 each); *E. coli* takes 0.77 of the *E. coli* / *E. fergusonii* pair at `-N 3`; on the e2e fixture 51 *E. coli* / *S. flexneri* reads tie, and *S. flexneri* is estimated at 0.037 against 0.05. Two columns already flag this: `identical_16s_species` (the database's view, D24) and `shared_read_fraction` (the sample's view, D27). Neither reports what the sample does support, the pair's sum.

F57 asks for a per-sample table in which species that the sample cannot separate are reported together. The species table stays as it is (the prompt's step 2).

## What a group can be built from

### (a) The database's identical-16S table (`identical_16s.tsv`, D24)

The file has two notions, and they behave very differently:

- **Components** (`group`, `group_size`): species of a genus linked through shared sequences, directly or through a chain. They are large: in the merged r226 build *Pseudomonas* G2 has 137 species, *Streptococcus* G5 60, *Bacillus* G6 42; 18 of the 20 simulated species are in one. D24's "group L1 0.010 against species L1 0.131" was scored over these components (`evaluate.py`, from the Phase 7a audit's groups). Scored that way here (truth and estimate summed over components, unreported members included): L1 0.0098, close to D24's figure, but 95% of the abundance goes into groups and 10.9 of the 13.8 species per sample that the species table had right (within 25% of the truth) disappear into a group. A component is close to a genus (Phase 7a genus L1 0.008) and is not a usable row.
- **Direct links** (`identical_16s_species`): the species that share at least one identical sequence with a given species. Smaller, but 18 of the 20 simulated species still have one (*P. aeruginosa* 4, *E. faecalis* 6, *S. mitis* 20). Grouping reported species that are directly linked: L1 0.061, 45% of the abundance in groups, 3.4 correct species per sample hidden. It groups *K. pneumoniae* with *K. variicola*, which the sample splits correctly (0.485 vs 0.499 in Phase 7a): one shared copy among several does not make two species inseparable.

The database's view is the same in every sample, so groups built from it line up across samples. Its weakness is that sharing one copy is not the same as being inseparable: it cannot tell the *Klebsiella* pair (separable) from the *E. coli* pair (not).

### (b) Reads the sample's fit shares between two rows

RUN_VI has the final responsibilities φ (VB; `em_estimator.responsibilities` for EM, D2). For two reported table rows a and b (rows as the species table forms them, `species_row_keys`, so `--binomial_names` merges first, F46):

  f_a(b) = Σ_r φ_ra · [φ_rb ≥ 0.01] / E_a,   E_a = Σ_r φ_ra

the share of a's expected reads that sit on reads b takes at least 1% of. Rows a and b are linked when **max(f_a(b), f_b(a)) ≥ τ**, and groups are the connected components of the links. The maximum, not the minimum, because the case to catch is a small row whose estimate rests on reads it shares with a large one (*S. flexneri*: f = 0.77; *E. coli* beside it: 0.35).

This is not `shared_read_fraction`. That column counts a read as shared when all other rows together take 0.01 of it; f is per partner. For a row with one partner they nearly agree (e2e *S. flexneri*: 0.771 both), with several partners f is smaller (the D37 correction about "the same set of reads" applies here too).

Preliminary distribution over the 92 simulations (`pair_sharing.tsv`; the larger f of each reported pair):

| Pair (both reported) | even R10.4.1 | even R9.4.1 | log R10.4.1 | log R9.4.1 |
|---|---|---|---|---|
| *E. coli* / *E. fergusonii* | 0.46 (0.28-0.55) | 0.33 (0.10-0.50) | 0.30 | 0.31 |
| *B. subtilis* / *B. spizizenii* | 0.30 (0.25-0.56) | 0.33 (0.27-0.46) | (1 sample) | (1 sample) |
| *S. pneumoniae* / *S. mitis* | 0.33 (0.18-0.38) | 0.25 (0.18-0.42) | — | — |
| *K. pneumoniae* / *K. variicola* | 0.06 (0.02-0.08) | 0.09 (0.04-0.16) | 0.22 | 0.29 |
| *S. aureus* / *S. argenteus* | never both reported | | | |

Median (range); "—": the second member is not reported (log communities place it 5-9 ranks lower). Every other reported pair (1,505 pairs): median 0.002, 95th percentile 0.29, 99th 0.50; most of the high values are rows of 1-3 expected reads at 1k reads beside a community species (37 of the 93 above 0.25 have a member with fewer than 5 expected reads), and 3 cross a genus (*Akkermansia muciniphila* / *Lepagella* sp009775355 0.48). At τ = 0.25, 171 of 275 design-pair observations link and 93 other pairs; at τ = 0.3, 118 and 62. The pairs whose split is wrong in the species table (*E. coli*, *Bacillus*, *Streptococcus* in even communities; mean split error 0.03-0.06) sit at 0.25-0.46, and *Klebsiella* in even communities (split error 0.007) at 0.06-0.08.

Grouping by (b) alone (truth summed over the reported members only):

| τ | species L1 | group L1 | abundance in groups | correct species hidden / sample |
|---|---|---|---|---|
| none | 0.132 | — | 0 | 0 |
| 0.20 | | 0.069 | 0.34 | 1.8 |
| 0.25 | | 0.075 | 0.29 | 1.5 |
| 0.30 | | 0.091 | 0.21 | 1.1 |
| 0.40 | | 0.116 | 0.06 | 0.25 |

Its structural limit: **(b) can only link rows that are reported.** *S. argenteus* is pruned to 0 in every even simulation, so no in-sample rule built on φ can form an *S. aureus* / *S. argenteus* row, and that pair alone leaves about 0.10 of L1 on the even communities. Ties in the likelihood do not rescue it: the share of *S. aureus*'s reads whose best score *S. argenteus* equals is 0.24 at R10.4.1 and 0.05 at R9.4.1 (`ties.tsv`), because at `-N 3` (D21) the identical copy's alignment is often not among the alignments minimap2 keeps. Only the database knows that *S. argenteus* could be there.

The other cost of (b): groups depend on the sample, so two samples can group the same species differently, and the per-sample group tables cannot be lined up row by row.

On the e2e fixture (no identical sequences in its 11-sequence database, so this is the pure (b) case): *S. flexneri* f = 0.771 (16.8 of its 21.8 expected reads), *E. coli* 0.351; *S. pneumoniae* / *S. mitis*, which the fixture splits correctly (0.150 / 0.050 against 0.15 / 0.05), 0.101. At τ = 0.25 or 0.3 the table would have one group, *E. coli* + *S. flexneri*, at 0.300 against a truth of 0.30.

### (c) Both

Rows from (b); each row annotated with the database's direct relatives of its members that are not reported in the sample (from `identical_16s_species`). The annotation says "these species could be in this row too"; it does not merge reported rows, so it cannot group the *Klebsiella* pair.

Two ways to score it against the truth: **strict** (the truth summed over the reported members only, as for (b)) and **lenient** (the annotated relatives' truth counts towards the row; a relative named by two rows joins them for scoring). Preliminary (`group_rules3.tsv`):

| τ | strict L1 | lenient L1 | lenient, even | lenient, log |
|---|---|---|---|---|
| no rows from (b), annotation only | 0.132 | 0.021 | 0.031 | 0.010 |
| 0.25 | 0.075 | 0.0096 | 0.009 | 0.010 |
| 0.30 | 0.091 | 0.011 | 0.012 | 0.010 |

The lenient score is generous by construction (a row of *S. aureus* with *S. argenteus* in its annotation gets *S. argenteus*'s truth), but it is the honest summary of what the table would tell a reader: "*S. aureus*, 0.10, which may include *S. argenteus*, *S. schweitzeri* or *S. singaporensis*".

## Recommendation for D39

(c): groups from in-sample sharing (b) with **τ = 0.25** (0.3 as the alternative; the benchmark decides), annotated with the unreported database relatives (a, direct links only, never components), written to a separate per-sample table, **off by default** behind one parameter. Reasons: (b) is the only signal that separates the *Klebsiella* case from the *E. coli* case; (a) is the only one that can name *S. argenteus*; components are too coarse to be rows.

## The group table: `sample_outputs/<sample>_groups.tsv`

A partition of the species table: every species row is in exactly one group (a row that links to nothing is a group of one), then the `Below detection` and `Unassigned` rows as in the species table, so the group table sums to 1 (abundance) and to the filtered reads (counts) as the species table does.

| Column | Meaning |
|---|---|
| `group` | `G1`, `G2`, … by decreasing abundance; the placeholder rows keep their labels |
| `members` | The species rows in the group, named as in the species table, `; `-separated, by decreasing abundance. This is the label: a list, because a lowest common rank throws away what the reader most needs (which species), and (b) can cross genera |
| `n_members` | Number of species rows |
| `abundance` | Sum of the members' abundances (exactly the species table's values summed) |
| `estimated_counts` | With `--keep_counts`: sum of the members' `estimated_counts` (so the group table's counts add up as the species table's do) |
| `lowest_common_rank`, `lowest_common_taxon` | The lowest rank at which every member's lineage agrees (`species` for a group of one), and the taxon there |
| `max_shared` | The largest f between two members (empty for a group of one): how strongly the group is held together |
| `shared_read_fraction` | The group's share of expected reads split with rows **outside** the group (the D27 rule, judged between groups). It falls when the group absorbs the partner it shared with: for e2e *E. coli* + *S. flexneri* it would be about 0 |
| `identical_16s_relatives` | Species of the members' genera that share an identical 16S with a member in the database and are not reported in the sample, `; `-separated (empty without `identical_16s.tsv`, as the species column is under `--db_index` without a FASTA) |
| `ci_lower`, `ci_upper` | With `--credible_intervals` (VB): from the summed α, Beta(Σ_group α, A − Σ_group α), as summed rows already get. The Phase 7a toy case (`intervals_toy.tsv`) showed the pair's sum is the identifiable quantity, so these should cover better than the members' intervals; the benchmark measures it |

The members' own rows stay in the species table unchanged.

### Where it is computed

The group-level `shared_read_fraction` needs each read's φ split between groups, which pairwise sums cannot give without double counting. So RUN_VI forms the groups: `estimate_abundance.py` (given τ, the taxonomy and `--binomial-names`, as it already is for F46) computes f between rows from φ, forms the components, and writes them into the abundance JSON as a new key, only when the option is on (`"groups": [{"rows": [tax IDs], "expected_reads", "shared_reads", "max_shared"}]`; the JSON stays byte-identical to dev's without it). WRITE_OUTPUT (`write_abundance_table.py --groups`) turns it into `<sample>_groups.tsv` with the lineages, counts, intervals and the relatives from `identical_16s.tsv`. No per-read data reaches WRITE_OUTPUT. The cost is one pass over φ (R × K entries) and an S × S matrix over reported rows (S is tens), negligible beside the fit.

### Other outputs

- **Species table, `collapse-taxonomy`, `combine-outputs`, consolidated tables**: unchanged. They work on the species table. `collapse-taxonomy` to genus already reports what a genus-level reader needs.
- **No consolidated group table.** Groups from (b) differ between samples, so their rows cannot be lined up; a cross-sample union (components over all samples' links) is possible later, but it would move groups away from what each sample shows. The group table is per sample only, and the consolidated tables keep species rows.
- **Per-read file (Phase 12, `--read_assignments`)**: unchanged. A group's reads are the union of its members' lines; their probabilities summed over the members give the group's expected reads (less what falls below 0.01). Lines stay per species row, so the file does not depend on τ.
- **EM**: the same computation from `em_estimator.responsibilities` (D2); no intervals (as now).
- **Phase 15 (copy numbers)**: a group's copy number is the abundance-weighted mean of its members' (the Phase 15 prompt says so); `members` names the rows it needs.

### Parameters

`--ambiguity_groups` (boolean, default false) writes the table. τ: either fixed at 0.25 with a script-only option for tests (as D37 fixed 0.01), or a second pipeline parameter `--ambiguity_min_shared` (default 0.25). I'd fix it: one parameter, and a value the benchmark sets.

## Questions for D39

1. **Definition**: (a) direct links, (b) in-sample sharing, or (c) (b) with (a) as annotation (recommended). Components are not offered as rows.
2. **Threshold** τ on max(f_a(b), f_b(a)): 0.25 (recommended) or 0.3, final value after the benchmark; fixed, or a parameter.
3. **Output**: the separate per-sample table above (recommended; the species table never changes), or rows in the species table.
4. **Default**: off (recommended until the benchmark), or on.
5. **Benchmark database**: the Phase 7a merged exact build, as Phase 13 ran (recommended: same runs, comparable to D24 and D38), plus the e2e fixture. The shipped default build keeps GTDB suffixes (`identical_16s.tsv` 5,215 species, against 4,335 merged), where suffixed siblings such as *L. monocytogenes* / *L. monocytogenes_B* would also group; checking that needs the Phase 7b `suffixes2` runs.

## Benchmark plan (step 3)

Through `estimate_abundance.py` and `write_abundance_table.py` of the NanoVI branch (NANOVI_REPO), on the stored Phase 13 score files: the 92 simulations (even and log, R9.4.1 and R10.4.1, 1k/10k at seeds 1-10 and 50k at seeds 1-3) and the 13 real sup mocks (`phase13/real/runs/*/lik/exact.npz`), τ ∈ {0.2, 0.25, 0.3, 0.4}, VB (EM on the 10k samples as a check). Scores: L1, precision and recall of the species table and of the group table, against the truth summed over the same groups, strict and lenient; reads moved into multi-member groups (`estimated_counts`); species that the species table had right (within 25% of the truth, a definition to keep fixed) that end up in a multi-member group; groups that cross a genus; coverage of the group intervals at 10k. On the e2e fixture, the *E. coli* / *S. flexneri* group (51 tied reads). Byte identity of the pipeline outputs with the option off.
