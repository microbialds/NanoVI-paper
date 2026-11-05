# 16S copy-number correction (Phase 15, F58, November 2026)

NanoVI's abundances are shares of 16S reads, which follow shares of 16S genes: a species
with ten copies of the gene per genome gives about ten times the reads per cell of a
species with one. This report compares sources of per-species copy numbers for the GTDB
r226 database, then scores NanoVI's optional cell-level table (`--cell_abundance`:
read-level abundance / copy number, renormalised) against cell-level truth on real mock
communities and on simulations whose reads are drawn in proportion to cells × copies. The
choices it informed are decision D40; the read-level table does not change and the
option is off by default.

Every table cited here is in `docs/benchmarks/data/` and the figure in
`docs/benchmarks/figures/`, with its source data under the same name in `data/`. Data,
databases, alignments and per-sample outputs are under
`~/pipelines/NanoVI_benchmarks/phase15` (D11) and are not in the repository.

## 1. Summary

- **Source.** GTDB r226's `ssu_count` (16S genes HMMER finds with Barrnap's models; the
  genome's records in `ssu_all`) is usable only for complete genomes: draft assemblies
  collapse their rRNA operons (72% count 0 or 1). The median over a species' `Complete
  Genome` assemblies covers 14.2% of the default build's 84,596 species, against 11.5%
  for rrnDB 5.10 matched by assembly accession, and the two agree on 99.1% of the 43,167
  genomes they share. It gives Zymo's own copy numbers for 19 of the 20 rows of the two mock tables.
- **Filling.** A species without a value takes the median of its genus, else of its
  family (D40): leave-one-out error 0.45 copies (genus, 91% within one copy) and 0.81
  (family, 82%). 62% of species are filled this way and 24% have no value; weighted by
  abundance, 59% of a stool cohort's assigned abundance has a species-level value, 98.6%
  of the real mocks'.
- **Real mocks** (12 with a cell-level truth): the cell table's species L1 against the
  cell truth is 0.367, against 0.438 for the read-level table (lower on 11 of 12). This
  only checks that the database's copy numbers agree with Zymo's: Zymo derives its cell
  shares from the same 16S shares and copy numbers, and the cell table's error (0.367)
  is the read-level table's error against the 16S truth (0.359) carried over.
- **Simulations** (40, reads ∝ cells × copies of one complete genome per species): cell
  table 0.130 against 0.343 for the read-level table (lower on 40 of 40); dividing by the
  simulated genome's own copy number gives 0.116, so strain-to-species differences in
  copy number add 0.014, and the read-level table's own error against its read truth is
  0.121. Genus L1: 0.027 (cell table), 0.265 (read table).
- **Filled values.** Masking the community species' own values and filling them from
  their genus raises L1 by 0.106 on the simulations and 0.033 on the real mocks: a genus
  median can be far from the species (*Pseudomonas* 6 for *P. aeruginosa*'s 4,
  *Enterococcus* 6 for *E. faecalis*'s 4), and the masked table is still better than the
  read-level table on every sample but one.

## 2. Sources (step 1)

`benchmarks/copy_number_sources.py` (`compare`, then `weighted`) reads the GTDB r226
`bac120_metadata_r226.tsv.gz` and `ar53_metadata_r226.tsv.gz` (732,475 genomes; already
under `phase7a/downloads/gtdb_r226`) and rrnDB 5.10 (release of 2025-04-07, NCBI taxonomy
of 2025-03-24: `rrnDB-5.10.tsv.zip`, 1.26 MB, 48,700 records; and
`rrnDB-5.10_pantaxa_stats_NCBI.tsv.zip`, 316 KB). Both builds of the Phase 13 GTDB r226
database are covered: the default (GTDB's suffixes kept, 84,596 species) and the merged
build (`--merge-gtdb-suffixes`, 81,173 species) that the Phase 7a-14 benchmarks run on.

**What `ssu_count` counts** (`copy_number_sources_assembly_level.tsv`). The metadata's
field description gives "Number of 16S rRNA genes, HMMER w/ Barrnap HMMs"; it equals the
genome's records in `ssu_all_r226.fna.gz` for every one of the 732,475 genomes. Only
complete genomes give a usable count:

| NCBI assembly level | Genomes | Mean `ssu_count` | `ssu_count` ≥ 1 |
|---|---|---|---|
| Complete Genome | 49,324 | 5.34 | 100% |
| Chromosome | 6,750 | 3.80 | 94% |
| Contig | 428,853 | 1.18 | 63% |
| Scaffold | 247,548 | 1.04 | 62% |

In the species that have complete genomes, 87% of their draft genomes count fewer copies
than the species' complete-genome median and 72% count 0 or 1. `Chromosome`-level
assemblies are left out too.

**Sources compared, per species** (`copy_number_sources_coverage.tsv`):

| Source | Default build | Merged build |
|---|---|---|
| GTDB complete genomes: median `ssu_count` | 11,971 (14.2%) | 11,124 (13.7%) |
| rrnDB 5.10 by assembly accession (median over the species' rrnDB genomes) | 9,748 (11.5%) | 9,040 (11.1%) |
| rrnDB 5.10 through NCBI taxonomy (pan-taxa median of the NCBI species most of the GTDB species' genomes carry) | 10,483 (12.4%) | 7,912 (9.7%) |
| NCBI's own annotation (`ncbi_ssu_count`), complete genomes | 10,783 (12.7%) | 9,997 (12.3%) |
| GTDB, every genome (drafts included) | 84,590 (99.99%) | 81,167 (99.99%) |

rrnDB's records are NCBI RefSeq assemblies: 48,485 of 48,700 are assembly accessions, and
43,344 of them are GTDB r226 genomes (43,167 complete), so rrnDB is a subset of the
complete genomes GTDB already has.

**Agreement** (`copy_number_sources_agreement.tsv`):

| Pair | n | Equal | Within 1 | Pearson |
|---|---|---|---|---|
| Per genome: GTDB `ssu_count` vs rrnDB (complete genomes) | 43,167 | 99.13% | 99.83% | 0.998 |
| Per genome: NCBI `ncbi_ssu_count` vs rrnDB | 43,136 | 99.74% | 99.98% | 1.000 |
| Per species: GTDB complete vs rrnDB by accession | 9,731 | 98.95% | 99.79% | 0.998 |
| Per species: GTDB complete vs rrnDB through NCBI taxonomy | 7,187 | 95.14% | 98.72% | 0.987 |
| Per species: GTDB every genome vs GTDB complete | 11,971 | 65.6% | 72.8% | 0.678 |

rrnDB through NCBI taxonomy copies one NCBI species' value to every GTDB species it is
split into (*Faecalibacterium prausnitzii* over 51 GTDB species), which is where most of
its disagreement comes from. Drafts would cover every species but under-count by 1.17
copies on average.

**Against Zymo.** The D6300 and D6331 manuals' Table 1 gives each strain's "16S Only" and
"Genome Copy" shares; their ratio is the strain's copy number (scaled to *E. coli* = 7:
4.06, 7.00, 7.04, 5.02, 3.99, 6.01, 5.98, 9.95 for D6300). The GTDB complete-genome
species value equals it, rounded, for 19 of the 20 rows of the two tables with a 16S
share above 1e-4 (`copy_number_values.tsv`); *Bifidobacterium adolescentis* is 4 in
every source and 5 in Zymo's table.

**Filling species without a value** (`copy_number_loo.tsv`). D40's rule: the median of
the species values of the species' genus, else of its family; ranks are compared with the
ranks above them, so a genus name used in two families is two genera; a species with
neither has none. Leave-one-out over the species with their own value (each filled from
the others):

| Build | Filled from | n | Mean abs. error (copies) | Median rel. error | Within 1 copy | Mean signed |
|---|---|---|---|---|---|---|
| default | genus | 9,958 | 0.452 | 0 | 90.9% | −0.03 |
| default | family | 1,537 | 0.814 | 0.10 | 81.5% | +0.09 |
| merged | genus | 9,269 | 0.507 | 0 | 89.5% | −0.04 |
| merged | family | 1,376 | 0.884 | 0.125 | 80.3% | +0.23 |

Order medians would be off by 1.00 copies (76% within one) and were left out. These
errors are an optimistic bound: the species with a value sit in well-sequenced genera,
and the species that need filling are in sparser ones. In the default build the rule
gives 11,971 species their own value, 31,049 a genus median, 21,629 a family median and
leaves 19,947 (23.6%) without one.

**Weight in real samples** (`copy_number_sources_weighted.tsv`; default build, GTDB
complete genomes, mean share of each sample's assigned abundance):

| Samples | Species value | Genus | Family | None (D40) |
|---|---|---|---|---|
| 13 real mocks (Phase 7b, default build) | 98.6% (min 81.6%) | 0.4% | 1.0% | 0% |
| 117 CRC and healthy stool samples (Phase 9) | 59.4% (min 16.3%) | 25.3% | 14.5% | 0.8% |
| 12 vaginal samples, BV study (Phase 9) | 84.7% | 7.8% | 7.6% | 0% |
| 20 vaginal samples, preprint cohort (Phase 9) | 85.6% | 11.4% | 3.0% | 0% |

D40 (asked with these numbers, 2026-11-05): GTDB complete genomes; genus, then family
median; a separate per-sample table. `main.py build-database` (Emu-style input) takes a
user table or rrnDB's pan-taxa file looked up by NCBI taxid.

## 3. How the benchmark ran

- **Code.** NanoVI commit `92e3fbd` (branch `roadmap/p15-copy-number`) for every step
  (`NANOVI_REPO` pointed at `phase15/br_src`, a `git archive` of that commit):
  `build_gtdb_db.py --gtdb-metadata`, `estimate_abundance.py` and
  `write_abundance_table.py --copy-numbers --cell-output`. NanoVI-paper scripts at
  `317b571`: `benchmarks/copy_number_sources.py` and `benchmarks/copy_number.py`.
- **Database.** The merged GTDB r226 build rebuilt by that commit with
  `--merge-gtdb-suffixes --gtdb-metadata` (`phase15/db/merged_cn`): `species_taxid.fasta`,
  `taxonomy.tsv`, `identical_16s.tsv`, `cross_genus_identical.tsv` and
  `merged_species.tsv` are byte-identical to the Phase 13 `merged_keep` build (and the
  taxonomy to Phase 7a's exact database; `copy_number.py check-db`), plus
  `copy_numbers.tsv` (11,124 species values, 29,664 genus, 20,469 family, 19,916 none).
- **Tools.** fastplong 0.7.0, minimap2 2.31 + samtools 1.24 at `-N 3` against the exact
  index (the pipeline since Phase 7b; `harness.cross_run`), VB at the defaults. Badread
  0.4.2. Python 3.11.17 (numpy 2.4.6, pandas 3.0.6, scipy 1.17.1).
- **Real mocks.** The 12 with a cell-level truth: ONT Zymo 2025.09 MSPlus sup (6
  replicates) and the Emu and preprint Zymo runs SRR10391201 and SRR23636353 (D6300 /
  D6305 composition), and ONT Zymo 2025.09 GutMock sup (4 replicates, D6331). The score
  files Phase 13 stored for the exact database (`phase13/real/runs/<sample>/lik/exact.npz`)
  were fitted again with VB. The Emu gut mock has no cell-level truth and is left out.
- **Truth.** `benchmarks/truth/zymo_d6305_16s.tsv` and `zymo_d6331_16s.tsv` gain
  `genome_copy` and `cell_abundance`: Table 1's "Genome Copy" and "Cell Number" columns
  (D6300 manual Ver. 1.3.0, revised 2023-09-21; D6331 manual Ver. 1.2.0, revised
  2023-01-25; files.zymoresearch.com), renormalised over the species with a 16S gene (the
  yeasts left out, D6331's five *E. coli* strains summed). Zymo computes them from the
  genomic DNA composition (cells = DNA / genome size / ploidy; 16S = genome copies × 16S
  copies per genome), so for these bacteria the cell truth is the 16S truth divided by
  Zymo's copy numbers. The MSPlus spike-ins are present only (no share) and do not enter
  the L1.
- **Simulations** (`copy_number.py simulate`; `phase15/sims`). Phase 7a's 20-species
  community (five near-neighbour pairs and ten others), each species from one RefSeq
  `Complete Genome` assembly (seed 15) whose 16S copies in `ssu_all` are all at least
  1,400 bp, so the copies Badread draws from are the genome's `ssu_count`; genomes with a
  copy also found in another genus are skipped, except for *S. argenteus*, all of whose
  complete genomes (16 RefSeq, 1 GenBank) have one (`sims/references/genomes.tsv`). Cells even (5% each) or
  log-staggered (Phase 7a's log weights), 1,000 and 10,000 reads, R9.4.1 and R10.4.1
  Badread models, seeds 1-5: 40 samples. Each species' reads are a multinomial draw with
  p ∝ cells × copies; the truth is the realised cell shares, reads_i / copies_i
  renormalised. Four simulated genomes differ from their species' database value
  (*B. spizizenii* 8 vs 10, *C. difficile* 14 vs 12, *L. monocytogenes* 5 vs 6,
  *S. argenteus* 5 vs 6); the other 16 equal it.
- **Arms.** `read`: the species table (read-level); `db`: the cell table with the
  database's `copy_numbers.tsv`; `oracle`: the cell table with the truth species' own copy
  numbers (the simulated genome's count; for the mocks the Zymo-implied value, rounded,
  and the database's for the three D6331 species below 1e-4, whose rounded shares give no
  ratio); `masked`: the cell table with the truth species' own values removed (all of a
  community's at once) and filled again by D40's rule. Other species keep the database's
  values in every arm.
- **Metrics.** Species and genus L1 and Bray-Curtis (L1 / 2): `evaluate.compare` on the
  simulations, `evaluate.compare_real` on the mocks (the estimate renormalised over the
  species with a known share, as in every earlier phase).
- **Run.** 2026-11-05: genomes and references chosen by 16:11, the 40 samples' reads
  written by 16:24 (16 Badread processes); then 52 samples, 4 workers × 4 minimap2
  threads (16 cores, D20), 17:13-18:41. No failures.

## 4. Real mock communities

Species L1 against the cell truth (`copy_number_summary.tsv`, `copy_number_paired.tsv`;
means over samples):

| Samples | Read table vs 16S truth | Read table | Cell table (`db`) | `oracle` | `masked` |
|---|---|---|---|---|---|
| D6300 MSPlus (6) | 0.180 | 0.299 | 0.179 | 0.179 | 0.232 |
| D6300 SRR10391201 | 0.428 | 0.406 | 0.388 | 0.388 | 0.400 |
| D6300 SRR23636353 | 0.353 | 0.378 | 0.386 | 0.386 | 0.384 |
| D6331 GutMock (4) | 0.610 | 0.670 | 0.640 | 0.650 | 0.656 |
| All 12 | 0.359 | 0.438 | 0.367 | 0.371 | 0.400 |

The cell table is lower than the read table on 11 of 12 samples (mean −0.071); the
exception is SRR23636353 (+0.008), where the read-level estimate's own errors happen to
offset the copy-number difference. On D6300 the database's copy numbers equal Zymo's, so
`db` and `oracle` are identical; on D6331 the database's 4 for *B. adolescentis* (Zymo:
5) gives the lower L1 on all four replicates (−0.010). Genus L1: 0.439 (read table) and
0.350 (cell table).

The cell table's error is the read-level table's error carried over: 0.367 against the
read table's 0.359 against the 16S truth, on average (per sample the two differ by
−0.040 to +0.033). The GutMock's large L1 is the read-level estimate's (Phase 13/14):
*Fusobacterium nucleatum* and *Methanocatella smithii* are not detected in any arm, and
*Veillonella rogosae* is overestimated (0.27 of the known-share abundance in the read
table against a cell truth of 0.20, 0.34 in the cell table: dividing by its 4 copies
raises it further; `copy_number_species.tsv`). This panel cannot show more than that the copy numbers are right, because the
truth was made with them.

## 5. Simulations

Species L1 against the realised cell shares (means over 5 seeds per row):

| Cells, chemistry, reads | Read table vs read truth | Read table | Cell table (`db`) | `oracle` | `masked` |
|---|---|---|---|---|---|
| even, R10.4.1, 10k | 0.152 | 0.406 | 0.182 | 0.153 | 0.256 |
| even, R10.4.1, 1k | 0.156 | 0.412 | 0.189 | 0.158 | 0.259 |
| even, R9.4.1, 10k | 0.276 | 0.496 | 0.284 | 0.259 | 0.353 |
| even, R9.4.1, 1k | 0.297 | 0.513 | 0.300 | 0.279 | 0.365 |
| log, R10.4.1, 10k | 0.009 | 0.228 | 0.011 | 0.009 | 0.161 |
| log, R10.4.1, 1k | 0.011 | 0.226 | 0.013 | 0.011 | 0.158 |
| log, R9.4.1, 10k | 0.029 | 0.231 | 0.027 | 0.025 | 0.165 |
| log, R9.4.1, 1k | 0.037 | 0.236 | 0.032 | 0.031 | 0.167 |
| All 40 | 0.121 | 0.343 | 0.130 | 0.116 | 0.236 |

The cell table is lower than the read table on all 40 samples (mean −0.214; genus L1
0.265 → 0.027). With the genome's own copy numbers (`oracle`) it is 0.116, the
read-level table's own error (0.121) carried over; the database's species values add
0.014 (higher than `oracle` on 40 of 40), from the four species whose simulated genome
differs from the species median (`copy_number_species.tsv`: *C. difficile* 14 vs 12,
mean |error| 0.0049 against 0.0004 with its own value; *L. monocytogenes* 5 vs 6, 0.0045
against 0.0005; *B. spizizenii* 8 vs 10, 0.0108 against 0.0072). The even communities'
L1 stays high in every arm because of the near-neighbour pairs the read-level estimate
splits wrongly (*S. argenteus* is pruned beside *S. aureus* in every sample, as in Phase
14).

## 6. Filled values

`masked` against `db` (the same samples with every community species' own value replaced
by its genus median; no family fill was needed):

- Simulations: +0.106 species L1 (higher on 40 of 40); the largest per-species errors are
  *Pseudomonas aeruginosa* (genus median 6 for its 4: mean |error| 0.039 against 0.001),
  *Escherichia coli* (genus median 7, the same value, but other species' shifts move its
  share: 0.030 against 0.009), *Enterococcus faecalis* (6 for 4: 0.015 against 0.001) and
  *Bifidobacterium longum* (3 for 4: 0.009 against 0.0004).
- Real mocks: +0.033 (higher on 11 of 12): *P. aeruginosa*, *E. faecalis* and *Roseburia
  hominis* take a genus median of 6 for their 4, *B. fragilis* 5 for 6, *B. adolescentis*
  3 for 4.
- Even so, the masked table is lower than the read-level table on 51 of 52 samples
  (simulations −0.108, real mocks −0.038).

So a genus median brings most of the correction but can miss a well-known species by
two copies; the `copy_number_level` column of the cell table marks those rows, and the
leave-one-out errors of section 2 give their typical size.

## 7. Limits

- Copy numbers are per species; strains differ (*E. coli*'s complete genomes are 7 in
  3,948 of 4,006 cases, but some have 1 or up to 19), which `oracle` against `db`
  measures on the simulations (+0.014).
- The correction divides the read-level estimate: it cannot repair species the estimate
  misses or splits wrongly (near-neighbour pairs, D24/D39), and it moves their error
  with them.
- 24% of the default build's species have no copy number and are left out of the cell
  table; the comment line gives their read-level share. In the stool cohort that is
  0.8% of the assigned abundance, while 40% of the abundance rests on a genus or family
  median.
- The real-mock panel checks the copy numbers against Zymo's, not an independent cell
  count; a mock with measured cell counts (flow cytometry, qPCR) would test more.

## 8. Commands

```bash
B=~/pipelines/NanoVI_benchmarks
export NANOVI_REPO=$B/phase15/br_src   # git archive of NanoVI 92e3fbd
PY=$B/phase7a/env_py/bin/python
cd benchmarks
# Step 1: sources (downloads in $B/phase15/downloads)
$PY copy_number_sources.py compare && $PY copy_number_sources.py weighted
# The database with copy numbers
G=$B/phase7a/downloads/gtdb_r226
(cd $B/phase15/db && $PY $NANOVI_REPO/bin/build_gtdb_db.py --ssu $G/ssu_all_r226.fna.gz --db-name gtdb_r226 \
    --output-dir merged_cn --merge-gtdb-suffixes --gtdb-metadata $G/bac120_metadata_r226.tsv.gz $G/ar53_metadata_r226.tsv.gz)
$PY copy_number.py check-db
# Simulations (Badread 0.4.2 on PATH), the runs, the report's tables and figure
PATH=$B/phase7a/env_tools/bin:$PATH $PY copy_number.py simulate --threads 16
$PY copy_number.py run --workers 4
$PY copy_number.py analyse
```

The ported step 1 script reproduces the first run's `summary.json`, `weighted.json` and
species tables byte for byte.

NanoVI's review fixes after the benchmark (commit `8ae1a65`: blank or NA species names,
duplicated tax IDs, short metadata lines, tax IDs with no abundance, a species name
shared by several tax IDs counted once in a median) change none of these results: the
merged build's `copy_numbers.tsv` rebuilt with them is byte-identical, and so are the
species and cell tables of four samples rerun with them (two real, two simulated).
