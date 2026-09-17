# Revised manuscript: supplied-table checks

Run from the repository root with Python 3.10 or later:

```sh
python manuscript-reproducibility/verify_derived_tables.py
```

Only Python's standard library is needed. Paths are relative to the script. The
script does not alter the supplied tables. An optional `--output NEW_REPORT.json`
writes a new report and refuses to overwrite an existing one.

## Figure 4 and Table S24

Panels A–C describe statistics calculated after averaging ranks for each target
over 100 draws. Panel D summarizes paired differences calculated separately
within each of those draws. Nonlinear aggregation means that these are different
estimands. The script independently calculates rank correlations and fractional
tie-aware coverage curves from the supplied 17-target mean-rank tables. It also
recalculates means and 2.5th–97.5th percentile ranges from the supplied 100-draw
metric tables and checks every condition against the estimand summary.

The range in panel D describes the conditional draw distribution. It is not a
confidence interval for either the mean difference or the target-mean statistic.
Table S24 lists both aggregation orders. All six conditional ranges for the
guide-aware minus technical correlation comparison include zero; the ordering
check does not establish a stable guide-aware gain in these panels.

The original 100 × 17 × 17 score matrices are not included in this compact
payload. This script verifies aggregation of supplied derived outputs; it does
not rebuild the ranks from original score matrices or rerun inference.

## Cell-ID audit tables: current Supplementary Methods S11

Tables S25–S28 retain the distinctions between earlier query-specific exclusion,
later matched-outcome reference sampling, and external globally separated cell
roles. The current section is S11 after removal of the former S10; table numbers
are unchanged.

- Table S25: the script recalculates counts, proportions, ranges and means from
  the 79,200 supplied target × direction × draw overlap records. It does not
  recalculate the union of unique starting-cell IDs.
- Table S26: checks the supplied earlier query-specific exclusion records for
  zero start/reference and start/treated overlap and membership in the original
  control pool. This design excludes a query's own starts; it is not a global
  exclusion of every target's starting cells.
- Table S27: checks the supplied external role-overlap records and native budgets.
- Table S28: checks that common and version-specific gallery sizes agree.

These are checks of reported intersection counts and table structure. Fresh
intersection calculations require the frozen original cell metadata, native-start
indices, reference-draw identities and gallery definitions; they are not
reconstructed from these tables. Full raw-data processing, tokenization and model
inference require the original public datasets, the specified model checkpoints
and the full study input layout. The compact interface fixtures in `examples/`
do not replace those inputs.

`input_manifest.json` records the hashes of the ten supplied CSVs. They retain
their corresponding source-data provenance and terms; the repository code
licence does not relicense all input or model-derived assets. See the repository's
`THIRD_PARTY_NOTICES.md` for GSE335887, GSE178317 and GSE153056 attribution. The
overlap tables also include AG07657 public microglial data from
[GSE311359](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311359).
