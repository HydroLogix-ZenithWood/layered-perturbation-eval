# Input schema

All file locations in a counts configuration are resolved relative to that
configuration. Neither the generic module nor fixture requires absolute paths.

| File/key | Contract |
|---|---|
| counts | scipy sparse CSR NPZ; rows match cells.csv, columns match features.csv; finite nonnegative integer counts |
| cells | CSV with unique nonmissing cell_id; nonmissing guide_id and target_id; explicit boolean is_control and source_excluded_guide; positive integer source_guide_n; original full-feature library total column |
| features | CSV with unique feature_id and gene_symbol; original ensembl_id can also be retained |
| full_count_total_column | Name of the original library-total column; values must be positive and at least subset row sums |
| dataset_id | Stable input identifier used in deterministic sampling |
| source_panel | All source perturbation gene symbols excluded from RNA coordinates, not only selected targets |
| target_ids | Frozen ordered candidate list selected without effect scores |
| min_source_cells_per_guide | Default 20, eligibility of the original labelled guide population |
| min_sampled_cells_per_guide | Default 32, fixed pool supporting two disjoint n=16 arms |
| pca_fit_cells_per_control_guide | Default 20; permanently reserved controls |
| reference_cells_per_control_guide | Default 16; two control identities per arm |
| pca_dimensions | Default 64, capped by fit sample/coordinate rank; use 0 to omit PCA |
| coordinate_scope | Human-readable full/partial coordinate provenance; example is a restricted 512-coordinate fixture |
| provenance | Public accession, source experiment, hashes, selection rule and redistribution notes |

## Source eligibility in version 0.2.1

`is_control` and `source_excluded_guide` accept only CSV `true`/`false`, with
case ignored and optional outer whitespace. `True`/`False` produced by pandas
are accepted. Blank/missing values, `0`/`1`, `yes`/`no` and other spellings fail
with the column and an affected cell/data-row location. Both columns are required
for preparation and validation, including on excluded or otherwise unused cells.
There is no implicit false default or string-truthiness conversion.

Both flags and the positive integer `source_guide_n` must be constant within each
`guide_id`. Resolve inconsistent source annotations before execution. A guide
marked `source_excluded_guide=true` can remain in the preserved metadata and count
array but cannot contribute any PCA-fitting, reference or treated cell. The same
`min_source_cells_per_guide` threshold applies to control and treated guides.
The original source count is distinct from the smaller sampled cell count.

`prepare` rejects insufficient eligible control-guide capacity after applying the
source flags, source-count threshold and fitting/reference budgets. It saves the
threshold in `design.json`. `validate` independently rejects excluded or below-
threshold guides in any of the three used roles, even if the bundle was prepared
elsewhere or metadata was changed later. A legacy bundle that omits the threshold
is checked at the documented default 20; the validation output explicitly records
`source_min_cells_origin=legacy_default_20`. This fallback must not be used to
claim a different historical custom threshold was verified.

The target qualification behavior is unchanged. `prepare` can omit a requested
target with fewer than two eligible guides, records its eligibility, and writes
the qualified gallery to `design.json`. Reconcile the requested and resulting
lists before freezing a task. An unexpected change to an already frozen gallery
is a stop condition, not permission to improve scores by narrowing it.

Prepared bundle contract: `manifest.json`, `design.json`, `cells.csv`,
`features.csv`, `rna.npz`, optional `representations/*.npy`, independent PCA fit
records and eligibility ledger. A supplied model representation must be joined to
explicit cell IDs with `align_rows` before saving it in bundle order; do not infer
matching identity from equal row counts. Prediction arrays require additional
model-native intervention and starting-cell records and are not implied by a
representation artifact.

The local scientific study additionally preserves full metadata and complete
source hashes for five original sources and three external transduction replicates.
The example's source subsets retain full-cell totals, exact source
cell row indices and original gene symbols. No fake counts or generated cells are
inserted into this example.

## Frozen prediction configuration

`score-predictions` accepts a JSON configuration with `dataset_id`,
`gallery_targets`, `predictions`, `observed` and `evaluation_provenance`.
`gallery_targets` is an ordered JSON list of unique string target IDs or a path
to such a JSON file. It is the complete frozen candidate gallery, not a result
of filtering by model score.

Both `predictions` and `observed` contain `values` (2D NPY array), `targets` and
`features` (unique string-ID lists or relative JSON file paths), `output_space`,
`readout` and nonempty `provenance`. Rows are targets and columns are coordinates.
The two output-space/readout strings and the two coordinate ID sets must match
exactly. Coordinate ordering can differ and is explicitly aligned. Observed
target rows must cover the exact gallery. Prediction target rows can be a subset;
missing predictions remain NA query rows in the full gallery. An empty prediction
target list paired with a zero-row
array is valid and explicitly records complete unavailability. Extra prediction
targets, duplicate IDs, mismatched axes or incompatible spaces fail validation
before any output is created. Nonfinite or zero directions remain unavailable.

Output includes `score_matrix.npz`, `per_target.csv`, `input_provenance.json`,
`report.json` and `report.md`. Input hashes, target/coordinate alignment and the
fixed gallery are saved. Cosine correspondence is evaluated in the declared space;
an embedding delta is not interchangeable with an RNA delta. The command does
not verify a model's unprovided provenance, execute inference, pool separate
models into a common leaderboard or establish biological causality.
