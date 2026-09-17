# Layered evaluation numerical core

This core separates measured response retrieval from model prediction. It never
interprets low rank as absence of biological information, nor technical repeats
as donors. It is a research implementation, version 0.2.1. Source eligibility is
validated before fitting or scoring. The numerical kernels are unchanged from
version 0.2.0.

## Public Python interface

```python
from layered_eval import (
    align_rows, assert_disjoint, aggregate_responses, cosine_matrix,
    summarize_scores, score_rows, response_weights, cosine_from_grams,
    conditional_label_null, seeded_rng,
)
from layered_eval.bundle import validate_bundle
from layered_eval.predictions import score_predictions

# Always join models to explicit source cell IDs before aggregation.
y = align_rows(model_values, model_cell_ids, required_cell_ids)
responses = aggregate_responses(y, target_to_guide_to_indices,
                               control_guide_to_indices)
scores = cosine_matrix(query_responses, reference_responses)
summary = summarize_scores(scores)

# Explicit prediction/observed target and coordinate IDs are required.
# This evaluates supplied frozen arrays; it does not run a prediction model.
report = score_predictions('prediction_config.json', 'new_output_directory')
```

A score matrix is query target × candidate target with the identical target order
on both axes, or explicit `correct_columns` if using another candidate order.
Rows with ANY missing candidate score remain unavailable, avoiding a changing
candidate gallery. Zero directions have undefined cosine. Exact ties receive
mid-rank and fractional top-1; all tied K×K scores have rank 0.5 and top-1 1/K.
The conditional label null refuses any nonfinite matrix. It describes only the
fixed target panel and is not a biological confidence interval.

`aggregate_responses` weights guide centroids equally, rather than treating
large guides as more important. `response_weights` is an algebraically identical
sparse matrix operation. `cosine_from_grams` accelerates repeated sampling in high
dimensional RNA; squared norms <=1e-12 are treated as numerical zero. All actual
source/representation tasks were checked against direct response arrays.

## Portable bundle contract

- `manifest.json`: schema_version, cell_count and `representations`, each with a
  relative file path, type (`sparse_rna` or `array`) and dimensions.
- `cells.csv`: unique cell_id; guide_id; target_id; explicit boolean is_control
  and source_excluded_guide; original source_guide_n; source, experiment,
  background and known technical-unit metadata. Flags and source counts must be
  consistent within each guide. See SCHEMA.md for strict true/false parsing.
- `features.csv`: unique feature_id and original gene symbol/Ensembl provenance.
- `design.json`: ordered targets, target_groups[target][guide] = integer cell
  indices, control_groups[guide], pca_fit_indices and sampling settings.
- `rna.npz`: scipy CSR matrix on the explicit cell and feature axes.
- `representations/*.npy`: per-cell frozen model/PCA readouts in cell order.

Run `PYTHONPATH=src python -m layered_eval validate path/to/bundle` to check unique
IDs, target/group identity, representation shapes/values, source eligibility and
mutually disjoint PCA fitting, evaluated controls and treated cells. Excluded
guides cannot contribute to any of these used roles. This verifies computational
identity and design, not source annotation truth or biological validity.

Source-specific study adapters are distinct from this portable evaluator and are
distributed in the separate scientific archive. The portable commands and real
fixtures are documented in the package root README.md; no absent study scripts
are required to execute those fixtures.

## Verification

`PYTHONPATH=src python -m unittest layered_eval.test_core -v` checks ID reordering,
missing/duplicated IDs, ties, correctly remapped target labels, zero/NA directions,
invalid nulls, split leakage, guide-equal aggregation and exact dense/Gram agreement
including different-cell-source matrices.

`python -m unittest layered_eval.test_predictions -v` additionally checks actual
target/coordinate reordering, retained missing predictions, duplicate or mismatched
axes, zero and nonfinite directions, and incompatible output spaces. Both sets
total 13 tests. The two real-data fixtures are executed after installation outside
the scientific project; the prediction fixture preserves a previously evaluated
real 26-target matrix to numerical precision and exact target ranks.

The version 0.2.1 source-eligibility regression suite is in `tests/` at the package
root. Run `python -m unittest discover -s tests -v` using the installed version.
