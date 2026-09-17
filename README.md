# Layered perturbation evaluation 0.2.1

This research toolkit combines the repaired source-eligibility interface
with real count and supplied-prediction examples and a standard-library native
input-contract validator. Version 0.2.1 enforces excluded-guide flags and original
source-count thresholds for controls as well as treated guides, and rejects
missing, malformed or inconsistent boolean metadata. The scoring kernels,
normalization formula, random draws and native prediction scorer are unchanged.
See SCHEMA.md and CHANGELOG.md for the exact contract and verification boundary.

It is an execution and identity-check demonstration, not an additional biological
validation or a universal native-prediction engine. The `manuscript-reproducibility` directory provides portable checks of the
revised manuscript's supplied derived tables. Full study matrices and model
inference adapters are separate. The original evaluation
code is provided under the MIT licence in LICENSE. Source data, derived arrays and models retain their applicable
terms; see THIRD_PARTY_NOTICES.md.

## Install in an isolated environment

Python 3.10 or later is required. From this extracted package directory:

```sh
python3.10 -m venv /tmp/layered-eval-clean
/tmp/layered-eval-clean/bin/python -m pip install -r requirements-tested.txt
/tmp/layered-eval-clean/bin/python -m pip install --no-deps .
```

Use the new environment's Python or `layered-eval` entry point. The examples need
no model weights, GPU, credentials or network once dependencies are installed.
Every output directory below must be new; the interface refuses overwrites.

## Run the real count example

```sh
/tmp/layered-eval-clean/bin/layered-eval prepare examples/public_microglia_fixture/config.json --output /tmp/le-counts-bundle
/tmp/layered-eval-clean/bin/layered-eval validate /tmp/le-counts-bundle
/tmp/layered-eval-clean/bin/layered-eval evaluate /tmp/le-counts-bundle --output /tmp/le-counts-result --iterations 4 --n 16
```

This is the unchanged GSE335887 fixture with 464 real cells, four targets, four
control guides and 512 nonpanel coordinates. Its original full-feature library
totals are used for normalization. It is an execution subset, not an additional
biological validation. Eligibility flags may be retained for unused cells, but
excluded guides cannot contribute to PCA fitting or evaluated responses.

The measured evaluator materializes a cell-by-cell Gram matrix, so its memory
cost grows quadratically with source cell count. The scientific sources contained
approximately 1,800–2,700 cells per bundle. Atlas-scale execution is not established.

## Run the existing WTC11 native-vector example

```sh
/tmp/layered-eval-clean/bin/layered-eval score-predictions examples/public_native_prediction_fixture/config.json --output /tmp/le-wtc11-native
```

The unchanged GSE178317 fixture contains 26 targets in 768-dimensional Geneformer
V2 CLS space. Its historical start budgets vary and reported start/reference
overlap is retained. Do not label it as satisfying the stronger THP-1 role contract.

## Validate a native task's input contract

From this directory, with Python 3.10 or later:

```sh
python validate_native_contract.py examples/native_role_contract/config.json --output /tmp/native-contract-result.json
```

The command checks source identity against metadata, exact declared target and cell identities, native guide/cell budgets, declared nonspecial target-token visibility, and disjoint native/reference/PCA roles. Metadata row order may change; identity must not. It rejects duplicates, missing targets, special or invisible target tokens, unequal start budgets, and reference overlap. It does not authenticate a token's mapping to a biological gene: that requires the separately frozen official dictionary audit. It cannot infer biological validity, experimental independence beyond declared metadata, model training exposure, or whether token deletion is a valid physical perturbation.

The example uses the actual first external source's preselected 10 targets, two control guides and eight cells per guide per target. Only sequences needed for native starts are included; all cell roles are present. These input-contract facts were frozen before new native response scoring. Original source and model terms remain applicable; packaging does not assign new rights.

## Score supplied response vectors

The `src/layered_eval` package retains `prepare`, `validate`, `evaluate`, and `score-predictions`; see SCHEMA.md for contracts. The qualified external native response example is included in examples/external_native_fixture after completed inference and numerical audits. It supplies target, neighbour-operation, and observed vectors from a metadata-selected first draw. No weights, GPU, credentials, or network are needed to score already supplied arrays. The measured and predicted coordinate axes must match exactly, and missing predictions never shrink the candidate gallery.

For this numerical core only, the previously tested dependencies are listed in requirements-tested.txt. No environment upgrade is required if those dependencies are already available. To run directly set PYTHONPATH=src and call `python -m layered_eval score-predictions CONFIG --output NEW_DIRECTORY`. Refusing an existing output directory preserves earlier results.

## Scope of reuse

These examples demonstrate the stated interfaces on a real new source. They do not reproduce the full raw-data-to-model-inference pipeline on another machine, establish independent human replication, or guarantee arbitrary model compatibility. The original evaluation code is covered by the MIT licence in LICENSE. See RELEASE_STATUS.md for the current distribution status.

## Reproduce the new external matrices

```sh
PYTHONPATH=src python -m layered_eval score-predictions examples/external_native_fixture/target_config.json --output /tmp/external-target-score
PYTHONPATH=src python -m layered_eval score-predictions examples/external_native_fixture/neighbour_config.json --output /tmp/external-neighbour-score
```

Each command keeps the same ten-target gallery and 768-dimensional qualified CLS output space. Expected matrices and source hashes are provided in the example. The candidate neighbours are computational controls, not known biological negative interventions.

## Run numerical and source-eligibility checks

From the extracted package directory, using the installed package in the clean
environment (do not set PYTHONPATH for this installed-package acceptance check):

```sh
/tmp/layered-eval-clean/bin/python -m unittest layered_eval.test_core layered_eval.test_predictions -v
/tmp/layered-eval-clean/bin/python -m unittest discover -s tests -v
```

The original 13 numerical checks are preserved. The new regression suite checks
excluded controls, excluded treated guides, source-count eligibility, each used
role in a tampered bundle, explicit boolean parsing, within-guide consistency and
unchanged results when an additional excluded guide is present. Altered metadata
and a cloned extra guide occur only in temporary test inputs, never as added
scientific data. Historical v0.2.0 execution records are retained under provenance;
they do not certify the repaired version.

## Check the revised manuscript tables

```sh
python manuscript-reproducibility/verify_derived_tables.py
```

This standard-library check distinguishes the two Figure 4 estimands, checks
their supplied summaries, and verifies the internal consistency of the
version-specific overlap tables. It does not rerun model inference or reconstruct
cell-ID intersections from raw metadata. See `manuscript-reproducibility/README.md`.
