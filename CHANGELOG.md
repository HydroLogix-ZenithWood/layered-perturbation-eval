# Changes in version 0.2.1

This patch repairs source-eligibility handling. The frozen scientific results
are unchanged.

- `prepare_counts` excludes flagged control guides before allocating PCA and
  reference pools and applies `min_source_cells_per_guide` to both control and
  treated guides. If insufficient clean controls remain, it rejects preparation
  before writing an output bundle.
- Preparation and validation share strict metadata checks. Both boolean fields
  must be explicit true/false, with case and outer whitespace normalized. Missing
  columns/values, numeric flags, other spellings and inconsistent guide-level
  flags or original source counts produce actionable errors.
- `validate_bundle` independently rejects excluded or below-threshold guides in
  PCA-fitting, reference and treated roles. New designs save the threshold;
  legacy designs use the documented default 20 and label that fallback explicitly.
- Real v1 count/native examples are included unchanged alongside the external
  THP-1 examples. Both original third-party notices and historical execution
  provenance are retained. The distribution adds the MIT licence for original code; source-data and third-party terms
  remain unchanged.
- Regression tests reproduce the excluded-control failure and test the corrected
  boundary behavior. Test-only altered metadata and the duplicated extra guide
  are not scientific evidence or new biological samples.

The core cosine/rank/tie/NA kernels, random sampling, normalization formula and
supplied-native-vector scoring code are unchanged. The patch is not a model
update or a new scientific analysis. The separate release audit records fresh
environment execution, fixture equivalence and the impact check on the eight
frozen scientific bundles. Earlier checks in `provenance/` apply only to their
historical versions and must not be reported as the version 0.2.1 acceptance run.
