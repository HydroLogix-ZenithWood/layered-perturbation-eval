# Submission software copy: source attribution and licensing status

## Included public-data examples

The examples in this companion derive from the Papalexi THP-1 ECCITE-seq source, Gene Expression Omnibus accession [GSE153056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153056), RNA sample [GSM4633614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4633614). The local partition is `papalexi2021_rep1`, one of three independently transduced partitions in the same THP-1 background. These are not participant records or three independent studies.

`examples/native_role_contract` contains the actual source's frozen cell-role metadata, the declared ten-target contract, and token sequences only for the selected native-start cells. There are 39 unique start cells and 160 cell-by-target uses in this first partition; a cell may be reused across different target operations. Declared source and model identities are preserved.

When included after successful inference and numerical qualification, `examples/external_native_fixture` contains ten measured response vectors and ten vectors for each of target deletion and neighbouring-token deletion, all in the same 768-dimensional Geneformer V2 CLS space, with target/feature IDs, expected score matrices and provenance hashes. These are derived coordinates, not gene-expression measurements. The fixture uses the metadata-selected first draw; it is not selected for favourable prediction performance and supplies no new biological replication.

Model revision: `ctheodoris/Geneformer`, commit `1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5`, available from the [official repository](https://huggingface.co/ctheodoris/Geneformer/tree/1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5). Original Geneformer code, model and input-data terms continue to apply. No checkpoint or model weight is included. Neighbouring-token deletion is a computational comparison, not a known biologically inactive intervention.

Original investigators and source repositories must be attributed. Public availability is not represented as a blanket licence to redistribute every input or derivative. This research companion grants no new rights over source data, model assets or their derivatives. Applicable source-data and model terms continue to govern reuse and redistribution.

## Numerical software and original core

Version 0.2.1 repairs the source-eligibility preparation and bundle validation of
the previous v0.2.0 package. Its numerical scoring kernels and separate standard-
library native-input validator remain unchanged. Tested numerical dependency
versions are recorded in `requirements-tested.txt`: NumPy 1.26.4, SciPy 1.11.4,
pandas 2.1.4, scikit-learn 1.3.2 and threadpoolctl 3.6.0. Dependencies are not
vendored; their original copyright and licence notices govern their distributions.
The input validator alone does not require these numerical libraries.

The unchanged real measured-RNA and WTC11 native examples from the previous v1
package are now included alongside the THP-1 examples. The count fixture derives
from [GSE335887](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335887),
`cite_6tf`, and preserves observed counts, source labels and full-feature totals.
The WTC11 native fixture derives from
[GSE178317](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178317),
`microglia_pool`, and existing qualified Geneformer V2 inference. It contains 26
derived 768-dimensional predicted and observed-response vectors. Neither fixture
adds biological replication or includes model weights. The source-data and dependency provisions of the original notices are retained
in `provenance/V1_THIRD_PARTY_NOTICES.md` and
`provenance/V2_THIRD_PARTY_NOTICES.md`; only their original-code licensing status
and release wording are updated in this release copy. Public availability does not grant a
blanket redistribution licence; each source's attribution and applicable terms
remain in force.

## New implementation

The original evaluation code is provided under the top-level MIT LICENSE, with copyright attributed to 2026 layered-perturbation-eval contributors. Third-party, source-data, derived-array and model terms retain their own scope. The scReliability author function belongs only to the separate scientific evidence archive, which preserves its MIT licence and fixed commit; it is not a runtime dependency or vendored function in this companion.

This package supports input-contract validation and evaluation of supplied arrays. It does not independently rerun all model inference, demonstrate unseen pretraining membership, imply journal submission or establish manuscript-specific author declarations.
