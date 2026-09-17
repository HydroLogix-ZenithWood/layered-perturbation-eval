# Third-party provenance and licensing status

Release-copy note: source and dependency provisions below preserve the historical
package scope. The original-code licensing and release-status wording has been
updated for this distribution; the sealed source notices are unchanged.

## Demonstration data

The small example derives from public Gene Expression Omnibus accession
GSE335887, experiment `cite_6tf`:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335887

It preserves observed raw counts for an explicitly selected subset, the original
full-feature cell totals, guide/target labels and source-selection metadata.
Original paper/data-source attribution is required. Public accessibility is not
represented as a blanket redistribution license. No new license is assigned to
the example data; applicable source terms continue to govern redistribution.
This distribution contains no complete third-party model, checkpoint or original
full experiment archive.

The native-prediction example derives from GSE178317 (`microglia_pool` WTC11),
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178317, and the qualified
Geneformer V2-104M frozen inference already present in the scientific project.
It contains 26 derived 768-dimensional predicted vectors and 26 measured-response
vectors, target/coordinate IDs and provenance hashes. These latent coordinates
are not gene-expression measurements. Its purpose is executable evaluation of
real supplied predictions, not redistribution of the foundation-model weights or
a new independent biological result. Geneformer code/model and source-data terms
remain applicable to reuse and redistribution of derived data.

## Numerical software

Direct runtime dependencies and locally verified versions are recorded in
`requirements-tested.txt`: NumPy 1.26.4, SciPy 1.11.4, pandas 2.1.4, scikit-learn
1.3.2 and threadpoolctl 3.6.0. Installed package metadata identifies the project
licenses as BSD-style/BSD 3-Clause. Dependencies and their full copyright/license
notices are not vendored here and remain governed by their own distributions.
Package building uses setuptools and wheel under their respective licenses.

## New implementation

The original evaluation code is provided under the top-level MIT LICENSE, with copyright attributed to 2026 layered-perturbation-eval contributors. Third-party, source-data, derived-array and model terms retain their own scope.
