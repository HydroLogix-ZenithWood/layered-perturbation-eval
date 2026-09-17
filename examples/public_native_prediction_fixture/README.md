# Real native prediction execution fixture

This contains actual frozen Geneformer V2 CLS target-deletion response vectors and actual measured microglial response vectors from the existing public study, with every eligible target and the first frozen A-arm draw. It is an execution demonstration, not an additional independent biological validation.

Every vector uses the same 768-coordinate model output basis; these coordinates are latent model dimensions, not genes. Config provenance records original model-array and draw hashes. Native starting cells vary by target and overlap the observed reference to the extent explicitly recorded in `native_start_coverage.csv`. Counts of native starting cells are not matched to the observed response budget.

Run `layered-eval score-predictions examples/public_native_prediction_fixture/config.json --output /tmp/native-report`. This scores supplied frozen outputs and does not load weights or execute native perturbation inference.
