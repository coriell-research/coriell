## Submission

This is a new submission.

## Test environments

* local Ubuntu 24.04 (Pop!_OS), R 4.6.2
* win-builder, R-devel (2026-11-01 r90825 ucrt)
* win-builder, R-release -- TODO
* macOS builder (release) -- TODO
* R-hub (linux, windows, macos, nosuggests) -- TODO

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes for reviewers

* Several packages in Suggests are distributed through Bioconductor
  (csaw, DelayedArray, DelayedMatrixStats, edgeR, fgsea, limma, metapod,
  PCAtools, qsmooth, SummarizedExperiment). All are used conditionally:
  functions check `requireNamespace()`, examples use `\examplesIf`, and
  tests use `skip_if_not_installed()`.
* `panther_go()` queries the PANTHER web service. Its example only runs in
  interactive sessions, the tests mock all HTTP responses (the one live test
  is skipped on CRAN), and the function returns `NULL` with a message when
  the service is unavailable.
