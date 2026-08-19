# Shared test helpers. testthat sources helper-*.R before every test file, so
# these live here rather than being duplicated per file (nm_path / nm_path_t).

nm_path <- function(file) {
  system.file("testmodels", "nonmem", file, package = "ferxtranslate",
              mustWork = TRUE)
}

r2_path <- function(file) {
  system.file("testmodels", "nlmixr2", file, package = "ferxtranslate",
              mustWork = TRUE)
}

# Every bundled NONMEM model. Filtered by extension so a README, dataset or
# .lst dropped into the directory does not get fed to nonmem2rx().
.bundled_nm_models <- function() {
  list.files(system.file("testmodels", "nonmem", package = "ferxtranslate"),
             pattern = "\\.(ctl|mod)$", full.names = TRUE)
}

# Base-R scratch directory; withr is not a declared dependency of this package,
# and adding one for a two-line helper would put a WARNING in
# `R CMD check --as-cran` ("'::' or ':::' import not declared"). R clears
# tempdir() at session end, so no explicit cleanup is needed.
#
# Lives here rather than in one test file because more than one file needs it:
# defined inside test-translate.R it was invisible to test-rxui_to_ir.R, which
# is the kind of breakage that shows up as an opaque "could not find function".
tmp_ctl_dir <- function() {
  d <- tempfile("ferxtr-")
  dir.create(d)
  d
}
