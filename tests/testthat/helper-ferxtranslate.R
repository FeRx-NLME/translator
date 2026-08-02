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
