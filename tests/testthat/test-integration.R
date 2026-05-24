# Integration tests -- require nonmem2rx and/or rxode2.
# All tests are gated with skip_if_not_installed(); they are skipped in CI
# unless the relevant packages are available.  On first run with the packages
# present, expect_snapshot() creates the _snaps/ files.  Subsequent runs
# compare against those snapshots.

# -- helpers ------------------------------------------------------------------

nm_path <- function(file) {
  system.file("testmodels", "nonmem", file, package = "ferxtranslate",
              mustWork = TRUE)
}

r2_path <- function(file) {
  system.file("testmodels", "nlmixr2", file, package = "ferxtranslate",
              mustWork = TRUE)
}

# -- NONMEM models ------------------------------------------------------------

test_that("1-cpt oral NONMEM: snapshot + no unsupported", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("1cpt_oral.ctl"))
  expect_snapshot(cat(result$ferx_text))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "one_cpt_oral", fixed = TRUE)
})

test_that("2-cpt oral with covariates: snapshot + no unsupported", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("2cpt_oral_cov.ctl"))
  expect_snapshot(cat(result$ferx_text))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "two_cpt_oral", fixed = TRUE)
})

test_that("2-cpt IV bolus: infers two_cpt_iv_bolus", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("2cpt_iv.ctl"))
  expect_snapshot(cat(result$ferx_text))
  expect_match(result$ferx_text, "two_cpt_iv_bolus", fixed = TRUE)
})

test_that("3-cpt IV: emits ERROR for unsupported three_cpt_iv_bolus", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("3cpt_iv.ctl"))
  expect_length(result$unsupported, 1L)
  expect_match(result$unsupported[1], "three_cpt", fixed = FALSE)
})

test_that("ODE warfarin: full $DES path, [odes] section present", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("ode_warfarin.ctl"))
  expect_snapshot(cat(result$ferx_text))
  expect_match(result$ferx_text, "[odes]",      fixed = TRUE)
  expect_match(result$ferx_text, "d/dt(depot)", fixed = TRUE)
  expect_length(result$unsupported, 0L)
})

test_that("block omega: block_omega line in output", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("block_omega.ctl"))
  expect_snapshot(cat(result$ferx_text))
  expect_match(result$ferx_text, "block_omega", fixed = TRUE)
})

test_that("IOV model: kappa line + iov_column in fit_options", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("iov.ctl"))
  expect_snapshot(cat(result$ferx_text))
  expect_match(result$ferx_text, "kappa",      fixed = TRUE)
  expect_match(result$ferx_text, "iov_column", fixed = TRUE)
})

# -- nlmixr2 models -----------------------------------------------------------

test_that("1-cpt oral nlmixr2: snapshot + one_cpt_oral", {
  skip_if_not_installed("rxode2")
  fn     <- source(r2_path("1cpt_oral_nlmixr2.R"))$value
  result <- nlmixr2_to_ferx(fn)
  expect_snapshot(cat(result$ferx_text))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "one_cpt_oral", fixed = TRUE)
})

test_that("ODE nlmixr2: d/dt expressions produce [odes] section", {
  skip_if_not_installed("rxode2")
  fn     <- source(r2_path("ode_nlmixr2.R"))$value
  result <- nlmixr2_to_ferx(fn)
  expect_snapshot(cat(result$ferx_text))
  expect_match(result$ferx_text, "[odes]",      fixed = TRUE)
  expect_match(result$ferx_text, "d/dt(depot)", fixed = TRUE)
  expect_length(result$unsupported, 0L)
})

# -- output writing -----------------------------------------------------------

test_that("nm_to_ferx writes file when output path given", {
  skip_if_not_installed("nonmem2rx")
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))
  nm_to_ferx(nm_path("1cpt_oral.ctl"), output = path)
  expect_true(file.exists(path))
  expect_match(paste(readLines(path), collapse = "\n"), "[parameters]", fixed = TRUE)
})
