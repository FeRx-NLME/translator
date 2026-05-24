# Numerical concordance tests: verify that translated .ferx models fit data and
# recover known true parameters within tolerance.
#
# These tests are gated on ferx being installed and are skipped on CI (too slow,
# require the ferx binary in PATH). Run locally with devtools::test().
#
# Bundled datasets in inst/testdata/ were simulated from the translated models
# using ferx_simulate() at the theta initial values (which are the "true" values
# for the simulation). See data-raw/generate_concordance_data.R for the script.
#
# Acceptance criteria (per plans/v0.1-implementation.md Section 10):
#   structural thetas: within 15% of truth
#   random-effect variances: within 20% of truth

skip_if_not_installed("ferx")
skip_on_ci()

library(ferxtranslate)
library(ferx)

# Helper: translate a bundled NONMEM model and write .ferx to a temp file
.translate_to_tmp <- function(model_name) {
  ctl <- system.file(file.path("testmodels/nonmem", model_name),
                     package = "ferxtranslate")
  result <- nm_to_ferx(ctl)
  ferx_file <- tempfile(fileext = ".ferx")
  writeLines(result$ferx_text, ferx_file)
  ferx_file
}

# Helper: bundled concordance dataset path
.conc_data <- function(name) {
  system.file(file.path("testdata", name), package = "ferxtranslate")
}

# ---------------------------------------------------------------------------
# 1-cpt oral (linCmt → one_cpt_oral pk macro)
#   True thetas: TVCL=0.134, TVV=8.1, TVKA=1.0
#   Simulated with 100 subjects; large omega_KA (0.4) means KA needs wider
#   tolerance -- TVCL and TVV are the more sensitive translation targets.
# ---------------------------------------------------------------------------

test_that("1-cpt oral: TVCL and TVV recover within 15% of truth", {
  ferx_file <- .translate_to_tmp("1cpt_oral.ctl")
  data_file  <- .conc_data("1cpt_oral_concordance.csv")
  fit <- ferx_fit(ferx_file, data_file,
                  method     = "focei",
                  covariance = FALSE,
                  verbose    = FALSE)

  # True values from theta initials in 1cpt_oral.ctl
  expect_lt(abs(fit$theta["TVCL"] / 0.134 - 1), 0.15,
            label = "TVCL relative error")
  expect_lt(abs(fit$theta["TVV"]  / 8.1   - 1), 0.15,
            label = "TVV relative error")
})


# ---------------------------------------------------------------------------
# 2-cpt IV bolus (linCmt → two_cpt_iv_bolus pk macro)
#   True thetas: TVCL=5.0, TVV1=20.0, TVQ=8.0, TVV2=60.0
#   Small IIV (omega <= 0.10); all four PK params expected within 10%.
# ---------------------------------------------------------------------------

test_that("2-cpt IV: all structural thetas recover within 10% of truth", {
  ferx_file <- .translate_to_tmp("2cpt_iv.ctl")
  data_file  <- .conc_data("2cpt_iv_concordance.csv")
  fit <- ferx_fit(ferx_file, data_file,
                  method     = "focei",
                  covariance = FALSE,
                  verbose    = FALSE)

  # True values from theta initials in 2cpt_iv.ctl
  expect_lt(abs(fit$theta["TVCL"] / 5.0  - 1), 0.10, label = "TVCL")
  expect_lt(abs(fit$theta["TVV1"] / 20.0 - 1), 0.10, label = "TVV1")
  expect_lt(abs(fit$theta["TVQ"]  / 8.0  - 1), 0.10, label = "TVQ")
  expect_lt(abs(fit$theta["TVV2"] / 60.0 - 1), 0.10, label = "TVV2")
})

test_that("2-cpt IV: omega estimates recover within 10% of truth", {
  ferx_file <- .translate_to_tmp("2cpt_iv.ctl")
  data_file  <- .conc_data("2cpt_iv_concordance.csv")
  fit <- ferx_fit(ferx_file, data_file,
                  method     = "focei",
                  covariance = FALSE,
                  verbose    = FALSE)

  # 2-cpt IV omegas are larger (0.05-0.10) and recover within 1% in practice;
  # 10% tolerance gives a comfortable margin while still catching translation bugs.
  expect_lt(abs(fit$omega["ETA_CL", "ETA_CL"] / 0.10 - 1), 0.10, label = "omega_CL")
  expect_lt(abs(fit$omega["ETA_V1", "ETA_V1"] / 0.10 - 1), 0.10, label = "omega_V1")
  expect_lt(abs(fit$omega["ETA_Q",  "ETA_Q"]  / 0.08 - 1), 0.10, label = "omega_Q")
})
