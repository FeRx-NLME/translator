# Numerical concordance tests: verify that translated .ferx models fit data and
# recover known true parameters within tolerance.
#
# These tests are gated on ferx being installed and are skipped on CI (too slow,
# require the ferx binary in PATH). Run locally with devtools::test().
#
# Two sets of tests:
#
# Self-contained (always run when ferx is installed):
#   Datasets in inst/testdata/*_concordance.csv were simulated from the
#   translated models using ferx_simulate() at the theta initial values.
#   See data-raw/generate_concordance_data.R for the generation script.
#
# amp.sim (also gated on amp.sim being installed):
#   Truth values come from published NONMEM reference estimates stored in
#   amp.sim's PK.1CMT.ORAL.ext file. The bundled dataset was simulated at
#   those reference values via ferx_simulate(). This tests the linCmt
#   translation path (one_cpt_oral pk macro) against an external NONMEM
#   reference rather than our own model initials.
#   Note: the amp.sim ODE path (pk_1cmt_oral.mod / ADVAN6) is not tested
#   here because the S2=V compartment scaling is not yet translated.
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

# Helper: print a deviation table and return pct errors invisibly.
# estimated: named numeric (fit$theta or fit$omega diagonal)
# reference: named numeric of true/reference values with matching names
# label:     printed above the table (e.g. model name)
.report_deviations <- function(estimated, reference, label = "") {
  if (nzchar(label)) message("\n-- ", label, " --")
  pct <- vapply(names(reference), function(nm) {
    est <- if (nm %in% names(estimated)) estimated[[nm]] else NA_real_
    (est / reference[[nm]] - 1) * 100
  }, numeric(1))
  rows <- data.frame(
    param    = names(reference),
    reference = unname(reference),
    estimated = unname(estimated[names(reference)]),
    pct_error = round(pct, 1),
    stringsAsFactors = FALSE
  )
  message(paste(capture.output(print(rows, row.names = FALSE)), collapse = "\n"))
  invisible(pct)
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

  ref <- c(TVCL = 0.134, TVV = 8.1)
  .report_deviations(fit$theta, ref, "1-cpt oral thetas")
  expect_lt(abs(fit$theta["TVCL"] / ref["TVCL"] - 1), 0.15, label = "TVCL relative error")
  expect_lt(abs(fit$theta["TVV"]  / ref["TVV"]  - 1), 0.15, label = "TVV relative error")
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

  ref <- c(TVCL = 5.0, TVV1 = 20.0, TVQ = 8.0, TVV2 = 60.0)
  .report_deviations(fit$theta, ref, "2-cpt IV thetas")
  expect_lt(abs(fit$theta["TVCL"] / ref["TVCL"] - 1), 0.10, label = "TVCL")
  expect_lt(abs(fit$theta["TVV1"] / ref["TVV1"] - 1), 0.10, label = "TVV1")
  expect_lt(abs(fit$theta["TVQ"]  / ref["TVQ"]  - 1), 0.10, label = "TVQ")
  expect_lt(abs(fit$theta["TVV2"] / ref["TVV2"] - 1), 0.10, label = "TVV2")
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
  omega_diag <- c(ETA_CL = fit$omega["ETA_CL", "ETA_CL"],
                  ETA_V1 = fit$omega["ETA_V1", "ETA_V1"],
                  ETA_Q  = fit$omega["ETA_Q",  "ETA_Q"])
  ref_omega  <- c(ETA_CL = 0.10, ETA_V1 = 0.10, ETA_Q = 0.08)
  .report_deviations(omega_diag, ref_omega, "2-cpt IV omegas")
  expect_lt(abs(omega_diag["ETA_CL"] / ref_omega["ETA_CL"] - 1), 0.10, label = "omega_CL")
  expect_lt(abs(omega_diag["ETA_V1"] / ref_omega["ETA_V1"] - 1), 0.10, label = "omega_V1")
  expect_lt(abs(omega_diag["ETA_Q"]  / ref_omega["ETA_Q"]  - 1), 0.10, label = "omega_Q")
})

# ---------------------------------------------------------------------------
# amp.sim benchmark: compare against published NONMEM reference estimates
#   Model : pk_1cmt_oral_ampsim.ctl (ADVAN2 linCmt; mirrors amp.sim
#           PK.1CMT.ORAL IIV structure: ETA on KA and CL only, V fixed)
#   Truth : amp.sim PK.1CMT.ORAL.ext final estimates (NONMEM FOCEI run)
#   Data  : simulated at reference parameter values via ferx_simulate()
#           (NM.theoph.02B.csv is not bundled in the amp.sim package)
#
# This validates the one_cpt_oral pk macro translation against an external
# NONMEM reference, not just against our own model initials.
# ---------------------------------------------------------------------------

test_that("amp.sim: 1-cpt oral thetas recover within 10% of NONMEM reference", {
  skip_if_not_installed("amp.sim")

  # Load published NONMEM reference estimates from amp.sim
  ext_file <- system.file("example_models/PK.1CMT.ORAL.ext", package = "amp.sim")
  ext      <- read.table(ext_file, header = TRUE, skip = 1)
  ref      <- ext[ext$ITERATION == -1000000000, ]

  ferx_file <- .translate_to_tmp("pk_1cmt_oral_ampsim.ctl")
  data_file  <- .conc_data("ampsim_1cpt_oral_concordance.csv")
  fit <- ferx_fit(ferx_file, data_file,
                  method     = "focei",
                  covariance = FALSE,
                  verbose    = FALSE)

  ref_nm <- c(KA = ref$THETA1, CL = ref$THETA2, V = ref$THETA3)
  .report_deviations(fit$theta, ref_nm, "amp.sim 1-cpt oral thetas vs NONMEM reference")
  expect_lt(abs(fit$theta["KA"] / ref_nm["KA"] - 1), 0.10, label = "KA vs amp.sim ref")
  expect_lt(abs(fit$theta["CL"] / ref_nm["CL"] - 1), 0.10, label = "CL vs amp.sim ref")
  expect_lt(abs(fit$theta["V"]  / ref_nm["V"]  - 1), 0.10, label = "V vs amp.sim ref")
})

# ---------------------------------------------------------------------------
# ODE path: pk_1cmt_oral.mod (ADVAN6 with S2=V scaling)
#   True thetas: KA=0.1, CL=2.0, V=1.0 (theta initials)
#   Tests that [scaling] obs_scale=V divides amount by V before comparing to
#   concentration data -- without scaling, IPRED >> DV and fit diverges.
# ---------------------------------------------------------------------------

test_that("ODE 1-cpt oral with S2=V: structural thetas recover within 15% of truth", {
  ferx_file <- .translate_to_tmp("pk_1cmt_oral.mod")
  data_file  <- .conc_data("ode_1cpt_oral_concordance.csv")
  fit <- ferx_fit(ferx_file, data_file,
                  method     = "focei",
                  covariance = FALSE,
                  verbose    = FALSE)

  ref <- c(KA = 0.1, CL = 2.0, V = 1.0)
  .report_deviations(fit$theta, ref, "ODE 1-cpt oral thetas")
  expect_lt(abs(fit$theta["KA"] / ref["KA"] - 1), 0.15, label = "KA")
  expect_lt(abs(fit$theta["CL"] / ref["CL"] - 1), 0.15, label = "CL")
  expect_lt(abs(fit$theta["V"]  / ref["V"]  - 1), 0.15, label = "V")
})

# ===========================================================================
# FEATURE-GAP SKELETONS
# Each block below corresponds to a ferx-core feature that is not yet
# supported. The test is skipped with a labelled message so:
#   - The skip list in CI output documents the outstanding gaps.
#   - When ferx-core ships the feature, remove the skip(), add a test
#     model to inst/testmodels/nonmem/, generate concordance data with
#     data-raw/generate_concordance_data.R, and wire up the assertions.
# ===========================================================================

# ---------------------------------------------------------------------------
# 3-cpt infusion (ADVAN12 TRANS4 with zero-order input)
#   ferx-core gap: three_cpt_infusion pk macro not yet implemented.
#   When added: create 3cpt_infusion.ctl, generate concordance CSV,
#   assert CL/V1/Q/V2/Q2/V2 within 10%.
# ---------------------------------------------------------------------------
test_that("FEATURE GAP [three_cpt_infusion]: 3-cpt infusion concordance", {
  skip("ferx-core gap: three_cpt_infusion pk macro not yet supported")
})

# ---------------------------------------------------------------------------
# Model-defined infusion rate / duration (R1= / D1= in $PK)
#   ferx-core gap: rate parameter in pk macro (R1/D1 passthrough).
#   When added: create r1_infusion.ctl, generate concordance CSV,
#   assert rate/duration params within 15%.
# ---------------------------------------------------------------------------
test_that("FEATURE GAP [R1_D1_infusion]: model-defined infusion rate/duration", {
  skip("ferx-core gap: R1/D1 model-defined infusion rate not yet supported")
})

# ---------------------------------------------------------------------------
# Multiple DVIDs / endpoints (e.g. PK+PD with DVID column)
#   ferx-core gap: multi-endpoint / DVID support.
#   When added: create pkpd_dvid.ctl, generate concordance CSV per endpoint,
#   assert PK and PD params within 15%.
# ---------------------------------------------------------------------------
test_that("FEATURE GAP [multiple_dvid]: multi-endpoint DVID concordance", {
  skip("ferx-core gap: multiple DVIDs / endpoints not yet supported")
})

# ---------------------------------------------------------------------------
# MIXTURE models
#   ferx-core gap: mixture model support.
#   When added: create mixture.ctl, generate concordance CSV,
#   assert mixture proportions and subpop params within 20%.
# ---------------------------------------------------------------------------
test_that("FEATURE GAP [mixture_model]: MIXTURE model concordance", {
  skip("ferx-core gap: MIXTURE models not yet supported")
})

# ---------------------------------------------------------------------------
# IOV with block kappa (off-diagonal IOV omega)
#   ferx-core gap: block kappa (currently only diagonal kappas emitted).
#   When added: extend iov.ctl to have correlated IOV, update concordance CSV,
#   assert block kappa elements within 20%.
# ---------------------------------------------------------------------------
test_that("FEATURE GAP [iov_block_kappa]: IOV with block kappa concordance", {
  skip("ferx-core gap: block kappa (off-diagonal IOV omega) not yet supported")
})

# ---------------------------------------------------------------------------
# Transit compartment absorption
#   ferx-core gap: transit compartment pk macro or ODE pattern.
#   When added: create transit.ctl (e.g. ADVAN6 transit chain), generate
#   concordance CSV, assert KTR/MTT/KA/CL/V within 15%.
# ---------------------------------------------------------------------------
test_that("FEATURE GAP [transit_compartments]: transit absorption concordance", {
  skip("ferx-core gap: transit compartment absorption not yet supported")
})
