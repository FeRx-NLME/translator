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
#   Note: the amp.sim benchmark uses the linCmt path only. The ODE variant
#   (pk_1cmt_oral.mod / ADVAN6, with S2=V scaling) is exercised separately by
#   the ODE concordance test below, not against the amp.sim reference.
#
# Acceptance criteria (per plans/v0.1-implementation.md Section 10):
#   structural thetas: within 15% of truth
#   random-effect variances: within 20% of truth

skip_if_not_installed("ferx")
# These Tier-4 tests exercise the ferx engine. In CI they run ONLY in the
# dedicated "engine" job, which installs a pinned ferx and sets
# FERXTRANSLATE_ENGINE_TESTS=true; the fast PR job (no ferx) skips them. Locally
# they run whenever ferx is installed -- same as the old skip_on_ci() behaviour.
# Force a local run of just this tier with:
#   FERXTRANSLATE_ENGINE_TESTS=true Rscript -e 'devtools::test(filter="concordance")'
if (tolower(Sys.getenv("CI")) %in% c("true", "1") &&
    !identical(Sys.getenv("FERXTRANSLATE_ENGINE_TESTS"), "true"))
  skip("engine (Tier-4) tests run only in the CI 'engine' job (pinned ferx)")

library(ferxtranslate)
library(ferx)

# Helper: translate a bundled NONMEM model and write .ferx to a temp file
.translate_to_tmp <- function(model_name) {
  ctl <- system.file(file.path("testmodels/nonmem", model_name),
                     package = "ferxtranslate")
  result <- suppressWarnings(nm_to_ferx(ctl))
  ferx_file <- tempfile(fileext = ".ferx")
  writeLines(result$ferx_text, ferx_file)
  ferx_file
}

# Helper: bundled concordance dataset path
# Named diagonal of the fitted omega. unlist() on a matrix is a no-op, so
# positional indexing into it silently means "first and last cell", not the two
# etas -- and drops ETA_CL entirely the moment a third eta is added.
.omega_diag <- function(fit) {
  om <- fit$omega
  if (is.matrix(om)) {
    d <- diag(om)
    nm <- rownames(om) %||% colnames(om)
    if (!is.null(nm)) names(d) <- nm
    return(d)
  }
  unlist(om)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

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
# 1-cpt oral (linCmt -> one_cpt_oral pk macro)
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
# 2-cpt IV (linCmt -> two_cpt_iv pk macro)
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

test_that("2-cpt IV: omega estimates reproduce the reference fit", {
  ferx_file <- .translate_to_tmp("2cpt_iv.ctl")
  data_file  <- .conc_data("2cpt_iv_concordance.csv")
  fit <- ferx_fit(ferx_file, data_file,
                  method     = "focei",
                  covariance = FALSE,
                  verbose    = FALSE)

  # Unlike the structural thetas (well identified, so they recover the nominal
  # simulation truth within 10%), omega variances from 50 subjects carry a
  # sampling SE of ~omega*sqrt(2/50) ~ 20%, so the ML estimate genuinely
  # departs from the nominal 0.10/0.10/0.08/0.05 used to simulate. Asserting
  # against the nominal truth would be statistically unsound at this N.
  #
  # Instead assert that translate + fit reproduces the *reference fit* -- the
  # ML omegas a known-good run yields on this fixed dataset (ferx 0.1.x FOCEI,
  # mu-referenced). This is a deterministic regression check on the eta-to-
  # parameter wiring: a swapped/missing eta or wrong IIV structure shifts these
  # by >2x and trips the 10% tolerance, while normal cross-platform numerical
  # noise (<1%) does not. Nominal simulation values, for provenance:
  #   ETA_CL 0.10, ETA_V1 0.10, ETA_Q 0.08, ETA_V2 0.05.
  omega_diag <- c(ETA_CL = fit$omega["ETA_CL", "ETA_CL"],
                  ETA_V1 = fit$omega["ETA_V1", "ETA_V1"],
                  ETA_Q  = fit$omega["ETA_Q",  "ETA_Q"],
                  ETA_V2 = fit$omega["ETA_V2", "ETA_V2"])
  ref_omega  <- c(ETA_CL = 0.08535, ETA_V1 = 0.10094,
                  ETA_Q  = 0.05482, ETA_V2 = 0.02688)
  .report_deviations(omega_diag, ref_omega, "2-cpt IV omegas (vs reference fit)")
  expect_lt(abs(omega_diag["ETA_CL"] / ref_omega["ETA_CL"] - 1), 0.10, label = "omega_CL")
  expect_lt(abs(omega_diag["ETA_V1"] / ref_omega["ETA_V1"] - 1), 0.10, label = "omega_V1")
  expect_lt(abs(omega_diag["ETA_Q"]  / ref_omega["ETA_Q"]  - 1), 0.10, label = "omega_Q")
  expect_lt(abs(omega_diag["ETA_V2"] / ref_omega["ETA_V2"] - 1), 0.10, label = "omega_V2")
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

  # Theta names carry the TV prefix: the source names them KA/CL/V, which would
  # shadow the identically named individual parameters, so the translator
  # renames the thetas (see .deshadow_theta_names()).
  ref_nm <- c(TVKA = ref$THETA1, TVCL = ref$THETA2, TVV = ref$THETA3)
  .report_deviations(fit$theta, ref_nm, "amp.sim 1-cpt oral thetas vs NONMEM reference")
  expect_lt(abs(fit$theta["TVKA"] / ref_nm["TVKA"] - 1), 0.10, label = "KA vs amp.sim ref")
  expect_lt(abs(fit$theta["TVCL"] / ref_nm["TVCL"] - 1), 0.10, label = "CL vs amp.sim ref")
  expect_lt(abs(fit$theta["TVV"]  / ref_nm["TVV"]  - 1), 0.10, label = "V vs amp.sim ref")
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
  # Start ETA_CL's variance away from the truth it must recover, so that
  # "returned its initial value" and "recovered the truth" are distinguishable.
  start_om  <- 0.005
  txt <- sub("omega ETA_CL ~ [^\n]+", paste0("omega ETA_CL ~ ", start_om),
             paste(readLines(ferx_file), collapse = "\n"))
  expect_match(txt, paste0("omega ETA_CL ~ ", start_om), fixed = TRUE)
  writeLines(txt, ferx_file)
  fit <- ferx_fit(ferx_file, data_file,
                  method     = "focei",
                  covariance = FALSE,
                  verbose    = FALSE)

  # KA and CL are renamed away from the individual parameters that would shadow
  # them; V has no ETA, so it stays a plain theta. Before that fix, ETA_CL had
  # no effect on this model at all -- K20 = CL/V read the theta.
  ref <- c(TVKA = 0.1, TVCL = 2.0, V = 1.0)
  .report_deviations(fit$theta, ref, "ODE 1-cpt oral thetas")
  expect_lt(abs(fit$theta["TVKA"] / ref["TVKA"] - 1), 0.15, label = "KA")
  expect_lt(abs(fit$theta["TVCL"] / ref["TVCL"] - 1), 0.15, label = "CL")
  expect_lt(abs(fit$theta["V"]    / ref["V"]    - 1), 0.15, label = "V")

  # Regression guard for the theta-shadowing defect. CL reaches the ODEs only
  # through the derived K20 = CL/V, which lives in [individual_parameters] --
  # the one block where a theta named CL would shadow the individual CL. When it
  # did, ETA_CL had no gradient and its omega came back untouched at its initial
  # value.
  #
  # Asserting recovery of the truth is NOT enough on its own: the model's $OMEGA
  # initials ARE the simulation truth (0.01, 0.02), so the broken case returns
  # exactly the asserted value and passes with a deviation of ~1e-15. The fit is
  # therefore started away from the truth, so only a live gradient can reach it.
  om <- .omega_diag(fit)
  expect_named(om, c("ETA_KA", "ETA_CL"))
  expect_lt(abs(om[["ETA_KA"]] / 0.01 - 1), 0.35, label = "omega ETA_KA")
  expect_lt(abs(om[["ETA_CL"]] / 0.02 - 1), 0.35, label = "omega ETA_CL")
  # A dead parameter cannot move off its start; a live one must.
  expect_gt(abs(om[["ETA_CL"]] / start_om - 1), 0.20, label = "omega ETA_CL moved")
})

# ===========================================================================
# ENGINE VALIDATION SWEEP (Tier 4)
# Validates every bundled model's emitted .ferx with the engine and fails on
# error-severity diagnostics. The engine-free half of this sweep -- translate
# every model, fail on a crash, report $unsupported -- lives in
# test-integration.R so it runs on every PR without a Rust build.
# ===========================================================================
test_that("engine accepts the emitted .ferx for every bundled model", {
  skip_if_not_installed("nonmem2rx")

  models <- .bundled_nm_models()
  expect_gt(length(models), 0L)

  rows <- lapply(models, function(path) {
    result <- tryCatch(
      suppressWarnings(suppressMessages(nm_to_ferx(path, strict = FALSE))),
      error = function(e) NULL
    )
    if (is.null(result)) return(NULL)   # crashes are the integration test's job
    diags <- result$validation$diagnostics
    if (is.null(diags) || nrow(diags) == 0) return(NULL)
    data.frame(model    = basename(path),
               severity = as.character(diags$severity),
               finding  = paste0(diags$code, ": ", .one_line(diags$message)),
               stringsAsFactors = FALSE)
  })
  found <- do.call(rbind, rows)

  if (is.null(found) || nrow(found) == 0) {
    message("engine validation sweep: no diagnostics across ", length(models),
            " models")
  } else {
    message("\nengine validation sweep (", nrow(found), " diagnostic(s) across ",
            length(models), " models):")
    message(paste(capture.output(print(found, row.names = FALSE)), collapse = "\n"))
  }

  # Engine warnings are reported; engine errors fail. paste0() recycles
  # zero-length arguments against the length-1 separator, so guard on nrow().
  fatal <- character()
  if (!is.null(found) && nrow(found) > 0) {
    err <- found[which(found$severity == "error"), , drop = FALSE]
    if (nrow(err) > 0) fatal <- paste0(err$model, " -- ", err$finding)
  }
  expect_equal(fatal, character())
})
