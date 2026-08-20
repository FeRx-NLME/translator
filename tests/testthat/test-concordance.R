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
  # BOTH etas, not just one: pk_1cmt_oral.mod declares $OMEGA .01/.02, which are
  # also the simulation truth, so a dead parameter returning its untouched
  # initial is indistinguishable from a live one recovering the truth.
  start_ka <- 0.004
  start_cl <- 0.005
  txt <- paste(readLines(ferx_file), collapse = "\n")
  txt <- sub("omega ETA_KA ~ [^\n]+", paste0("omega ETA_KA ~ ", start_ka), txt)
  txt <- sub("omega ETA_CL ~ [^\n]+", paste0("omega ETA_CL ~ ", start_cl), txt)
  expect_match(txt, paste0("omega ETA_KA ~ ", start_ka), fixed = TRUE)
  expect_match(txt, paste0("omega ETA_CL ~ ", start_cl), fixed = TRUE)
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
  expect_gt(abs(om[["ETA_KA"]] / start_ka - 1), 0.20, label = "omega ETA_KA moved")
  expect_gt(abs(om[["ETA_CL"]] / start_cl - 1), 0.20, label = "omega ETA_CL moved")
})

test_that("issue #25: the dose reaches the right compartment when $MODEL and d/dt order differ", {
  # The only tier that can catch this. cmt_order_gap.ctl declares DUMMY as
  # $MODEL compartment 2 with no DADT, so nonmem2rx materialises `d/dt(DUMMY) =
  # 0` and places it FIRST: d/dt order is [DUMMY, DEPOT, CENTRAL], $MODEL order
  # is [DEPOT, DUMMY, CENTRAL]. ferx numbers compartments by position in
  # `states=[...]`, so emitting d/dt order put this dataset's CMT=1 dose into
  # DUMMY, whose derivative is zero.
  #
  # This fixture is unusually discriminating and that is deliberate: the broken
  # case does not merely estimate badly, it makes IPRED identically 0 for every
  # subject at every time, because no drug ever enters the system. No choice of
  # starting values can rescue that, which is why the structural thetas are
  # asserted against nominal truth here without the perturbation the omega
  # guards need. Measured on ferx 0.3.0: 0.000000 at every observation time
  # before the fix, and the reference profile after it.
  ferx_file <- .translate_to_tmp("cmt_order_gap.ctl")
  data_file <- .conc_data("cmt_order_gap_concordance.csv")

  # The emitted numbering is what the dataset is read against, so pin it here
  # too -- a fit that passed with the states in some other order would mean the
  # dataset had been regenerated to match the bug.
  txt <- paste(readLines(ferx_file), collapse = "\n")
  expect_match(txt, "states=[DEPOT, DUMMY, CENTRAL]", fixed = TRUE)

  # Every theta started at HALF its true value, and this is not optional.
  #
  # An earlier version of this test asserted nominal truth from the model's own
  # initials and was a tautology, which a sabotage run caught: with the dose in
  # DUMMY there is no gradient at all, so the optimiser returns every theta
  # exactly where it started -- and cmt_order_gap.ctl's $THETA initials ARE the
  # simulation truth. Measured, broken: theta = (3, 50, 1.2), deviation ~0, test
  # green. Starting away from the truth makes "never moved" and "recovered it"
  # different observations. Same trap as the omega guard below, in the tier
  # above it. See CLAUDE.md's discriminating-fixture rule.
  start <- c(TVCL = 1.5, TVV = 25.0, TVKA = 0.6)
  for (nm in names(start)) {
    pat <- sprintf("theta %s\\([^)]+\\)", nm)
    if (!grepl(pat, txt))
      stop("no theta '", nm, "' in the emitted .ferx -- has the naming changed?")
    txt <- sub(pat, sprintf("theta %s(%s, 0.0, 1e15)", nm, start[[nm]]), txt)
  }
  writeLines(txt, ferx_file)

  fit <- ferx_fit(ferx_file, data_file,
                  method     = "focei",
                  covariance = FALSE,
                  verbose    = FALSE)

  ref <- c(TVCL = 3.0, TVV = 50.0, TVKA = 1.2)
  .report_deviations(fit$theta, ref, "issue #25 compartment numbering")
  expect_lt(abs(fit$theta["TVCL"] / ref["TVCL"] - 1), 0.15, label = "CL")
  expect_lt(abs(fit$theta["TVV"]  / ref["TVV"]  - 1), 0.15, label = "V")
  expect_lt(abs(fit$theta["TVKA"] / ref["TVKA"] - 1), 0.15, label = "KA")
  # A dead parameter cannot move off its start; a live one must. This is what
  # separates a fit from a frozen optimiser, and it is the assertion that fails
  # first when the dose stops reaching the system.
  for (nm in names(start))
    expect_gt(abs(fit$theta[[nm]] / start[[nm]] - 1), 0.20,
              label = paste0("theta ", nm, " moved off its start"))
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

test_that("the engine sweep reports non-error diagnostics without failing", {
  # Kills the mutation removing the nrow(err) guard. paste0() recycles
  # zero-length arguments against the length-1 separator, so a frame with only
  # warning rows yielded " -- " instead of character(0) and failed the sweep on
  # exactly the case it is meant to print. The live sweep cannot exercise this
  # while the corpus is clean, so drive the expression directly.
  fatal_for <- function(found) {
    fatal <- character()
    if (!is.null(found) && nrow(found) > 0) {
      err <- found[which(found$severity == "error"), , drop = FALSE]
      if (nrow(err) > 0) fatal <- paste0(err$model, " -- ", err$finding)
    }
    fatal
  }
  warn_only <- data.frame(model = "m.ctl", severity = "warning",
                          finding = "W_X: something", stringsAsFactors = FALSE)
  expect_equal(fatal_for(warn_only), character())
  expect_equal(fatal_for(NULL), character())

  mixed <- rbind(warn_only,
                 data.frame(model = "n.ctl", severity = "error",
                            finding = "E_PARSE: bad", stringsAsFactors = FALSE))
  expect_equal(fatal_for(mixed), "n.ctl -- E_PARSE: bad")
})

# -- accidental dose attributes (issue #17) -----------------------------------
#
# The one tier that can see this defect at all. It is invisible to every other:
# the wrong model VALIDATES clean, `$unsupported` is empty, and no engine
# diagnostic mentions bioavailability -- so only comparing predictions against a
# spelling that cannot be a dose attribute separates the two readings.
test_that("a parameter named like a dose attribute is not applied to the dose", {
  skip_if_not_installed("rxode2")

  # An ordinary elimination rate constant that happens to be called F1. In
  # nlmixr2 bioavailability is written `f(central) <- ...`, so a bare `F1 <- `
  # assignment is an ordinary parameter and nothing in the source asks for a
  # dose attribute.
  m <- function() {
    ini({ tvv <- 50; tvf1 <- 0.1; eta.v ~ 0.04; prop.err <- 0.01 })
    model({
      v  <- tvv * exp(eta.v)
      F1 <- tvf1
      d/dt(central) <- -F1 * central
      cp <- central / v
      cp ~ prop(prop.err)
    })
  }
  emitted <- suppressWarnings(
    nlmixr2_to_ferx(rxode2::rxode2(m), validate = FALSE))$ferx_text

  # The guard must have fired, or the rest of this test compares a model with
  # itself and passes for the wrong reason.
  expect_match(emitted, "F1_PAR", fixed = TRUE)

  tmpl <- rbind(
    data.frame(ID = 1L, TIME = 0,   DV = 0, EVID = 1L, AMT = 100, CMT = 1L, MDV = 1L),
    data.frame(ID = 1L, TIME = c(0.5, 2, 8, 24), DV = 0, EVID = 0L, AMT = 0,
               CMT = 1L, MDV = 0L))
  data_file <- tempfile(fileext = ".csv")
  write.csv(tmpl, data_file, row.names = FALSE, quote = FALSE)

  write_tmp <- function(txt) {
    f <- tempfile(fileext = ".ferx")
    writeLines(txt, f)
    f
  }
  sim_of <- function(txt)
    ferx_simulate(write_tmp(txt), data_file, n_sim = 1L, seed = 20260819)$DV_SIM

  # `KE` is the control: a name ferx cannot read as a dose attribute, and
  # otherwise the same file. `F1` is what this package emitted before the guard.
  fixed   <- sim_of(emitted)
  control <- sim_of(gsub("F1_PAR", "KE", emitted, fixed = TRUE))

  # The rename changed the spelling and nothing else.
  expect_equal(fixed, control, tolerance = 1e-10)

  # The discriminating half, and it CHANGED SHAPE when the engine pin moved to
  # ferx-r@9c97c13 -- worth reading before touching it.
  #
  # This used to assert that the un-renamed spelling still simulates and comes
  # back scaled by exactly the parameter's value: `broken / fixed == 0.1`,
  # because the engine applied `F1` a second time as bioavailability on the dose,
  # silently. That was the whole hazard -- a wrong answer with no diagnostic.
  #
  # ferx-core#993/#1003 (`E_DOSE_ATTR_DOUBLE_USE`) turned that into a PARSE
  # ERROR: a dose attribute both applied by the engine and read by the model is
  # now rejected outright on ODE models. So the broken spelling no longer
  # simulates at all, and the old assertion failed with a zero-length `broken`
  # rather than a ratio -- which is how the pin bump surfaced it.
  #
  # The property under test is unchanged: WITHOUT the rename this package would
  # emit a model that does not mean what the source said. Only the consequence
  # moved, from silent-wrong to loud -- which is strictly better and is the point
  # of the upstream change. Asserting the rejection keeps the test discriminating;
  # dropping the half entirely would leave `fixed == control` comparing a model
  # with itself.
  #
  # NOTE this now REQUIRES the pinned engine. On a ferx older than ferx-core
  # 25b5f473 the broken spelling still validates, and this expectation fails --
  # correctly, since the concordance tier is tied to the pinned build.
  broken_v <- ferx_model_validate(write_tmp(gsub("F1_PAR", "F1", emitted,
                                                 fixed = TRUE)))
  expect_false(isTRUE(broken_v$ok))
  # Not merely "invalid": invalid FOR THIS REASON. The engine's message ends with
  # the literal clause below, which ferx-core pins with a test of its own so it
  # cannot silently degrade to a generic parse error.
  expect_match(paste(broken_v$diagnostics$message, collapse = " | "),
               "reserved dose-attribute name", fixed = TRUE)
})

# -- error-model sigma ORDER (issue #6 defect 10, restated) -------------------
#
# The issue calls this "the emitted combined() arguments are transposed". Measured
# against ferx 0.3.0, that lever does not exist: `combined(EPS1, EPS2)` and
# `combined(EPS2, EPS1)` fit to the same OFV to every digit. ferx consumes a
# single-endpoint error model's sigmas POSITIONALLY from the declaration order and
# discards the names (`build_error_spec`, ferx-core model_parser.rs:11328).
#
# The transposition the issue describes is real. What causes it is the order of
# the `sigma` lines in [parameters], so that is what this guards. A test written
# against the combined() arguments would pass whether or not the bug is fixed --
# it was, until this was measured.
test_that("sigma declaration order follows the error roles, and is load-bearing", {
  ctl <- c(
    "$PROBLEM combined error, additive term written first",
    "$INPUT ID TIME AMT EVID MDV CMT DV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINES ADVAN13 TOL=9",
    "$MODEL COMP=(CENT, DEFDOSE, DEFOBS)",
    "$PK",
    "  KE = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  S1 = V",
    "$DES",
    "  DADT(1) = -KE*A(1)",
    "$ERROR",
    "  IPRED = A(1)/V",
    "  Y = IPRED + EPS(1) + IPRED*EPS(2)",
    "$THETA (0.001, 0.1, 1.0)",
    "$THETA (1.0, 10.0, 100.0)",
    "$OMEGA 0.09",
    "$SIGMA 4.0",      # EPS(1): additive,     SD 2.0
    "$SIGMA 0.0025",   # EPS(2): proportional, SD 0.05
    "$ESTIMATION METHOD=1 INTER MAXEVAL=9999")
  d <- tmp_ctl_dir()
  writeLines(ctl, file.path(d, "m.ctl"))
  emitted <- suppressWarnings(nm_to_ferx(file.path(d, "m.ctl"),
                                         validate = FALSE))$ferx_text

  # EPS(2) is the proportional term, so it must be DECLARED first -- the source
  # declares it second. The 40x gap between the two SDs means a transposition
  # cannot hide in rounding.
  decl <- grep("^\\s*sigma ", strsplit(emitted, "\n")[[1]], value = TRUE)
  expect_equal(trimws(decl), c("sigma EPS2 ~ 0.05 (sd)", "sigma EPS1 ~ 2.0 (sd)"))
  # The arguments are in role order too. Inert today, but it is what the file
  # means, and it is what keeps the file right if ferx-core ever binds by name.
  expect_match(emitted, "DV ~ combined(EPS2, EPS1)", fixed = TRUE)

  rows <- do.call(rbind, lapply(1:6, function(i) rbind(
    data.frame(ID = i, TIME = 0, AMT = 100, EVID = 1L, MDV = 1L, CMT = 1L, DV = 0),
    data.frame(ID = i, TIME = c(0.5, 2, 8, 24), AMT = 0, EVID = 0L, MDV = 0L,
               CMT = 1L, DV = c(9.5, 8.2, 4.5, 0.9)))))
  data_file <- file.path(d, "d.csv")
  write.csv(rows, data_file, row.names = FALSE, quote = FALSE)

  # maxiter = 0 holds the thetas at their initials, so the OFV is a deterministic
  # function of the model text alone.
  ofv_of <- function(txt) {
    txt <- sub("maxiter = 500", "maxiter = 0", txt, fixed = TRUE)
    txt <- sub("covariance = true", "covariance = false", txt, fixed = TRUE)
    f <- tempfile(tmpdir = d, fileext = ".ferx")
    writeLines(txt, f)
    invisible(utils::capture.output(fit <- suppressMessages(
      ferx_fit(f, data = data_file))))
    fit$ofv
  }
  correct <- ofv_of(emitted)
  swapped <- ofv_of(sub("sigma EPS2 ~ 0.05 (sd)\n  sigma EPS1 ~ 2.0 (sd)",
                        "sigma EPS1 ~ 2.0 (sd)\n  sigma EPS2 ~ 0.05 (sd)",
                        emitted, fixed = TRUE))
  # The discriminating half. If this ever stops holding, the declaration order has
  # stopped mattering and the assertion above has become decorative.
  expect_false(isTRUE(all.equal(correct, swapped, tolerance = 1e-6)))

  # And the transposed ARGUMENT order really is inert -- the fact that makes the
  # declaration order the only lever. Stated as a test so it cannot quietly
  # change under us.
  arg_swapped <- ofv_of(sub("combined(EPS2, EPS1)", "combined(EPS1, EPS2)",
                            emitted, fixed = TRUE))
  expect_equal(arg_swapped, correct, tolerance = 1e-9)
})

# ===========================================================================
# ENDPOINT DISPATCH (Tier 4, issue #6 defect 5)
# The readout and the error model both branch per observation. Neither is
# visible to a snapshot -- a dropped dispatch still emits a valid file that
# fits -- so this is the tier that can tell a working one from a lost one.
# ===========================================================================

# A forward evaluation, not a fit: `maxiter = 0` makes PRED a deterministic
# function of the theta initials, so two models differing only in their readout
# are directly comparable. PRED is used rather than IPRED because it is eta-free;
# a different readout changes the EBEs, so IPRED would differ for reasons that
# have nothing to do with the branch selected.
.eval_only <- function(txt) sub("covariance = true", "covariance = false",
                                sub("maxiter = \\d+", "maxiter = 0", txt))

.write_tmp <- function(txt) {
  f <- tempfile(fileext = ".ferx")
  writeLines(txt, f)
  f
}

.pred_of <- function(txt, data) {
  fit <- suppressWarnings(ferx_fit(.write_tmp(txt), data = data))
  list(pred = fit$sdtab$PRED, cmt = fit$sdtab$CMT, ofv = fit$ofv)
}

test_that("a FLAG dispatch evaluates both readouts and selects on FLAG", {
  txt  <- .eval_only(suppressWarnings(
    nm_to_ferx(system.file("testmodels/nonmem/qss_tmdd.mod",
                           package = "ferxtranslate")))$ferx_text)
  data <- .conc_data("qss_tmdd_dispatch.csv")
  expect_true(nzchar(data))
  expect_match(txt, "y = if (FLAG == 2) c_RTOT else CENT/VC", fixed = TRUE)

  base <- .pred_of(txt, data)
  # Each branch hardcoded for every row. CMT and FLAG are correlated in this
  # dataset by construction, so CMT indexes the rows the FLAG condition selects.
  conc <- .pred_of(sub("y = if \\(FLAG == 2\\) c_RTOT else CENT/VC",
                       "y = CENT/VC", txt), data)
  rtot <- .pred_of(sub("y = if \\(FLAG == 2\\) c_RTOT else CENT/VC",
                       "y = c_RTOT", txt), data)
  a <- base$cmt == 1L
  b <- base$cmt == 3L
  expect_true(any(a) && any(b))

  expect_equal(max(abs(base$pred[a] - conc$pred[a])), 0)
  expect_equal(max(abs(base$pred[b] - rtot$pred[b])), 0)
  # The wrong-branch controls. Without them the two zeros above are equally
  # consistent with the engine evaluating one arm for every row.
  expect_gt(max(abs(base$pred[b] - conc$pred[b])), 1)
  expect_gt(max(abs(base$pred[a] - rtot$pred[a])), 1)

  # FLAG, not CMT, is what selects. Inverting FLAG while leaving CMT alone must
  # move every prediction; if the engine were keying on CMT -- or ignoring the
  # condition -- this would change nothing, and the assertions above would pass
  # just the same on a dataset where the two agree.
  d <- utils::read.csv(data)
  d$FLAG <- ifelse(d$FLAG == 1L, 2L, 1L)
  flipped <- tempfile(fileext = ".csv")
  utils::write.csv(d, flipped, row.names = FALSE, quote = FALSE)
  expect_gt(max(abs(base$pred - .pred_of(txt, flipped)$pred)), 1)
})

test_that("a CMT dispatch evaluates both readouts and selects on CMT", {
  txt  <- .eval_only(suppressWarnings(
    nm_to_ferx(system.file("testmodels/nonmem/pkpd_cmt.mod",
                           package = "ferxtranslate")))$ferx_text)
  data <- .conc_data("pkpd_cmt_dispatch.csv")
  expect_true(nzchar(data))
  expect_match(txt, "y[CMT=1] = CENT/VC", fixed = TRUE)
  expect_match(txt, "y[CMT=2] = PD", fixed = TRUE)

  base <- .pred_of(txt, data)
  one  <- "  y\\[CMT=1\\] = CENT/VC\n  y\\[CMT=2\\] = PD"
  conc <- .pred_of(sub(one, "  y = CENT/VC", txt), data)
  resp <- .pred_of(sub(one, "  y = PD", txt), data)
  a <- base$cmt == 1L
  b <- base$cmt == 2L
  expect_true(any(a) && any(b))

  expect_equal(max(abs(base$pred[a] - conc$pred[a])), 0)
  expect_equal(max(abs(base$pred[b] - resp$pred[b])), 0)
  expect_gt(max(abs(base$pred[b] - conc$pred[b])), 1)
  expect_gt(max(abs(base$pred[a] - resp$pred[a])), 1)
})

test_that("the per-endpoint sigma names bind, so the branch mapping is load-bearing", {
  # Measured on 0.3.0: per-CMT and covariate-selected error models resolve their
  # sigma BY NAME, unlike the single-endpoint form, which consumes them
  # positionally and discards the names (ferx-core#1001). So swapping the two
  # names between branches must move the OFV -- and if it does not, the mapping
  # this phase derives is decorative and a transposed one would ship green.
  for (m in c("qss_tmdd.mod", "pkpd_cmt.mod")) {
    txt  <- .eval_only(suppressWarnings(
      nm_to_ferx(system.file(file.path("testmodels/nonmem", m),
                             package = "ferxtranslate")))$ferx_text)
    data <- .conc_data(sub("\\.mod$", "_dispatch.csv", m))
    # Line-scoped: EPS1 <-> EPS2 inside the [error_model] block and nowhere else,
    # so the [parameters] declarations keep their values and only the mapping
    # moves.
    ln  <- strsplit(txt, "\n")[[1]]
    blk <- which(ln == "[error_model]")
    end <- which(!nzchar(ln) & seq_along(ln) > blk)[1]
    rng <- (blk + 1L):(end - 1L)
    # Only the sigma NAMES. chartr() over the whole line also flipped the
    # `FLAG == 2` condition, which turns `if (FLAG == 2) EPS2 else EPS1` into
    # `if (FLAG == 1) EPS1 else EPS2` -- the same model spelled backwards, and
    # an unmoved OFV that reads as "the names do not bind".
    ln[rng] <- gsub("EPS_TMP_", "EPS",
                    gsub("EPS2", "EPS_TMP_1",
                         gsub("EPS1", "EPS_TMP_2", ln[rng])))
    swapped <- paste(ln, collapse = "\n")
    expect_false(identical(swapped, txt), label = paste(m, "swap changed the text"))
    expect_gt(abs(.pred_of(txt, data)$ofv - .pred_of(swapped, data)$ofv), 1e-6)
  }
})
