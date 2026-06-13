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

# Strip machine-specific installed paths from header comment so snapshots are
# portable across machines and CI.  Reduces e.g.
#   "# Translated from nonmem: /Library/.../1cpt_oral.ctl"
# to
#   "# Translated from nonmem: 1cpt_oral.ctl"
norm_snap <- function(txt) {
  sub("(# Translated from [^:]+: ).*/([^/\n]+)", "\\1\\2", txt, perl = TRUE)
}

# -- NONMEM models ------------------------------------------------------------

test_that("1-cpt oral NONMEM: snapshot + no unsupported", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("1cpt_oral.ctl"))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "one_cpt_oral", fixed = TRUE)
})

test_that("2-cpt oral with covariates: snapshot + no unsupported", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("2cpt_oral_cov.ctl"))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "two_cpt_oral", fixed = TRUE)
})

test_that("2-cpt IV: infers two_cpt_iv", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("2cpt_iv.ctl"))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "two_cpt_iv", fixed = TRUE)
})

test_that("3-cpt IV: translates to three_cpt_iv pk macro", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("3cpt_iv.ctl"))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "three_cpt_iv", fixed = TRUE)
})

test_that("ODE warfarin: full $DES path, [odes] section present", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("ode_warfarin.ctl"))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "[odes]",       fixed = TRUE)
  expect_match(result$ferx_text, "d/dt(DEPOT)",  fixed = TRUE)
  expect_length(result$unsupported, 0L)
})

test_that("block omega: block_omega line in output", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("block_omega.ctl"))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "block_omega", fixed = TRUE)
})

test_that("IOV model: KAPPA_CL emitted as omega + flattening warning (nonmem2rx treats IOV as IIV)", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("iov.ctl"))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "KAPPA_CL", fixed = TRUE)
  # nonmem2rx flattens the ETA-coded IOV to IIV; the translator must warn so the
  # silent loss of occasion structure is visible to the user.
  expect_true(any(grepl("inter-occasion", result$warnings, fixed = TRUE)))
})

# -- nlmixr2 models -----------------------------------------------------------

test_that("1-cpt oral nlmixr2: snapshot + one_cpt_oral", {
  skip_if_not_installed("rxode2")
  fn     <- source(r2_path("1cpt_oral_nlmixr2.R"))$value
  result <- nlmixr2_to_ferx(fn)
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "one_cpt_oral", fixed = TRUE)
})

test_that("nlmixr2 source: a KAPPA-named IIV eta does NOT trigger the NONMEM-only IOV warning", {
  skip_if_not_installed("rxode2")
  fn     <- source(r2_path("iov_kappa_nlmixr2.R"))$value
  # suppressWarnings() silences rxode2's benign "non-mu referenced" parse note
  # for exp(eta.cl + kappa.cl); it does not touch result$warnings (the
  # translator's own channel), which is what the assertions below check.
  result <- suppressWarnings(nlmixr2_to_ferx(fn))
  # The eta is present in the IIV block (so the helper would match its name)...
  expect_match(result$ferx_text, "omega KAPPA_CL", fixed = TRUE)
  # ...but the flattening warning is nonmem2rx-specific and must stay silent here.
  expect_false(any(grepl("inter-occasion", result$warnings, fixed = TRUE)))
})

test_that("ODE nlmixr2: d/dt expressions produce [odes] section", {
  skip_if_not_installed("rxode2")
  fn     <- source(r2_path("ode_nlmixr2.R"))$value
  result <- nlmixr2_to_ferx(fn)
  expect_snapshot(cat(norm_snap(result$ferx_text)))
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

# -- amp.sim example models ---------------------------------------------------

test_that("amp.sim 1-cpt oral ODE: [odes] section + obs_cmt inferred", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("pk_1cmt_oral.mod"))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "[odes]",        fixed = TRUE)
  expect_match(result$ferx_text, "obs_cmt=",      fixed = TRUE)
  expect_match(result$ferx_text, "d/dt(",         fixed = TRUE)
  expect_match(result$ferx_text, "proportional",  fixed = TRUE)
  expect_length(result$unsupported, 0L)
})

test_that("pk_1cmt_oral.mod: S2=V scaling emits [scaling] obs_scale = V", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("pk_1cmt_oral.mod"))
  expect_match(result$ferx_text, "[scaling]",      fixed = TRUE)
  expect_match(result$ferx_text, "obs_scale = V",  fixed = TRUE)
  expect_true(any(grepl("S2 = V", result$warnings, fixed = TRUE)))
})

test_that("amp.sim PKPD indirect response: 4-state ODE + additive error", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("pkpd_ir.mod"))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "[odes]",   fixed = TRUE)
  expect_match(result$ferx_text, "obs_cmt=", fixed = TRUE)
  expect_match(result$ferx_text, "d/dt(",    fixed = TRUE)
  expect_match(result$ferx_text, "additive", fixed = TRUE)
  n_odes <- length(regmatches(result$ferx_text,
                              gregexpr("d/dt\\(", result$ferx_text))[[1]])
  expect_equal(n_odes, 4L)
})

test_that("pk_1cmt_oral_ampsim: fixed-effect V passthrough appears in pk macro", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("pk_1cmt_oral_ampsim.ctl"))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  # Fixed-effect V (no ETA) must appear as passthrough in [individual_parameters]
  # and be passed to the pk macro, otherwise ferx predicts zero concentration.
  expect_match(result$ferx_text, "V = V",              fixed = TRUE)
  expect_match(result$ferx_text, "v=V",                fixed = TRUE)
  expect_match(result$ferx_text, "one_cpt_oral",       fixed = TRUE)
  expect_length(result$unsupported, 0L)
})
