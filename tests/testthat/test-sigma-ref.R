# SIGMA(i,j) read in the model code (issue #46). nonmem2rx hands it over as a
# theta named sigma.i.j holding the initial variance, which used to be emitted as
# a free estimated theta: the residual variance declared twice, once as the
# `sigma` and once as an unlinked theta. ferx has no way to name a sigma in an
# expression, so the reference is dropped when nothing reads it and reported when
# something does.

.sigma_ctl <- function(error, v = "  V = THETA(2)", sigma = "$SIGMA 400",
                       y = "  Y = IPRED + EPS(1)") {
  ctl <- file.path(tmp_ctl_dir(), "g.ctl")
  writeLines(c(
    "$PROBLEM sigma reference", "$INPUT ID TIME DV AMT EVID MDV CMT", "$DATA d.csv IGNORE=@",
    "$SUBROUTINES ADVAN1 TRANS2",
    "$PK", "  CL = THETA(1)*EXP(ETA(1))", v, "  S1 = V",
    "$ERROR", "  IPRED = F", error, y,
    "$THETA", "  (0,3)", "  (0,30)", "$OMEGA 0.09", sigma, "$EST METHOD=1 INTER"), ctl)
  ctl
}

.sigma_run <- function(...) suppressWarnings(suppressMessages(
  nm_to_ferx(.sigma_ctl(...), validate = FALSE)))

.theta_lines <- function(txt) grep("^  theta ", strsplit(txt, "\n")[[1]], value = TRUE)

test_that("the weighted-residual idiom emits no theta for SIGMA(1,1)", {
  skip_if_not_installed("nonmem2rx")
  r <- .sigma_run(c("  W = SQRT(SIGMA(1,1))", "  IWRES = (DV - IPRED)/W"))
  expect_equal(length(.theta_lines(r$ferx_text)), 2L)
  expect_false(grepl("SIGMA_1_1", r$ferx_text, fixed = TRUE))
  # The dead W is not left behind for ferx to report as computed but never used.
  expect_false(grepl("sqrt(", r$ferx_text, fixed = TRUE))
  expect_false(any(grepl("^  W = ", strsplit(r$ferx_text, "\n")[[1]])))
  # The residual variance is declared exactly once.
  expect_equal(sum(grepl("^  sigma ", strsplit(r$ferx_text, "\n")[[1]])), 1L)
  expect_true(any(grepl("W reads a SIGMA(i,j) reference", r$warnings, fixed = TRUE)))
  expect_length(r$unsupported, 0L)
})

test_that("each of two referenced sigmas is dropped, and a live one is kept", {
  skip_if_not_installed("nonmem2rx")
  two <- c("  W1 = SQRT(SIGMA(1,1))", "  W2 = SQRT(SIGMA(2,2))")
  r <- .sigma_run(two, sigma = c("$SIGMA 0.04", "$SIGMA 400"),
                  y = "  Y = IPRED*(1 + EPS(1)) + EPS(2)")
  expect_false(grepl("SIGMA_", r$ferx_text, fixed = TRUE))
  expect_length(r$unsupported, 0L)

  # SIGMA(2,2) read by the volume: that one is live, SIGMA(1,1) is not.
  r <- .sigma_run(two, sigma = c("$SIGMA 0.04", "$SIGMA 400"),
                  v = "  V = THETA(2)*(1 + 0.5*SQRT(SIGMA(2,2)))",
                  y = "  Y = IPRED*(1 + EPS(1)) + EPS(2)")
  th <- .theta_lines(r$ferx_text)
  expect_false(any(grepl("SIGMA_1_1", th, fixed = TRUE)))
  expect_equal(sum(grepl("SIGMA_2_2", th, fixed = TRUE)), 1L)
})

test_that("a SIGMA(i,j) the model reads is a FIXED theta and an ERROR", {
  skip_if_not_installed("nonmem2rx")
  r <- .sigma_run("  W = SQRT(SIGMA(1,1))",
                  v = "  V = THETA(2)*(1 + 0.5*SQRT(SIGMA(1,1)))")
  th <- .theta_lines(r$ferx_text)
  sg <- grep("SIGMA_1_1", th, value = TRUE)
  # Never a free estimated theta: that is the sigma declared twice.
  expect_length(sg, 1L)
  expect_match(sg, "FIX", fixed = TRUE)
  expect_true(any(startsWith(r$warnings, "ERROR") &
                  grepl("SIGMA(1,1) is read by the translated model", r$warnings,
                        fixed = TRUE)))
  expect_true(any(grepl("SIGMA(1,1) reference", r$unsupported, fixed = TRUE)))
  # And it is announced where the reader will look.
  expect_match(r$ferx_text, "# WARNING: SIGMA(1,1)", fixed = TRUE)
})
