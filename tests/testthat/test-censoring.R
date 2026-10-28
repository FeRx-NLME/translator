# M3 censored-data likelihood in $ERROR (issue #49). A conditional that sets
# F_FLAG = 1 and Y = PHI(...) used to be dropped, so the censored records were
# fitted as ordinary measurements and every parameter came out biased. ferx
# scores them itself under `bloq_method = m3`, so the conditional is recognised
# and turned into that option.

.m3_ctl <- function(error, input = "ID TIME DV AMT EVID MDV CMT BLQ",
                    sigma = "$SIGMA 1 FIX") {
  ctl <- file.path(tmp_ctl_dir(), "m3.ctl")
  writeLines(c(
    "$PROBLEM m3", paste("$INPUT", input), "$DATA d.csv IGNORE=@",
    "$SUBROUTINES ADVAN2 TRANS2",
    "$PK", "  CL = THETA(1)*EXP(ETA(1))", "  V = THETA(2)", "  KA = THETA(3)", "  S2 = V",
    "$ERROR", "  IPRED = F", "  W = IPRED*THETA(4)", error,
    "$THETA", "  (0,3)", "  (0,30)", "  (0,2)", "  (0,0.3)", "$OMEGA 0.09", sigma,
    "$EST METHOD=1 INTER LAPLACIAN"), ctl)
  ctl
}

.m3_run <- function(...) suppressWarnings(suppressMessages(
  nm_to_ferx(.m3_ctl(...), validate = FALSE)))

.m3_then <- c("  LLOQ = 0.3", "  F_FLAG = 0", "  Y = IPRED + W*EPS(1)",
              "  IF (BLQ.EQ.1) THEN", "    F_FLAG = 1", "    Y = PHI((LLOQ-IPRED)/W)",
              "  ENDIF")

.lines <- function(txt) strsplit(txt, "\n")[[1]]

test_that("an IF-then M3 conditional becomes bloq_method = m3", {
  skip_if_not_installed("nonmem2rx")
  r <- .m3_run(.m3_then)
  expect_true(any(grepl("^  bloq_method = m3$", .lines(r$ferx_text))))
  # The error model is the observed-record one, intact.
  expect_true(grepl("DV ~ proportional(THETA4)", r$ferx_text, fixed = TRUE))
  # The action the user has to take is named, with the column the source tested.
  w <- grep("M3 censored-data likelihood", r$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "^WARN")
  expect_match(w, "named CENS", fixed = TRUE)
  expect_match(w, "`BLQ`", fixed = TRUE)
  expect_length(r$unsupported, 0L)
})

test_that("the M3 scaffolding leaves nothing behind in [individual_parameters]", {
  skip_if_not_installed("nonmem2rx")
  r <- .m3_run(.m3_then)
  ip <- .lines(r$ferx_text)
  ip <- ip[seq(which(ip == "[individual_parameters]"), which(ip == "[structural_model]"))]
  expect_false(any(grepl("LLOQ|F_FLAG|BLQ|phi", ip, ignore.case = TRUE)))
  expect_false(any(grepl("^(WARN|ERROR)", r$warnings) & grepl("never used", r$warnings)))
})

test_that("the IF/ELSE spelling is recognised and keeps the observed error model", {
  skip_if_not_installed("nonmem2rx")
  r <- .m3_run(c("  LLOQ = 0.3", "  F_FLAG = 0",
                 "  IF (BLQ.EQ.0) THEN", "    Y = IPRED*(1 + EPS(1))", "  ELSE",
                 "    F_FLAG = 1", "    Y = PHI((LLOQ-IPRED)/W)", "  ENDIF"),
               sigma = "$SIGMA 0.04")
  expect_true(any(grepl("^  bloq_method = m3$", .lines(r$ferx_text))))
  # Without recognition the whole [error_model] vanished (E_MISSING_BLOCK).
  expect_true(grepl("DV ~ proportional(EPS1)", r$ferx_text, fixed = TRUE))
})

test_that("an upper-limit 1 - PHI conditional is recognised and names the sign", {
  skip_if_not_installed("nonmem2rx")
  r <- .m3_run(c("  ULOQ = 9", "  F_FLAG = 0", "  Y = IPRED + W*EPS(1)",
                 "  IF (BLQ.EQ.1) THEN", "    F_FLAG = 1",
                 "    Y = 1 - PHI((IPRED-ULOQ)/W)", "  ENDIF"))
  expect_true(any(grepl("^  bloq_method = m3$", .lines(r$ferx_text))))
  expect_match(grep("M3 censored-data", r$warnings, value = TRUE), "-1 = above",
               fixed = TRUE)
})

test_that("a stream with no censoring is untouched (control)", {
  skip_if_not_installed("nonmem2rx")
  r <- .m3_run("  Y = IPRED + W*EPS(1)", input = "ID TIME DV AMT EVID MDV CMT")
  expect_false(grepl("bloq_method", r$ferx_text, fixed = TRUE))
  expect_false(any(grepl("M3 censored", r$warnings)))
})

test_that("an F_FLAG conditional that is not the M3 shape is an ERROR, never guessed at", {
  skip_if_not_installed("nonmem2rx")
  # Single-line form: sets the flag, defines no likelihood for it.
  r <- .m3_run(c("  F_FLAG = 0", "  Y = IPRED + W*EPS(1)", "  IF (BLQ.EQ.1) F_FLAG = 1"))
  expect_false(grepl("bloq_method", r$ferx_text, fixed = TRUE))
  expect_true(any(grepl("^ERROR \\| a conditional assigns F_FLAG", r$warnings)))
  expect_true(any(grepl("F_FLAG conditional", r$unsupported)))
  # An extra statement in the censored branch is also outside the recognised shape.
  r2 <- .m3_run(c("  F_FLAG = 0", "  Y = IPRED + W*EPS(1)", "  IF (BLQ.EQ.1) THEN",
                  "    F_FLAG = 1", "    Z = 2", "    Y = PHI((0.3-IPRED)/W)", "  ENDIF"))
  expect_false(grepl("bloq_method", r2$ferx_text, fixed = TRUE))
  expect_true(any(grepl("F_FLAG conditional", r2$unsupported)))
})

test_that("the M3 detector only runs on NONMEM sources", {
  e <- quote(if (BLQ == 1) { f_flag <- 1; y <- phi((LLOQ - ipred)/w) })
  off <- .extract_m3(list(e), enabled = FALSE)
  expect_false(off$found)
  expect_identical(off$lst, list(e))
  expect_true(.extract_m3(list(e), enabled = TRUE)$found)
})
