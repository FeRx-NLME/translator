# $OMEGA BLOCK(1) SAME: a shared occasion variance becomes a kappa.
#
# nonmem2rx keeps only `fix = TRUE` on the repeated rows, so these tests drive
# the raw-record reader, the override rewriter, and the whole pipeline. Every
# control stream below is written for the test; none is a bundled analysis.

# The bundled SAME fixture with whole-line substitutions applied (a line matches
# when it equals the pattern after trimming), written to a temp control stream.
# Stops if a substitution matches nothing: a silent no-op here would turn a test
# of a variant into a second test of the original.
.same_ctl <- function(..., name = "m.ctl") {
  txt <- readLines(nm_path("iov_same.ctl"))
  subs <- list(...)
  for (pat in names(subs)) {
    hit <- trimws(txt) == pat
    if (!any(hit)) stop("no line matches '", pat, "'", call. = FALSE)
    txt <- unlist(lapply(seq_along(txt), function(i)
      if (hit[i]) subs[[pat]] else txt[i]), use.names = FALSE)
  }
  f <- file.path(tmp_ctl_dir(), name)
  writeLines(txt, f)
  f
}

.same_run <- function(ctl) {
  skip_if_not_installed("nonmem2rx")
  suppressWarnings(suppressMessages(nm_to_ferx(ctl, validate = FALSE)))
}

.omega_lines <- function(r) grep("^\\s*omega ", strsplit(r$ferx_text, "\n")[[1]], value = TRUE)

# -- the $OMEGA reader --------------------------------------------------------

test_that(".nm_omega_record_info counts the etas of every record spelling", {
  n <- function(txt) .nm_omega_record_info(txt)$n
  expect_equal(n(" 0.1 0.2 0.3"), 3L)
  expect_equal(n(" (0.1)x3 (0.2 FIX) 0.3"), 5L)
  expect_equal(n(" 0 FIX"), 1L)
  expect_equal(n(" BLOCK(2) 0.1 0.01 0.2"), 2L)
  expect_equal(n(" BLOCK(1) SAME"), 1L)
  expect_equal(n(" DIAG(3) 1 2 3"), 3L)
  expect_true(.nm_omega_record_info(" BLOCK(1) SAME")$same)
  expect_false(.nm_omega_record_info(" BLOCK(1) 0.04")$same)
})

test_that(".extract_nm_same groups a block with the SAME records that repeat it", {
  g <- .extract_nm_same(nm_path("iov_same.ctl"))
  expect_length(g$groups, 1L)
  expect_equal(g$groups[[1]]$etas, 2:4)
  expect_equal(g$groups[[1]]$block_n, 1L)
  expect_equal(g$n_eta, 4L)
  expect_length(g$problems, 0L)
})

test_that(".extract_nm_same is NULL for a stream with no SAME, including in a comment", {
  ctl <- .same_ctl("$OMEGA BLOCK(1) SAME" = "$OMEGA BLOCK(1) 0.04 ; not a BLOCK(1) SAME")
  expect_null(.extract_nm_same(ctl))
  expect_null(.extract_nm_same(nm_path("1cpt_oral.ctl")))
})

test_that(".extract_nm_same keeps two SAME sets apart", {
  f <- file.path(tmp_ctl_dir(), "two.ctl")
  writeLines(c("$PROBLEM x", "$INPUT ID OCC", "$OMEGA 0.1",
               "$OMEGA BLOCK(1) 0.04", "$OMEGA BLOCK(1) SAME",
               "$OMEGA BLOCK(1) 0.05", "$OMEGA BLOCK(1) SAME",
               "$SIGMA 0.1"), f)
  g <- .extract_nm_same(f)
  expect_equal(lapply(g$groups, function(x) x$etas), list(2:3, 4:5))
})

test_that(".extract_nm_same reports a SAME it cannot attribute", {
  f <- file.path(tmp_ctl_dir(), "amb.ctl")
  writeLines(c("$PROBLEM x", "$OMEGA 0.1 0.2", "$OMEGA BLOCK(1) SAME"), f)
  expect_match(.extract_nm_same(f)$problems, "multi-eta diagonal")
  writeLines(c("$PROBLEM x", "$OMEGA BLOCK(1) 0.1", "$OMEGA BLOCK(2) SAME"), f)
  expect_match(.extract_nm_same(f)$problems, "BLOCK\\(2\\) SAME")
})

test_that(".nm_input_items drops DROPped items and keeps both halves of a synonym", {
  f <- file.path(tmp_ctl_dir(), "inp.ctl")
  writeLines(c("$PROBLEM x", "$INPUT ID TIME OCC=DROP DV=CONC",
               "  VISIT = OCC2", "$DATA d.csv"), f)
  expect_setequal(.nm_input_items(f), c("ID", "TIME", "DV", "CONC", "VISIT", "OCC2"))
})

# -- the override rewriter ----------------------------------------------------

test_that(".same_iov_rewrite collapses the overrides and keeps the first eta", {
  lst <- list(quote(cl <- theta1 * exp(eta1)),
              quote(if (OCC == 1) cl <- theta1 * exp(eta1 + eta2)),
              quote(if (OCC == 2) cl <- theta1 * exp(eta1 + eta3)),
              quote(v <- theta2))
  rw <- .same_iov_rewrite(lst, c("eta2", "eta3"))
  expect_true(rw$ok)
  expect_equal(rw$column, "OCC")
  expect_length(rw$lst, 2L)
  expect_equal(deparse(rw$lst[[1]]), "cl <- theta1 * exp(eta1 + eta2)")
  expect_equal(deparse(rw$lst[[2]]), "v <- theta2")
})

test_that(".same_iov_rewrite accepts the intermediate-variable spelling and `n == COL`", {
  lst <- list(quote(rxm.bov <- 0),
              quote(if (1 == OCC) rxm.bov <- eta2),
              quote(if (OCC == 2) { rxm.bov <- eta3 }),
              quote(cl <- theta1 * exp(eta1 + rxm.bov)))
  rw <- .same_iov_rewrite(lst, c("eta2", "eta3"))
  expect_true(rw$ok)
  expect_equal(deparse(rw$lst[[1]]), "rxm.bov <- eta2")
})

test_that(".same_iov_rewrite works without a default assignment", {
  lst <- list(quote(if (OCC == 1) cl <- theta1 * exp(eta2)),
              quote(if (OCC == 2) cl <- theta1 * exp(eta3)))
  expect_true(.same_iov_rewrite(lst, c("eta2", "eta3"))$ok)
})

# Each refusal below is a shape that would change the model if rewritten.
test_that(".same_iov_rewrite refuses every shape it cannot prove equivalent", {
  ov <- function(k, rhs) as.call(list(as.name("if"), call("==", as.name("OCC"), k),
                                      call("<-", as.name("cl"), rhs)))
  rhs1 <- quote(theta1 * exp(eta1 + eta2))
  rhs2 <- quote(theta1 * exp(eta1 + eta3))
  base <- list(quote(cl <- theta1 * exp(eta1)), ov(1, rhs1), ov(2, rhs2))
  expect_true(.same_iov_rewrite(base, c("eta2", "eta3"))$ok)

  # the default is not the override with the eta removed
  bad <- base; bad[[1]] <- quote(cl <- theta1 * 2 * exp(eta1))
  expect_match(.same_iov_rewrite(bad, c("eta2", "eta3"))$reason, "above the overrides")
  # the overrides differ in more than the eta
  bad <- base; bad[[3]] <- ov(2, quote(theta1 * exp(eta1 - eta3)))
  expect_match(.same_iov_rewrite(bad, c("eta2", "eta3"))$reason, "more than the eta")
  # one occasion has no override
  expect_match(.same_iov_rewrite(base[1:2], c("eta2", "eta3"))$reason, "each of the")
  # an override tests a second column
  bad <- base; bad[[3]] <- as.call(list(as.name("if"), quote(VISIT == 2),
                                        call("<-", quote(cl), rhs2)))
  expect_match(.same_iov_rewrite(bad, c("eta2", "eta3"))$reason, "more than one column")
  # two overrides test the same occasion
  bad <- base; bad[[3]] <- ov(1, rhs2)
  expect_match(.same_iov_rewrite(bad, c("eta2", "eta3"))$reason, "same occasion")
  # an override assigns a different variable
  bad <- base; bad[[3]] <- as.call(list(as.name("if"), quote(OCC == 2),
                                        call("<-", quote(v), rhs2)))
  expect_match(.same_iov_rewrite(bad, c("eta2", "eta3"))$reason, "more than one variable")
  # the eta is read outside the overrides
  bad <- c(base, list(quote(v <- theta2 * exp(eta2))))
  expect_false(.same_iov_rewrite(bad, c("eta2", "eta3"))$ok)
  # a statement sits between two overrides
  bad <- list(base[[1]], base[[2]], quote(w <- 1), base[[3]])
  expect_match(.same_iov_rewrite(bad, c("eta2", "eta3"))$reason, "contiguous")
  # an else branch
  bad <- base; bad[[3]] <- quote(if (OCC == 2) cl <- theta1 * exp(eta1 + eta3) else cl <- 1)
  expect_false(.same_iov_rewrite(bad, c("eta2", "eta3"))$ok)
  # a column that is not a data item
  expect_match(.same_iov_rewrite(base, c("eta2", "eta3"), c("ID", "TIME"))$reason,
               "not a data item")
})

test_that(".same_value separates expressions that deparse differently", {
  expect_true(.same_value(quote(a * exp(b + 0)), quote(a * exp(b))))
  expect_false(.same_value(quote(a * exp(b)), quote(a * exp(b) * 2)))
  expect_false(.same_value(quote(a * undefined_fn(b)), quote(a * undefined_fn(b))))
})

# -- the pipeline -------------------------------------------------------------

test_that("a BLOCK(1) SAME set translates to one kappa and an iov_column", {
  r <- .same_run(nm_path("iov_same.ctl"))
  txt <- r$ferx_text
  expect_match(txt, "kappa IOV_CL ~ 0.04", fixed = TRUE)
  expect_match(txt, "iov_column = OCC", fixed = TRUE)
  expect_match(txt, "CL = TVCL * exp(ETA_CL + IOV_CL)", fixed = TRUE)
  expect_false(grepl("FIX", txt, fixed = TRUE))
  expect_false(grepl("if (OCC", txt, fixed = TRUE))
  # exactly one IIV omega is left; the repeated etas are not declared at all
  expect_length(.omega_lines(r), 1L)
  expect_length(r$unsupported, 0L)
  expect_false(any(grepl("^(WARN|ERROR)", r$warnings)))
  expect_equal(sum(grepl("SAME set", r$warnings)), 1L)
  # the label starts with IOV, which the ETA-coded-IOV heuristic would flag if
  # the eta were still in the omega block
  expect_false(any(grepl("looks like inter-occasion", r$warnings, fixed = TRUE)))
})

test_that("the same stream with independent blocks stays three omegas and no kappa", {
  r <- .same_run(.same_ctl("$OMEGA BLOCK(1) SAME" = "$OMEGA BLOCK(1) 0.04"))
  expect_length(.omega_lines(r), 4L)
  expect_false(grepl("(?m)^\\s*kappa ", r$ferx_text, perl = TRUE))
  expect_false(grepl("iov_column =", r$ferx_text, fixed = TRUE))
  expect_false(grepl("FIX", r$ferx_text, fixed = TRUE))
})

test_that("FIX on the base block carries to the kappa", {
  r <- .same_run(.same_ctl("$OMEGA BLOCK(1)" = "$OMEGA BLOCK(1) FIX"))
  expect_match(r$ferx_text, "kappa IOV_CL ~ 0.04 FIX", fixed = TRUE)
})

test_that("two SAME sets become two kappas on one iov_column", {
  ctl <- .same_ctl(
    "V  = THETA(2)" = c("  V = THETA(2) * EXP(ETA(5))",
                        "  IF (OCC.EQ.1) V = THETA(2) * EXP(ETA(5) + ETA(6))",
                        "  IF (OCC.EQ.2) V = THETA(2) * EXP(ETA(5) + ETA(7))",
                        "  IF (OCC.EQ.3) V = THETA(2) * EXP(ETA(5) + ETA(8))"))
  txt <- readLines(ctl)
  txt <- sub("^\\$SIGMA", paste0("$OMEGA 0.05 ; ETA_V\n$OMEGA BLOCK(1) 0.03 ; IOV_V\n",
                                 "$OMEGA BLOCK(1) SAME\n$OMEGA BLOCK(1) SAME\n$SIGMA"),
             txt)
  writeLines(txt, ctl)
  r <- .same_run(ctl)
  expect_length(r$unsupported, 0L)
  expect_match(r$ferx_text, "kappa IOV_CL ~ 0.04", fixed = TRUE)
  expect_match(r$ferx_text, "kappa IOV_V ~ 0.03", fixed = TRUE)
  expect_equal(length(gregexpr("iov_column =", r$ferx_text, fixed = TRUE)[[1]]), 1L)
  expect_match(r$ferx_text, "V = TVV * exp(ETA_V + IOV_V)", fixed = TRUE)
})

test_that("a SAME set whose overrides cannot be proved equivalent is an ERROR, not a quiet FIX", {
  # the third occasion's override uses a different formula
  r <- .same_run(.same_ctl(
    "IF (OCC.EQ.3) CL = THETA(1) * EXP(ETA(1) + ETA(4))" =
      "  IF (OCC.EQ.3) CL = THETA(1) * EXP(ETA(1) - ETA(4))"))
  expect_true(any(grepl("^ERROR \\| \\$OMEGA BLOCK SAME", r$warnings)))
  expect_true(any(grepl("^\\$OMEGA BLOCK SAME", r$unsupported)))
  expect_false(grepl("(?m)^\\s*kappa ", r$ferx_text, perl = TRUE))
  # the unchanged output is still what nonmem2rx produced
  expect_equal(sum(grepl("0.04 FIX$", .omega_lines(r))), 2L)
})

test_that("a SAME set on a column that is not in $INPUT is an ERROR", {
  r <- .same_run(.same_ctl("$INPUT  ID TIME AMT DV OCC" = "$INPUT  ID TIME AMT DV OCC=DROP"))
  expect_true(any(grepl("not a data item", r$warnings)))
  expect_false(grepl("iov_column =", r$ferx_text, fixed = TRUE))
})

# nonmem2rx refuses these streams itself, so the pipeline never reaches the
# reader's `problems`; the hint is handed over directly to pin what happens if
# it ever does.
test_that("a SAME record the reader could not attribute is an ERROR, not a guess", {
  ini <- .mk_ini(list(neta1 = 1L, neta2 = 1L, name = "eta1", est = 0.1, fix = FALSE,
                      condition = "id"))
  out <- .apply_same_iov(ini, list(), list(
    groups = list(), n_eta = 1L,
    problems = "a SAME record after a multi-eta diagonal $OMEGA record"))
  expect_match(out$warnings, "^ERROR \\| \\$OMEGA BLOCK SAME")
  expect_match(out$unsupported, "multi-eta diagonal")
  expect_identical(out$ini, ini)
})

test_that("BLOCK(2) SAME is an ERROR: ferx kappas are diagonal", {
  f <- file.path(tmp_ctl_dir(), "b2.ctl")
  writeLines(c(
    "$PROBLEM block 2 same", "$INPUT ID TIME AMT DV OCC", "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN1 TRANS2", "$PK",
    "CL = THETA(1) * EXP(ETA(1) + ETA(2))", "V = THETA(2) * EXP(ETA(3))", "S1 = V",
    "$ERROR", "IPRED = F", "Y = IPRED * (1 + EPS(1))",
    "$THETA (0, 3) (0, 30)", "$OMEGA BLOCK(2) 0.04 0.01 0.04", "$OMEGA BLOCK(2) SAME",
    "$SIGMA 0.02", "$ESTIMATION METHOD=1 INTER"), f)
  r <- .same_run(f)
  expect_true(any(grepl("^ERROR \\| \\$OMEGA BLOCK SAME", r$warnings)))
  expect_true(any(grepl("covariances", r$unsupported)))
})

test_that("an eta count that disagrees with the parsed model is an ERROR", {
  ini <- .mk_ini(list(neta1 = 1L, neta2 = 1L, name = "eta1", est = 0.1, fix = FALSE,
                      condition = "id"))
  out <- .apply_same_iov(ini, list(), list(groups = list(), n_eta = 3L, problems = character()))
  expect_match(out$warnings, "cannot be matched to etas")
  expect_length(out$unsupported, 1L)
})

test_that("a model with no SAME record is untouched", {
  expect_null(.extract_nm_same(nm_path("2cpt_iv.ctl")))
  ini <- .mk_ini(list(neta1 = 1L, neta2 = 1L, name = "eta1", est = 0.1, fix = FALSE,
                      condition = "id"))
  out <- .apply_same_iov(ini, list(quote(a <- 1)), NULL)
  expect_identical(out$ini, ini)
  expect_length(out$warnings, 0L)
})
