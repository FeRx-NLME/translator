# translator#74: compartment-free nlmixr2 models, from Bill Denney's standalone
# translator in FeRx-NLME/ferx-r#721. One block per row of the issue's table.
# The engine-level checks -- ferx's OFV at nlmixr2's estimates, both weightings
# -- are Tier 4 and live in test-concordance.R.
#
# Every test here was checked against a broken implementation before it was
# committed (CLAUDE.md, "a fixture picked to show a behaviour must be able to
# show that behaviour being wrong"); the break each one catches is named in it.

skip_if_not_installed("rxode2")

.cf <- function(f, ...) suppressWarnings(suppressMessages(
  nlmixr2_to_ferx(f, validate = FALSE, ...)))

# The [section] of an emitted file, as its body lines without indentation.
.section <- function(txt, name) {
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  start <- which(lines == paste0("[", name, "]"))
  if (length(start) == 0L) return(character())
  rest <- lines[-seq_len(start)]
  end  <- which(grepl("^\\[", rest))[1]
  body <- if (is.na(end)) rest else rest[seq_len(end - 1L)]
  trimws(body[nzchar(trimws(body))])
}

# Bill's example model, with the data columns DOSE, W and OCC.
mbma <- function() {
  ini({
    tve0 <- 9
    tvemax <- 6
    led50 <- log(30)
    lt50 <- log(2)
    eta.e0 ~ 1
    iov.emax ~ 0.05 | OCC
    sdw <- 1
  })
  model({
    E0 <- tve0 + eta.e0
    EMAX <- tvemax * exp(iov.emax)
    ONSET <- 0
    if (TIME > 0) ONSET <- TIME / (exp(lt50) + TIME)
    EFF <- EMAX * DOSE / (exp(led50) + DOSE) * ONSET
    pred <- (E0 - EFF) * sqrt(W)
    pred ~ add(sdw)
  })
}

# -- row 1: compartment-free model --------------------------------------------
# Break caught: without .lower_compartment_free() there is no [structural_model]
# at all, `PRED = ...` sits in [individual_parameters] and $unsupported reports
# "structural model".

test_that("row 1: a compartment-free prediction becomes `y =` with no pk/ode line", {
  res <- .cf(mbma)
  expect_identical(.section(res$ferx_text, "structural_model"),
                   "y = (E0 - EFF) * sqrt(W)")
  expect_false(grepl("^\\s*(pk |ode\\()", res$ferx_text))
  expect_false(any(startsWith(.section(res$ferx_text, "individual_parameters"), "PRED")))
  expect_length(res$unsupported, 0)
  expect_false(any(grepl("No structural model", res$warnings)))
})

test_that("row 1: a prediction reading a kappa stays an individual parameter", {
  # ferx refuses a kappa in a Form C output expression; measured on ferx-r
  # df7ffa7: "Form C output expressions cannot reference the IOV parameter".
  # Break caught: inlining unconditionally emits `y = (TVE0 + IOV_E) ...`.
  f <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; iov.e ~ 0.05 | OCC; sdw <- 1 })
    model({ E0 <- tve0 + eta.e0; pred <- (E0 + iov.e) * sqrt(W); pred ~ add(sdw) })
  }
  res <- .cf(f)
  expect_identical(.section(res$ferx_text, "structural_model"), "y = PRED")
  expect_true("PRED = (E0 + IOV_E) * sqrt(W)" %in%
                .section(res$ferx_text, "individual_parameters"))
})

test_that("row 1: a prediction assigned twice or read elsewhere is not inlined", {
  # Break caught: inlining the LAST assignment drops the first, so
  # `y = PRED * 2` would read an undefined PRED.
  f <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; sdw <- 1 })
    model({ pred <- tve0 + eta.e0; pred <- pred * 2; pred ~ add(sdw) })
  }
  res <- .cf(f)
  expect_identical(.section(res$ferx_text, "structural_model"), "y = PRED")
  ip <- .section(res$ferx_text, "individual_parameters")
  expect_identical(ip, c("PRED = TVE0 + ETA_E0", "PRED = PRED * 2"))

  # Break caught: ignoring the `if` arm inlines the default, so the block that
  # overrides it assigns a PRED nothing reads. ferx validates that file
  # without complaint, so nothing downstream would notice (review #75, #5).
  g <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; sdw <- 1 })
    model({ pred <- tve0 + eta.e0; if (TIME > 0) pred <- tve0 * 2; pred ~ add(sdw) })
  }
  res <- .cf(g)
  expect_identical(.section(res$ferx_text, "structural_model"), "y = PRED")
  expect_true("PRED = TVE0 + ETA_E0" %in% .section(res$ferx_text, "individual_parameters"))
})

test_that("row 1: a prediction is not inlined past a reassignment of what it reads", {
  # Break caught (review #75, #2): inlining moves the prediction to the end of
  # the model, so `y = E0 - EFF` read the later `EFF = 0` and predicted E0.
  f <- function() {
    ini({ e0 <- 9; temax <- 2; eta.e ~ 1; sdw <- 1 })
    model({
      E0 <- e0 + eta.e; eff <- temax * DOSE; pred <- E0 - eff; eff <- 0
      pred ~ add(sdw)
    })
  }
  res <- .cf(f)
  expect_identical(.section(res$ferx_text, "structural_model"), "y = PRED")
  expect_identical(.section(res$ferx_text, "individual_parameters"),
                   c("E0 = TVE0 + ETA_E", "EFF = TEMAX * DOSE", "PRED = E0 - EFF",
                     "EFF = 0"))
  # A reassignment BEFORE the prediction does not stop the inlining.
  g <- function() {
    ini({ e0 <- 9; temax <- 2; eta.e ~ 1; sdw <- 1 })
    model({
      E0 <- e0 + eta.e; eff <- 0; eff <- temax * DOSE; pred <- E0 - eff
      pred ~ add(sdw)
    })
  }
  expect_identical(.section(.cf(g)$ferx_text, "structural_model"), "y = E0 - EFF")
})

# -- row 2: one-line if after a default ---------------------------------------
# Break caught: the default must come first. ferx-core#1857 refuses an
# else-less `if` whose variable has no earlier assignment, and emitting the
# block before its default (e.g. by sorting conditionals after assignments the
# other way) gives exactly that.

test_that("row 2: a one-line if is a block, after its default", {
  ip <- .section(.cf(mbma)$ferx_text, "individual_parameters")
  i_default <- which(ip == "ONSET = 0")
  i_if      <- which(ip == "if (TIME > 0) { ONSET = TIME/(exp(LT50) + TIME) }")
  expect_length(i_default, 1)
  expect_length(i_if, 1)
  expect_lt(i_default, i_if)
  # ...and before the statement that reads it.
  expect_lt(i_if, which(startsWith(ip, "EFF =")))
})

# -- row 3: iov ~ v | OCC ------------------------------------------------------
# Break caught: read as IIV it would be `omega IOV_EMAX` with no iov_column.
# The engine-level counterpart (a kappa emitted as an omega changes the OFV) is
# in test-concordance.R.

test_that("row 3: `iov ~ v | OCC` is a kappa with iov_column = OCC", {
  txt <- .cf(mbma)$ferx_text
  par <- .section(txt, "parameters")
  expect_true("kappa IOV_EMAX ~ 0.05" %in% par)
  expect_false(any(startsWith(par, "omega IOV_EMAX")))
  expect_true("iov_column = OCC" %in% .section(txt, "fit_options"))
  expect_true("EMAX = TVEMAX * exp(IOV_EMAX)" %in% .section(txt, "individual_parameters"))
})

# -- row 4: inverse-variance weighting -----------------------------------------

test_that("row 4: as written, the prediction keeps its sqrt(W)", {
  res <- .cf(mbma)
  expect_identical(.section(res$ferx_text, "error_model"), "DV ~ additive(SDW)")
  expect_identical(.section(res$ferx_text, "covariates"), character())
  expect_identical(res$weight_column, NA_character_)
})

test_that("row 4: weighting = 'weight' strips sqrt(W) and weights the error", {
  # Break caught: leaving `* sqrt(W)` in `y` double-weights (the OFV test in
  # test-concordance.R then misses by thousands); omitting [covariates] makes
  # ferx refuse the weight ("references undeclared covariate `W`").
  res <- .cf(mbma, weighting = "weight")
  expect_identical(.section(res$ferx_text, "structural_model"), "y = (E0 - EFF)")
  expect_identical(.section(res$ferx_text, "error_model"),
                   "DV ~ additive(SDW) weight = 1/sqrt(W)")
  # Every covariate, not only W: ferx warns (W_PARSE) about any a present
  # [covariates] block leaves out.
  expect_identical(.section(res$ferx_text, "covariates"),
                   c("DOSE continuous", "W continuous"))
  expect_identical(res$weight_column, "W")
})

test_that("row 4: sqrt(W) on the left of the product is recognised too", {
  f <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; sdw <- 1 })
    model({ E0 <- tve0 + eta.e0; pred <- sqrt(W) * E0; pred ~ add(sdw) })
  }
  res <- .cf(f, weighting = "weight")
  expect_identical(.section(res$ferx_text, "structural_model"), "y = E0")
})

test_that("row 4: weighting = 'weight' refuses what it cannot represent", {
  # Each of these would otherwise be scored as DV/sqrt(W) against an unchanged
  # prediction. A silent fallback to "as_written" is the break caught.
  no_sqrt <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; sdw <- 1 })
    model({ pred <- tve0 + eta.e0; pred ~ add(sdw) })
  }
  prop <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; pw <- 0.1 })
    model({ E0 <- tve0 + eta.e0; pred <- E0 * sqrt(W); pred ~ prop(pw) })
  }
  w_param <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; sdw <- 1 })
    model({ W <- 4; E0 <- tve0 + eta.e0; pred <- E0 * sqrt(W); pred ~ add(sdw) })
  }
  lincmt <- function() {
    ini({ tcl <- 1; tv <- 10; eta.cl ~ 0.1; add.sd <- 0.1 })
    model({ cl <- tcl * exp(eta.cl); v <- tv; linCmt() ~ add(add.sd) })
  }
  expect_error(.cf(no_sqrt, weighting = "weight"), "f \\* sqrt\\(W\\)")
  expect_error(.cf(prop, weighting = "weight"), "single additive")
  expect_error(.cf(w_param, weighting = "weight"), "not a data column")
  # A W assigned only inside `if` arms is a model variable, not a column: it
  # was declared in [covariates] and ferx_weighted_data() then aborted on a
  # column that does not exist (review #75, #4).
  w_in_if <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; sdw <- 1 })
    model({
      E0 <- tve0 + eta.e0
      if (TIME > 0) { WW <- 2 } else { WW <- 1 }
      pred <- E0 * sqrt(WW); pred ~ add(sdw)
    })
  }
  # A solver builtin is not a column: `sqrt(TIME)` declared TIME continuous
  # (review #75, #6).
  w_builtin <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; sdw <- 1 })
    model({ E0 <- tve0 + eta.e0; pred <- E0 * sqrt(TIME); pred ~ add(sdw) })
  }
  reassigned <- function() {
    ini({ tve0 <- 9; eta.e0 ~ 1; sdw <- 1 })
    model({ E0 <- tve0 + eta.e0; pred <- E0 * sqrt(W); E0 <- 0; pred ~ add(sdw) })
  }
  expect_error(.cf(w_in_if, weighting = "weight"), "'WW' in sqrt\\(WW\\) is a model")
  expect_error(.cf(w_builtin, weighting = "weight"), "'TIME' .* is a ferx builtin")
  expect_error(.cf(reassigned, weighting = "weight"), "reassigned after")
  expect_error(.cf(lincmt, weighting = "weight"), "compartment-free")
  expect_error(to_ferx(nm_path("1cpt_oral.ctl"), "nonmem", validate = FALSE,
                       weighting = "weight"), "only available for nlmixr2")
})

test_that("row 4: ferx_weighted_data divides DV back by sqrt(W) on observations", {
  d <- data.frame(ID = 1, TIME = 0:2, DV = c(0, 6, 9), EVID = c(1, 0, 0),
                  W = c(NA, 4, 9))
  out <- ferx_weighted_data(d, "W")
  expect_identical(out$DV, c(0, 3, 3))
  expect_identical(out$W, d$W)
  bad <- d; bad$W[2] <- 0
  expect_error(ferx_weighted_data(bad, "W"), "positive and finite")
  expect_error(ferx_weighted_data(d, "SE"), "must name a column")

  # No EVID: observation rows are found from MDV, then from AMT, so a dose row
  # with W missing is neither checked nor divided (review #75, #7).
  nm <- data.frame(ID = 1, TIME = 0:2, DV = c(0, 6, 9), AMT = c(100, 0, 0),
                   MDV = c(1, 0, 0), W = c(NA, 4, 9))
  expect_identical(ferx_weighted_data(nm, "W")$DV, c(0, 3, 3))
  expect_identical(ferx_weighted_data(nm[, names(nm) != "MDV"], "W")$DV, c(0, 3, 3))
  # MDV wins over AMT when both are present: an MDV = 1 row with AMT = 0 (a
  # reset or other non-observation record) is not divided.
  other <- nm; other$AMT[1] <- 0
  expect_identical(ferx_weighted_data(other, "W")$DV, c(0, 3, 3))
})

# -- row 5: constant expressions and fix() in ini() -----------------------------
# Break caught: emitting the expression text (`log(30)`) is a parse error in a
# theta initial, and dropping FIX turns a held parameter into an estimated one.

test_that("row 5: a constant expression is evaluated and fix() is FIX", {
  f <- function() {
    ini({
      led50 <- log(30)
      lt50 <- fix(log(2))
      eta.e ~ fix(0.5)
      sdw <- fix(1)
    })
    model({ E <- exp(led50) + exp(lt50) + eta.e; pred <- E; pred ~ add(sdw) })
  }
  par <- .section(.cf(f)$ferx_text, "parameters")
  th <- par[startsWith(par, "theta LED50(")]
  expect_equal(as.numeric(sub("^theta LED50\\(([^,]+),.*$", "\\1", th)), log(30),
               tolerance = 1e-14)
  expect_true(any(grepl("^theta LT50\\(0\\.69314718055994[0-9]*, .*, FIX\\)$", par)))
  expect_true("omega ETA_E ~ 0.5 FIX" %in% par)
  expect_true("sigma SDW ~ 1.0 FIX (sd)" %in% par)
})

# -- row 6: names ferx cannot tell apart ---------------------------------------
# Break caught: without .source_name_collisions() each of these merges two
# source names into one emitted name with no diagnostic. The `dose` case is the
# worst: it emits `DOSE = 3` and then reads it where the source read the data
# column.

test_that("row 6: names differing only in case or dots are refused by name", {
  case_vars <- function() {
    ini({ e0 <- 9; eta.e ~ 1; sdw <- 1 })
    model({ eff <- e0 + eta.e; EFF <- eff * 2; pred <- EFF; pred ~ add(sdw) })
  }
  dotted <- function() {
    ini({ e0 <- 9; eta.e ~ 1; sdw <- 1 })
    model({ a.b <- e0 + eta.e; a_b <- a.b * 2; pred <- a_b; pred ~ add(sdw) })
  }
  shadows_column <- function() {
    ini({ e0 <- 9; eta.e ~ 1; sdw <- 1 })
    model({ dose <- 3; E <- e0 + eta.e + DOSE * dose; pred <- E; pred ~ add(sdw) })
  }
  r1 <- .cf(case_vars)
  expect_true(any(grepl("^ERROR .*'eff' \\(model variable\\) and 'EFF' \\(model variable\\)",
                        r1$warnings)))
  expect_true(any(grepl("eff, EFF", r1$unsupported)))
  r2 <- .cf(dotted)
  expect_true(any(grepl("'a.b' \\(model variable\\) and 'a_b' \\(model variable\\)",
                        r2$warnings, fixed = FALSE)))
  r3 <- .cf(shadows_column)
  expect_true(any(grepl("'dose' \\(model variable\\) and 'DOSE' \\(data column\\)",
                        r3$warnings)))
  # The refusal reaches the emitted file, where the merged name appears.
  expect_match(r3$ferx_text, "# WARNING: 'dose' \\(model variable\\) and 'DOSE'")
})

test_that("row 6: collisions an owner already resolves are renamed, not refused", {
  # A theta beside the variable it defines is .deshadow_theta_names()'s job,
  # which renames the theta -- refusing it here would regress every source
  # written `cl <- cl * exp(eta.cl)`. Reassigning one name is not a collision.
  f <- function() {
    ini({ e0 <- 9; eta.e ~ 1; sdw <- 1 })
    model({ E0 <- e0 + eta.e; E0 <- E0 * 2; pred <- E0; pred ~ add(sdw) })
  }
  res <- .cf(f)
  expect_length(res$unsupported, 0)
  expect_true(any(grepl("^theta TVE0\\(", .section(res$ferx_text, "parameters"))))
  expect_false(any(grepl("the same name to ferx", res$warnings)))
})

test_that("row 6: a theta or eta under a data column's header is renamed", {
  # Break caught (review #75, #1): ferx resolves theta and eta before a
  # covariate, so `theta W` beside the column `W` made `sqrt(W)` read the
  # theta (OFV 1364.8 against 2206.4 on the engine), and `omega DOSE` beside
  # the column `DOSE` read the eta, both with empty $unsupported.
  th <- function() {
    ini({ e <- 1; w <- 2; eta.e ~ 0.1; sd <- 0.5 })
    model({ E <- e + eta.e; G <- E * w; pred <- G * sqrt(W); pred ~ add(sd) })
  }
  res <- .cf(th)
  expect_true(any(grepl("^theta TVW\\(", .section(res$ferx_text, "parameters"))))
  expect_identical(.section(res$ferx_text, "individual_parameters"),
                   c("E = TVE + ETA_E", "G = E * TVW"))
  expect_identical(.section(res$ferx_text, "structural_model"), "y = G * sqrt(W)")
  expect_true(any(grepl("^WARN .*theta 'W' has the name of the data column 'W'",
                        res$warnings)))
  et <- function() {
    ini({ e <- 1; dose ~ 0.1; sd <- 0.5 })
    model({ G <- e + dose + DOSE; pred <- G; pred ~ add(sd) })
  }
  res <- .cf(et)
  expect_true("omega DOSE_1 ~ 0.1" %in% .section(res$ferx_text, "parameters"))
  expect_identical(.section(res$ferx_text, "individual_parameters"), "G = E + DOSE_1 + DOSE")

  # Exact header, not case-folded: ferx binds a column case-sensitively, and a
  # theta `W` leaves a column `w` readable (OFV identical either way), so it is
  # not renamed.
  lower <- function() {
    ini({ e <- 1; W <- 2; eta.e ~ 0.1; sd <- 0.5 })
    model({ E <- e + eta.e; G <- E * W; pred <- G * sqrt(w); pred ~ add(sd) })
  }
  expect_true(any(grepl("^theta W\\(", .section(.cf(lower)$ferx_text, "parameters"))))
})

test_that("row 6: data columns differing only in case are two columns", {
  # Break caught (review #75, #3): the key upper-cased columns too, so `wt`
  # and `WT` were refused, while ferx reads them as two columns.
  f <- function() {
    ini({ e0 <- 1; eta.e ~ 0.1; sd <- 0.5 })
    model({ E0 <- e0 + eta.e + wt + WT; pred <- E0; pred ~ add(sd) })
  }
  res <- .cf(f)
  expect_length(res$unsupported, 0)
  expect_identical(.section(res$ferx_text, "individual_parameters"), "E0 = TVE0 + ETA_E + wt + WT")

  # The same rule decides a variable against a column: `dose` emits `DOSE`,
  # which a column `Dose` does not match. Measured on ferx-r main
  # (maxiter = 0): a variable `W` beside a column `w` gives the OFV of the
  # renamed model (1035.56); only beside a column `W` does it shadow (995.66).
  g <- function() {
    ini({ e0 <- 1; eta.e ~ 0.1; sd <- 0.5 })
    model({ dose <- 3; E0 <- e0 + eta.e + Dose * dose; pred <- E0; pred ~ add(sd) })
  }
  expect_length(.cf(g)$unsupported, 0)
})

test_that("row 6: an ODE state is left to the state sanitiser, not refused", {
  # Break caught (review #75, #12): a state read on a right-hand side was taken
  # for a data column, so a state `CENTR` beside a variable `centr` was refused
  # as "(data column)" although .sanitise_state_names() already renames the
  # state.
  f <- function() {
    ini({ tk <- 0.1; eta.k ~ 0.1; sd <- 0.5 })
    model({
      k <- tk * exp(eta.k); centr <- 2 * k
      d/dt(CENTR) <- -centr * CENTR
      cp <- CENTR; cp ~ add(sd)
    })
  }
  res <- .cf(f)
  expect_false(any(grepl("the same name to ferx", res$warnings)))
  expect_false(any(grepl("cannot tell apart", res$unsupported)))
})

test_that("row 6: a reserved dose-attribute name is renamed with a warning", {
  # `f` would be bioavailability in ferx, silently (CLAUDE.md, reserved names).
  f <- function() {
    ini({ e0 <- 9; eta.e ~ 1; sdw <- 1 })
    model({ f <- 2; E <- e0 + eta.e + f; pred <- E; pred ~ add(sdw) })
  }
  res <- .cf(f)
  ip <- .section(res$ferx_text, "individual_parameters")
  expect_true("F_PAR = 2" %in% ip)
  expect_false(any(grepl("^F =", ip)))
})

test_that("row 6: the collision check raises nothing on the bundled NONMEM corpus", {
  # NONMEM is case-insensitive and nonmem2rx canonicalises case, so every hit
  # here would be a false refusal. Break caught: a check that also counts the
  # raw iniDf names (`t.CL` beside the variable `cl`) fires on most of these.
  skip_if_not_installed("nonmem2rx")
  for (f in .bundled_nm_models()) {
    ui <- suppressWarnings(suppressMessages(nonmem2rx::nonmem2rx(f)))
    hit <- .source_name_collisions(ui$lstExpr, .norm_map_from_ini(ui$iniDf))
    expect_identical(hit$unsupported, character(), label = basename(f))
  }
})

# -- row 7: linear thetas are not log-packed -----------------------------------
# ferx log-packs a theta whose lower bound is >= 0 (ferx-core
# parameterization.rs, theta_packs_log). nlmixr2 estimates an unbounded theta
# on its own scale, so the emitted lower bound must stay negative. Measured on
# the example data: (-1e15, 1e15) and Bill Denney's (init - 100, init + 100)
# give identical free fits, so the existing unbounded box is kept rather than
# imposing a +-100 range nlmixr2 never had. The engine-level break -- a theta
# starting positive whose estimate is negative -- is in test-concordance.R.

test_that("row 7: an unbounded theta starting >= 0 keeps a negative lower bound", {
  par <- .section(.cf(mbma)$ferx_text, "parameters")
  for (nm in c("TVE0", "TVEMAX", "LED50", "LT50")) {
    th <- par[startsWith(par, paste0("theta ", nm, "("))]
    lo <- as.numeric(strsplit(sub("^theta [A-Z0-9_]+\\((.*)\\)$", "\\1", th), ", ")[[1]][2])
    expect_lt(lo, 0, label = paste(nm, "lower bound"))
  }
})

test_that("row 7: an explicit nlmixr2 lower bound is kept as written", {
  # The user declared it; ferx's log packing is its representation of that.
  f <- function() {
    ini({ e0 <- c(0, 9); eta.e ~ 1; sdw <- 1 })
    model({ E <- e0 + eta.e; pred <- E; pred ~ add(sdw) })
  }
  expect_true("theta E0(9.0, 0.0, 1e15)" %in% .section(.cf(f)$ferx_text, "parameters"))
})

# -- row 8: variances below ferx's guard ----------------------------------------

test_that("row 8: an estimated variance below the guard starts at 1e-5, a FIX one is kept", {
  # Break caught: without the floor ferx moves 1e-8 to its exp(-12) rail with
  # no diagnostic. Flooring the FIX one too is the opposite break: ferx holds a
  # FIX 1e-8 exactly (measured, OFV 7932.16 against 7929.92 at 1e-5).
  f <- function() {
    ini({
      e0 <- 9
      eta.a ~ 1e-8
      eta.b ~ fix(1e-8)
      eta.c ~ 7e-6
      iov.k ~ 1e-9 | OCC
      sdw <- 1
    })
    model({
      A <- e0 + eta.a + iov.k
      B <- A + eta.b
      E <- B + eta.c
      pred <- E
      pred ~ add(sdw)
    })
  }
  res <- .cf(f)
  par <- .section(res$ferx_text, "parameters")
  expect_true("omega ETA_A ~ 0.00001" %in% par)
  expect_true("omega ETA_B ~ 0.00000001 FIX" %in% par)
  expect_true("omega ETA_C ~ 0.000007" %in% par)   # above the 6.144e-6 guard
  expect_true("kappa IOV_K ~ 0.00001" %in% par)
  floored <- grep("below the smallest variance ferx can represent", res$warnings, value = TRUE)
  expect_length(floored, 2)
  expect_true(any(grepl("^WARN .*omega ETA_A", floored)))
  expect_true(any(grepl("^WARN .*kappa IOV_K", floored)))
})

test_that("row 8: the floor names an omega as the emitted file does", {
  # Break caught (review #75, #9): the floor ran before random-effect renaming,
  # so its WARN named `ETA_E` while the file floored `ETA_E_1`.
  f <- function() {
    ini({ e0 <- 9; eta.e ~ 0.1; ETA.E ~ 1e-8; sdw <- 1 })
    model({ E <- e0 + eta.e; E2 <- E + ETA.E; pred <- E2; pred ~ add(sdw) })
  }
  res <- .cf(f)
  expect_true("omega ETA_E_1 ~ 0.00001" %in% .section(res$ferx_text, "parameters"))
  floored <- grep("below the smallest variance", res$warnings, value = TRUE)
  expect_length(floored, 1)
  expect_match(floored, "omega ETA_E_1 starts at")
})

test_that("row 8: a block omega below the guard is reported, not moved", {
  oms <- list(list(type = "block", names = c("A", "B"),
                   values = c(1e-8, 0, 0.1), fixed = FALSE))
  out <- .floor_variances(oms, list())
  expect_identical(out$omegas, oms)
  expect_length(out$warnings, 1)
  expect_match(out$warnings, "block omega element A .* NOT moved")
})

# -- IR / emitter units ---------------------------------------------------------

test_that("an algebraic structural model needs exactly one y", {
  ir <- new_ferx_ir(structural = list(type = "algebraic", y = c("A", "B")))
  expect_error(validate_ferx_ir(ir), "single non-empty string")
  ir <- new_ferx_ir(structural = list(type = "algebraic", y = "A"),
                    scaling = list(y = "B"))
  expect_error(validate_ferx_ir(ir), "cannot carry")
})

test_that("a weight is refused on a proportional or per-CMT error model", {
  ir <- new_ferx_ir(error_model = list(list(dv = "DV", type = "proportional",
                                            params = "P", weight = "1/sqrt(W)")))
  expect_error(validate_ferx_ir(ir), "proportional")
  ir <- new_ferx_ir(error_model = list(
    list(dv = "DV", type = "additive", params = "A", cmt = 1L, weight = "W"),
    list(dv = "DV", type = "additive", params = "B", cmt = 2L)))
  expect_error(validate_ferx_ir(ir), "single-endpoint")
})

test_that("a fitted nlmixr2 object is translated from its finalUi", {
  # Break caught (review #75, #8): rxode2() on the fit itself returns no usable
  # iniDf. A stub with the fit's class stands in for nlmixr2est.
  ui  <- rxode2::rxode2(mbma)
  fit <- structure(list(finalUi = ui), class = c("nlmixr2FitCore", "list"))
  expect_identical(.cf(fit)$ferx_text, .cf(mbma)$ferx_text)
})

test_that("an IR with NULL covariates still emits", {
  # Break caught (review #75, #10): an IR built by hand, or saved before the
  # field existed, has no `covariates` and was refused by validate_ferx_ir().
  ir <- rxui_to_ir(rxode2::rxode2(mbma))
  ir$covariates <- NULL
  expect_match(emit_ferx(ir), "y = (E0 - EFF) * sqrt(W)", fixed = TRUE)
})

test_that("[covariates] is in the same order whatever the collation locale", {
  # Break caught (review #75, #11): sort() follows LC_COLLATE, which testthat
  # forces to C but a user session does not.
  f <- function() {
    ini({ e0 <- 9; eta.e ~ 1; sdw <- 1 })
    model({ E0 <- e0 + eta.e + alpha + Beta; pred <- E0 * sqrt(W); pred ~ add(sdw) })
  }
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
  set <- suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))
  skip_if_not(nzchar(set) && identical(sort(c("alpha", "Beta")), c("alpha", "Beta")),
              "no en_US.UTF-8 collation on this machine")
  expect_identical(.section(.cf(f, weighting = "weight")$ferx_text, "covariates"),
                   c("Beta continuous", "W continuous", "alpha continuous"))
})
