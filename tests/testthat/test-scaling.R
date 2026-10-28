# S<n> as an expression (issue #45). Every bundled model before this wrote
# `S<n> = <bare parameter>`, the one shape a first-identifier read gets right, so
# nothing else in the suite could see the right-hand side being truncated.

# One-compartment oral ODE whose observed compartment is 2 and whose readout
# divides the amount by S2. `s_line` is the statement under test.
.scale_ctl <- function(s_line, pk_extra = character(), theta_extra = character()) {
  ctl <- file.path(tmp_ctl_dir(), "m.ctl")
  writeLines(c(
    "$PROBLEM scale expression",
    "$INPUT ID TIME DV AMT EVID MDV CMT WT",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINES ADVAN6 TOL=6",
    "$MODEL", "COMP = (DEPOT)", "COMP = (CENT)",
    "$PK",
    "  KA = THETA(1)", "  CL = THETA(2)*EXP(ETA(1))", "  V  = THETA(3)",
    pk_extra, s_line,
    "$DES", "  DADT(1) = -KA*A(1)", "  DADT(2) = KA*A(1) - CL/V*A(2)",
    "$ERROR", "  Y = A(2)/S2*(1+EPS(1))",
    "$THETA", "  (0,1)", "  (0,3)", "  (0,30)", theta_extra,
    "$OMEGA 0.09", "$SIGMA 0.01", "$EST METHOD=1 INTER"), ctl)
  ctl
}

.scale_run <- function(...) suppressWarnings(suppressMessages(
  nm_to_ferx(.scale_ctl(...), validate = FALSE)))

.scaling_block <- function(txt) {
  ln <- strsplit(txt, "\n")[[1]]
  i  <- match("[scaling]", ln)
  if (is.na(i)) return(character())
  ln[i + 1L]
}

test_that(".extract_nm_scaling() returns the whole right-hand side", {
  mk <- function(...) {
    f <- file.path(tmp_ctl_dir(), "m.ctl")
    writeLines(c("$PROBLEM x", "$PK", ..., "$ERROR", "Y = F"), f)
    .extract_nm_scaling(f)
  }
  expect_equal(mk("S2 = V*MW")[["2"]], "V*MW")
  expect_equal(mk("S2 = MW*V")[["2"]], "MW*V")
  expect_equal(mk("S2 = V/1000")[["2"]], "V/1000")
  expect_equal(mk("S2 = 1000*V")[["2"]], "1000*V")
  expect_equal(mk("S1 = 1000")[["1"]], "1000")
  expect_equal(mk("  S2 = V * MW   ; unit conversion")[["2"]], "V * MW")
  # The shape every bundled model has stays what it was.
  expect_equal(mk("S2 = V")[["2"]], "V")
  # The product and the quotient are different scales.
  expect_false(identical(mk("S2 = V*MW")[["2"]], mk("S2 = V/MW")[["2"]]))
})

test_that(".resolve_scale_expr() keeps the operators and refuses what it cannot place", {
  th <- c("TVKA", "TVCL"); ip <- c("V", "MW", "KA")
  r <- .resolve_scale_expr
  expect_equal(r("V*MW", th, ip), "V * MW")
  expect_equal(r("V/MW", th, ip), "V/MW")
  expect_equal(r("(V+1)*MW", th, ip), "(V + 1) * MW")
  # Operand order is not semantics: the two spellings evaluate the same.
  env <- list(V = 2, MW = 5)
  expect_equal(eval(str2lang(r("V*MW", th, ip)), env),
               eval(str2lang(r("MW*V", th, ip)), env))
  # ...and a product is not a quotient.
  expect_false(isTRUE(all.equal(eval(str2lang(r("V*MW", th, ip)), env),
                                eval(str2lang(r("V/MW", th, ip)), env))))
  # NONMEM is case-insensitive, ferx names are not: resolve to the emitted name.
  expect_equal(r("v*mw", th, ip), "V * MW")
  expect_equal(r("TVCL*2", th, ip), "TVCL * 2")
  # A constant-only scale is a number, not an unresolvable variable.
  expect_identical(r("1000", th, ip), 1000)
  expect_identical(r("1/1000", th, ip), 0.001)
  # Everything below must be refused so the caller can report the whole statement.
  expect_null(r("V*WT/70", th, ip))      # covariate: not a theta or individual parameter
  expect_null(r("exp(V)", th, ip))       # function call
  expect_null(r("V*", th, ip))           # does not parse
  expect_null(r("0", th, ip))            # obs_scale must be strictly positive
  expect_null(r("-3", th, ip))
})

test_that("S<n> as an expression reaches [scaling] whole", {
  skip_if_not_installed("nonmem2rx")
  prod <- .scale_run("  S2 = V*MW", pk_extra = "  MW = 500")
  expect_equal(.scaling_block(prod$ferx_text), "  obs_scale = V * MW")
  # The INFO quotes the source statement, not a token of it.
  expect_true(any(grepl("S2 = V*MW detected -- emitting [scaling] obs_scale = V * MW",
                        prod$warnings, fixed = TRUE)))

  # Operand order and operator both change the emitted scale.
  expect_equal(.scaling_block(.scale_run("  S2 = MW*V", pk_extra = "  MW = 500")$ferx_text),
               "  obs_scale = MW * V")
  expect_equal(.scaling_block(.scale_run("  S2 = V/MW", pk_extra = "  MW = 500")$ferx_text),
               "  obs_scale = V/MW")
  expect_equal(.scaling_block(.scale_run("  S2 = V/1000")$ferx_text),
               "  obs_scale = V/1000")
  expect_equal(.scaling_block(.scale_run("  S2 = 1000*V")$ferx_text),
               "  obs_scale = 1000 * V")

  # A bare number is an ordinary numeric scale, not an unresolvable variable.
  num <- .scale_run("  S2 = 1000")
  expect_length(num$unsupported, 0L)
  expect_equal(trimws(sub("^.*=", "", .scaling_block(num$ferx_text))), "1000.0")

  # The shape every bundled model has is unchanged.
  expect_equal(.scaling_block(.scale_run("  S2 = V")$ferx_text), "  obs_scale = V")
})

test_that("an S<n> expression that cannot be placed is reported whole and loudly", {
  skip_if_not_installed("nonmem2rx")
  for (src in c("V*WT/70", "EXP(V)")) {
    r <- .scale_run(paste0("  S2 = ", src))
    expect_equal(.scaling_block(r$ferx_text), character(), label = src)
    # The whole statement, so the message names what the source wrote.
    expect_true(any(startsWith(r$warnings, "ERROR") &
                    grepl(paste0("S2 = ", src, " for the observed"), r$warnings,
                          fixed = TRUE)), label = src)
    expect_true(any(grepl(paste0("not resolvable: ", src), r$unsupported,
                          fixed = TRUE)), label = src)
  }
})

# S<n> on a solved system (issue #44). The pk macros already return A/v, so only
# what S<n> carries beyond the volume belongs in [scaling]. Before, the solved
# path never read S<n> at all.
.solved_ctl <- function(s_lines, adv = "ADVAN1 TRANS2", pk_extra = character(),
                        ka = FALSE) {
  ctl <- file.path(tmp_ctl_dir(), "s.ctl")
  writeLines(c(
    "$PROBLEM solved scale", "$INPUT ID TIME DV AMT EVID MDV CMT", "$DATA d.csv IGNORE=@",
    paste("$SUBROUTINES", adv),
    "$PK", if (ka) "  KA = THETA(3)", "  CL = THETA(1)*EXP(ETA(1))", "  V  = THETA(2)",
    pk_extra, s_lines,
    "$ERROR", "  Y = F*(1+EPS(1))",
    "$THETA", "  (0,3)", "  (0,30)", if (ka) "  (0,1)",
    "$OMEGA 0.09", "$SIGMA 0.01", "$EST METHOD=1 INTER"), ctl)
  ctl
}

.solved_run <- function(...) suppressWarnings(suppressMessages(
  nm_to_ferx(.solved_ctl(...), validate = FALSE)))

test_that("S<n> on a solved system emits only what is left after the volume", {
  skip_if_not_installed("nonmem2rx")
  blk <- function(s, ...) .scaling_block(.solved_run(s, ...)$ferx_text)
  # The volume divides out of these.
  expect_equal(blk("  S1 = V/1000"), "  obs_scale = 0.001")
  expect_equal(blk("  S1 = 1000*V"), "  obs_scale = 1000.0")
  expect_equal(blk("  S1 = V*MW", pk_extra = "  MW = 500"), "  obs_scale = MW")
  expect_equal(blk("  S1 = MW*V", pk_extra = "  MW = 500"), "  obs_scale = MW")
  expect_equal(blk("  S1 = V/MW", pk_extra = "  MW = 500"), "  obs_scale = 1/MW")
  # The shape every bundled model has: nothing is left, so nothing is emitted.
  expect_equal(blk("  S1 = V"), character())
  # The volume does not divide out: the ratio is emitted as written.
  expect_equal(blk("  S1 = 1000"), "  obs_scale = 1000/V")
  expect_equal(blk("  S1 = V + 1"), "  obs_scale = (V + 1)/V")

  r <- .solved_run("  S1 = V/1000")
  expect_true(any(grepl("S1 = V/1000 detected on a solved system", r$warnings,
                        fixed = TRUE)))
  expect_length(r$unsupported, 0L)
})

test_that("S<n> on a solved system is read for the central compartment only", {
  skip_if_not_installed("nonmem2rx")
  # Behind a depot the observed compartment is 2.
  expect_equal(.scaling_block(.solved_run("  S2 = V/1000", adv = "ADVAN2 TRANS2",
                                          ka = TRUE)$ferx_text),
               "  obs_scale = 0.001")
  # S1 scales the depot, which is not observed: no scaling, and it says so.
  r <- .solved_run("  S1 = V/1000", adv = "ADVAN2 TRANS2", ka = TRUE)
  expect_equal(.scaling_block(r$ferx_text), character())
  expect_true(any(startsWith(r$warnings, "WARN") &
                  grepl("but not for the observed compartment 2", r$warnings,
                        fixed = TRUE)))
})

test_that("an S<n> on a solved system that cannot be placed is reported whole", {
  skip_if_not_installed("nonmem2rx")
  r <- .solved_run("  S1 = V*WT/70")
  expect_equal(.scaling_block(r$ferx_text), character())
  expect_true(any(startsWith(r$warnings, "ERROR") &
                  grepl("S1 = V*WT/70 for the observed", r$warnings, fixed = TRUE)))
  expect_true(any(grepl("not resolvable: V*WT/70", r$unsupported, fixed = TRUE)))
})

test_that(".residual_scale() divides the volume out of a scale", {
  r <- .residual_scale
  expect_null(r("V", "V"))
  expect_identical(r("V/1000", "V"), 0.001)
  expect_identical(r("V * 1000", "V"), 1000)
  expect_equal(r("V * MW", "V"), "MW")
  expect_equal(r("MW * V", "V"), "MW")
  expect_equal(r("(V * MW)/2", "V"), "0.5 * MW")
  expect_equal(r(1000, "V"), "1000/V")
  expect_equal(r("MW", "V"), "MW/V")
  # The volume in a denominator is not cancelled by dividing again.
  expect_equal(r("MW/V", "V"), "MW/(V * V)")
})
