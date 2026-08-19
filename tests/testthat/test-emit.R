# -- .fmt_num -----------------------------------------------------------------

test_that(".fmt_num handles Inf and -Inf without producing invalid text", {
  expect_equal(.fmt_num(Inf),  "1e15")
  expect_equal(.fmt_num(-Inf), "-1e15")
  expect_equal(.fmt_num(0.134), "0.134")
  expect_equal(.fmt_num(10L),   "10.0")
})

# -- helpers ------------------------------------------------------------------

warfarin_1cpt_ir <- function(...) {
  new_ferx_ir(
    source_format = "nonmem",
    source_file   = "1cpt_oral.ctl",
    thetas = list(
      list(name = "TVCL", init = 0.134, lower = 0.001, upper = 10.0),
      list(name = "TVV",  init = 8.1,   lower = 0.1,   upper = 500.0),
      list(name = "TVKA", init = 1.0,   lower = 0.01,  upper = 50.0)
    ),
    omegas = list(
      list(type = "diagonal", names = "ETA_CL", values = 0.07),
      list(type = "diagonal", names = "ETA_V",  values = 0.02),
      list(type = "diagonal", names = "ETA_KA", values = 0.40)
    ),
    sigmas = list(
      list(name = "PROP_ERR", value = 0.01, scale = "sd")
    ),
    indiv_params = list(
      list(lhs = "CL", rhs = "TVCL * exp(ETA_CL)"),
      list(lhs = "V",  rhs = "TVV  * exp(ETA_V)"),
      list(lhs = "KA", rhs = "TVKA * exp(ETA_KA)")
    ),
    structural = list(
      type    = "pk_macro",
      pk_call = "one_cpt_oral",
      pk_args = list(cl = "CL", v = "V", ka = "KA")
    ),
    error_model = list(
      list(dv = "DV", type = "proportional", params = "PROP_ERR")
    ),
    fit_options = list(method = "foce", maxiter = 300L, covariance = TRUE),
    ...
  )
}

# -- snapshot tests -----------------------------------------------------------

test_that("1-cpt oral IR emits correct .ferx (snapshot)", {
  expect_snapshot(cat(emit_ferx(warfarin_1cpt_ir())))
})

test_that("2-cpt oral with covariates emits correct .ferx (snapshot)", {
  ir <- new_ferx_ir(
    source_format = "nonmem",
    source_file   = "2cpt_oral_cov.ctl",
    thetas = list(
      list(name = "TVCL",       init = 5.0,  lower = 0.1,  upper = 100.0),
      list(name = "TVV1",       init = 50.0, lower = 1.0,  upper = 500.0),
      list(name = "TVQ",        init = 10.0, lower = 0.1,  upper = 100.0),
      list(name = "TVV2",       init = 100.0,lower = 1.0,  upper = 500.0),
      list(name = "TVKA",       init = 1.2,  lower = 0.01, upper = 10.0),
      list(name = "THETA_WT",   init = 0.75, lower = 0.01, upper = 5.0),
      list(name = "THETA_CRCL", init = 0.50, lower = 0.01, upper = 5.0)
    ),
    omegas = list(
      list(type = "diagonal", names = "ETA_CL", values = 0.10),
      list(type = "diagonal", names = "ETA_V1", values = 0.10),
      list(type = "diagonal", names = "ETA_Q",  values = 0.05),
      list(type = "diagonal", names = "ETA_V2", values = 0.05),
      list(type = "diagonal", names = "ETA_KA", values = 0.15)
    ),
    sigmas = list(
      list(name = "PROP_ERR", value = 0.02, scale = "sd")
    ),
    indiv_params = list(
      list(lhs = "CL", rhs = "TVCL * (WT / 70)^THETA_WT * (CRCL / 100)^THETA_CRCL * exp(ETA_CL)"),
      list(lhs = "V1", rhs = "TVV1 * (WT / 70)^THETA_WT * exp(ETA_V1)"),
      list(lhs = "Q",  rhs = "TVQ  * exp(ETA_Q)"),
      list(lhs = "V2", rhs = "TVV2 * exp(ETA_V2)"),
      list(lhs = "KA", rhs = "TVKA * exp(ETA_KA)")
    ),
    structural = list(
      type    = "pk_macro",
      pk_call = "two_cpt_oral",
      pk_args = list(cl = "CL", v1 = "V1", q = "Q", v2 = "V2", ka = "KA")
    ),
    error_model = list(
      list(dv = "DV", type = "proportional", params = "PROP_ERR")
    ),
    fit_options = list(method = "focei", maxiter = 500L, covariance = TRUE)
  )
  expect_snapshot(cat(emit_ferx(ir)))
})

# -- section-level checks -----------------------------------------------------

test_that("block omega emits block_omega line", {
  ir  <- new_ferx_ir(
    omegas = list(
      list(type = "block", names = c("ETA_CL", "ETA_V"), values = c(0.07, 0.02, 0.02))
    )
  )
  out <- emit_ferx(ir)
  expect_match(out, "block_omega (ETA_CL, ETA_V) = [0.07, 0.02, 0.02]",
               fixed = TRUE)
})

test_that("diagonal omega emits omega line", {
  ir  <- new_ferx_ir(
    omegas = list(list(type = "diagonal", names = "ETA_CL", values = 0.07))
  )
  expect_match(emit_ferx(ir), "omega ETA_CL ~ 0.07", fixed = TRUE)
})

test_that("sigma on variance scale emits no (sd) suffix", {
  ir  <- new_ferx_ir(
    sigmas = list(list(name = "PROP_ERR", value = 0.0001, scale = "var"))
  )
  out <- emit_ferx(ir)
  expect_match(out, "sigma PROP_ERR ~ 0.0001", fixed = TRUE)
  expect_false(grepl("(sd)", out, fixed = TRUE))
})

test_that("sigma on sd scale appends (sd)", {
  ir  <- new_ferx_ir(
    sigmas = list(list(name = "PROP_ERR", value = 0.01, scale = "sd"))
  )
  expect_match(emit_ferx(ir), "sigma PROP_ERR ~ 0.01 (sd)", fixed = TRUE)
})

test_that("ODE structural model emits [odes] and ode() call", {
  ir <- new_ferx_ir(
    structural = list(
      type     = "ode",
      obs_cmt  = "central",
      states   = c("depot", "central")
    ),
    odes = list(
      list(state = "depot",   rhs = "-KA * depot"),
      list(state = "central", rhs = "KA * depot / V - (CL / V) * central")
    )
  )
  out <- emit_ferx(ir)
  expect_match(out, "[odes]",                              fixed = TRUE)
  expect_match(out, "ode(obs_cmt=central, states=[depot, central])", fixed = TRUE)
  expect_match(out, "d/dt(depot) = -KA * depot",           fixed = TRUE)
})

test_that("diffusion section emits between [odes] and [error_model]", {
  ir <- new_ferx_ir(
    structural  = list(type = "ode", obs_cmt = "central", states = c("depot", "central")),
    odes        = list(list(state = "depot",   rhs = "-KA * depot"),
                       list(state = "central", rhs = "KA * depot / V - CL/V * central")),
    diffusion   = list(list(state = "central", value = 0.01)),
    error_model = list(list(dv = "DV", type = "proportional", params = "PROP_ERR"))
  )
  out   <- emit_ferx(ir)
  pos_d <- regexpr("[diffusion]",   out, fixed = TRUE)
  pos_e <- regexpr("[error_model]", out, fixed = TRUE)
  expect_true(pos_d > 0 && pos_d < pos_e)
  expect_match(out, "central ~ 0.01", fixed = TRUE)
})

test_that("IOV model emits kappa line and iov_column in fit_options", {
  ir <- new_ferx_ir(
    kappas = list(list(name = "KAPPA_CL", value = 0.04)),
    fit_options = list(method = "foce", covariance = FALSE, iov_column = "OCC")
  )
  out <- emit_ferx(ir)
  expect_match(out, "kappa KAPPA_CL ~ 0.04",  fixed = TRUE)
  expect_match(out, "iov_column = OCC",        fixed = TRUE)
})

test_that("scaling section emits obs_scale", {
  ir  <- new_ferx_ir(scaling = list(obs_scale = 1000))
  out <- emit_ferx(ir)
  expect_match(out, "[scaling]",        fixed = TRUE)
  expect_match(out, "obs_scale = 1000", fixed = TRUE)
})

test_that("scaling section appears before [fit_options]", {
  ir <- new_ferx_ir(
    scaling     = list(obs_scale = 1000),
    fit_options = list(method = "foce")
  )
  out   <- emit_ferx(ir)
  pos_s <- regexpr("[scaling]",     out, fixed = TRUE)
  pos_f <- regexpr("[fit_options]", out, fixed = TRUE)
  expect_true(pos_s > 0 && pos_s < pos_f)
})

test_that("unsupported feature emits # WARNING: comment", {
  ir  <- new_ferx_ir(unsupported = c("MIXTURE model"))
  out <- emit_ferx(ir)
  expect_match(out, "# WARNING: MIXTURE model", fixed = TRUE)
})

test_that("combined error model emits both params", {
  ir <- new_ferx_ir(
    error_model = list(
      list(dv = "DV", type = "combined", params = c("ADD_ERR", "PROP_ERR"))
    )
  )
  expect_match(emit_ferx(ir), "DV ~ combined(ADD_ERR, PROP_ERR)", fixed = TRUE)
})

test_that("fit_options logical covariance formats as true/false", {
  ir <- new_ferx_ir(fit_options = list(method = "focei", covariance = TRUE))
  expect_match(emit_ferx(ir), "covariance = true", fixed = TRUE)

  ir2 <- new_ferx_ir(fit_options = list(method = "focei", covariance = FALSE))
  expect_match(emit_ferx(ir2), "covariance = false", fixed = TRUE)
})

test_that("method appears before iov_column in fit_options", {
  ir  <- new_ferx_ir(fit_options = list(iov_column = "OCC", method = "foce"))
  out <- emit_ferx(ir)
  pos_m <- regexpr("method",     out)
  pos_i <- regexpr("iov_column", out)
  expect_true(pos_m < pos_i)
})

# -- ordered statement lists (issue #6 phase 5a) ------------------------------

# Minimal ODE IR whose [odes] block is whatever `stmts` says. Deliberately not
# built on warfarin_1cpt_ir(): this tier is about statement RENDERING, and a
# fixture carrying an unrelated pk macro would obscure which block is under test.
stmt_ode_ir <- function(stmts, indiv = list(list(lhs = "KE", rhs = "TVKE")),
                        inits = list()) {
  new_ferx_ir(
    source_format = "nonmem",
    thetas = list(list(name = "TVKE", init = 0.1, lower = 0.001, upper = 10)),
    omegas = list(list(type = "diagonal", names = "ETA_KE", values = 0.09)),
    sigmas = list(list(name = "PROP", value = 0.15, scale = "sd")),
    indiv_params = indiv,
    structural = list(type = "ode", obs_cmt = "CENTRAL", states = "CENTRAL"),
    initial_conditions = inits,
    odes = stmts,
    error_model = list(list(dv = "DV", type = "proportional", params = "PROP")))
}

odes_block <- function(txt) {
  ln <- strsplit(txt, "\n")[[1]]
  i  <- grep("^\\[odes\\]", ln)
  rest <- ln[(i + 1):length(ln)]
  e <- grep("^\\[", rest)
  if (length(e)) rest <- rest[seq_len(e[1] - 1)]
  # Sections are joined with a blank line, so the last one trails an empty
  # string. Drop it here rather than writing it into every expectation.
  while (length(rest) && !nzchar(rest[length(rest)])) rest <- rest[-length(rest)]
  rest
}

test_that("an entry with no kind still renders as it did before statement lists", {
  # The compatibility rule the hand-built IRs in this file and test-ir.R rely on:
  # a bare list(state, rhs) in odes is a d/dt, a bare list(lhs, rhs) in
  # indiv_params is an assignment. Neither carries a `kind` field to misread.
  txt <- emit_ferx(stmt_ode_ir(list(list(state = "CENTRAL", rhs = "-KE * CENTRAL"))))
  expect_match(txt, "  d/dt(CENTRAL) = -KE * CENTRAL", fixed = TRUE)
  expect_match(txt, "  KE = TVKE", fixed = TRUE)
})

test_that("a single-statement if is braced and rendered inline", {
  # Braces are not stylistic. Measured against ferx 0.3.0, the unbraced form is
  # `E_PARSE: Expected `{` after if-condition` -- and NONMEM writes its commonest
  # conditional unbraced (`IF (DSC.LT.0.0) DSC = 0.0`), so this is the shape that
  # matters most.
  txt <- emit_ferx(stmt_ode_ir(list(
    list(kind = "if", cond = "CT < 0",
         then = list(list(kind = "assign", lhs = "CT", rhs = "0"))),
    list(state = "CENTRAL", rhs = "-KE * CENTRAL"))))
  expect_match(txt, "  if (CT < 0) { CT = 0 }", fixed = TRUE)
})

test_that("an if/else with one statement per arm renders inline", {
  txt <- emit_ferx(stmt_ode_ir(list(
    list(kind = "if", cond = "BB >= 0",
         then  = list(list(kind = "assign", lhs = "CF", rhs = "0.5 * (BB + DD)")),
         else_ = list(list(kind = "assign", lhs = "CF", rhs = "2 * KSS"))),
    list(state = "CENTRAL", rhs = "-KE * CENTRAL"))))
  expect_match(txt, "  if (BB >= 0) { CF = 0.5 * (BB + DD) } else { CF = 2 * KSS }",
               fixed = TRUE)
})

test_that("a multi-statement arm renders multi-line with indented body", {
  txt <- emit_ferx(stmt_ode_ir(list(
    list(kind = "if", cond = "CT > 100",
         then  = list(list(kind = "assign", lhs = "AA",  rhs = "CT * 2"),
                      list(kind = "assign", lhs = "SCL", rhs = "AA / (AA + 1)")),
         else_ = list(list(kind = "assign", lhs = "SCL", rhs = "0.5"))),
    list(state = "CENTRAL", rhs = "-KE * CENTRAL * SCL"))))
  expect_equal(odes_block(txt),
               c("  if (CT > 100) {",
                 "    AA = CT * 2",
                 "    SCL = AA / (AA + 1)",
                 "  } else {",
                 "    SCL = 0.5",
                 "  }",
                 "  d/dt(CENTRAL) = -KE * CENTRAL * SCL"))
})

test_that("a nested if is never collapsed onto one line", {
  # Inlining a nested `if` produces `if (a) { if (b) { x = 1 } else ... } ...`,
  # which parses but is unreadable at the exact place a reader most needs to
  # follow the branching.
  txt <- emit_ferx(stmt_ode_ir(list(
    list(kind = "if", cond = "CT >= 0",
         then = list(list(kind = "if", cond = "CT > 100",
                          then  = list(list(kind = "assign", lhs = "SCL", rhs = "0.5")),
                          else_ = list(list(kind = "assign", lhs = "SCL", rhs = "1.0"))))),
    list(state = "CENTRAL", rhs = "-KE * CENTRAL * SCL"))))
  expect_equal(odes_block(txt),
               c("  if (CT >= 0) {",
                 "    if (CT > 100) { SCL = 0.5 } else { SCL = 1.0 }",
                 "  }",
                 "  d/dt(CENTRAL) = -KE * CENTRAL * SCL"))
})

test_that("statement order is preserved exactly, init() first", {
  # The correctness property of the whole block. ferx has NO use-before-def check
  # in [odes]: an intermediate below the d/dt line that reads it stays valid,
  # reads a stale slot, and collapses PRED to a constant with no diagnostic. A
  # renderer that sorted or grouped would be silently wrong.
  txt <- emit_ferx(stmt_ode_ir(
    list(list(kind = "assign", lhs = "CT",  rhs = "CENTRAL / V"),
         list(kind = "assign", lhs = "SCL", rhs = "CT + 1"),
         list(state = "CENTRAL", rhs = "-KE * CENTRAL * SCL"),
         list(kind = "assign", lhs = "ZZ", rhs = "1")),
    inits = list(list(state = "CENTRAL", rhs = "0.0"))))
  expect_equal(odes_block(txt),
               c("  init(CENTRAL) = 0.0",
                 "  CT = CENTRAL / V",
                 "  SCL = CT + 1",
                 "  d/dt(CENTRAL) = -KE * CENTRAL * SCL",
                 "  ZZ = 1"))
})

test_that("an illegal name declared inside an if branch is still caught", {
  # The census aborts, so a walk that stopped at the top level would not report a
  # weaker result -- it would report nothing and let the name reach the engine.
  expect_error(
    emit_ferx(stmt_ode_ir(list(
      list(kind = "if", cond = "CT < 0",
           then = list(list(kind = "assign", lhs = "c.BAD", rhs = "0"))),
      list(state = "CENTRAL", rhs = "-KE * CENTRAL")))),
    "c.BAD")
})

test_that("an unknown statement kind is refused rather than skipped", {
  # Silently dropping it would emit a model missing a statement the IR carried,
  # which is the whole class of defect issue #6 is about.
  expect_error(
    emit_ferx(stmt_ode_ir(list(list(kind = "loop", lhs = "X", rhs = "1")))),
    "Unknown statement kind")
})
