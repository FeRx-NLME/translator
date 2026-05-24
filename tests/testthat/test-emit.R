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
