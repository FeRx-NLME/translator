# -- mock helpers -------------------------------------------------------------

# Build a minimal mock rxUI list from explicit iniDf + lstExpr.
mock_ui <- function(iniDf, lstExpr = list()) {
  list(iniDf = iniDf, lstExpr = lstExpr)
}

# Build a minimal iniDf row list; rbind-friendly.
theta_row <- function(name, est, lower = 0, upper = Inf, fix = FALSE) {
  data.frame(ntheta = 1L, neta1 = NA_integer_, neta2 = NA_integer_,
             name = name, lower = lower, est = est, upper = upper,
             fix = fix, err = NA_character_, condition = NA_character_,
             stringsAsFactors = FALSE)
}
eta_row <- function(name, est, neta1, neta2 = neta1, condition = "id") {
  data.frame(ntheta = NA_integer_, neta1 = neta1, neta2 = neta2,
             name = name, lower = -Inf, est = est, upper = Inf,
             fix = FALSE, err = NA_character_, condition = condition,
             stringsAsFactors = FALSE)
}
sigma_row <- function(name, est, err = "prop") {
  data.frame(ntheta = NA_integer_, neta1 = NA_integer_, neta2 = NA_integer_,
             name = name, lower = 0, est = est, upper = Inf,
             fix = FALSE, err = err, condition = NA_character_,
             stringsAsFactors = FALSE)
}

# Construct a d/dt assignment as nonmem2rx/rxode2 would emit:
#   d/dt(STATE) <- rhs
# R parses d/dt(STATE) as `/`(d, dt(STATE)), wrapped in a `<-` call.
ddt <- function(state, rhs) {
  lhs <- as.call(list(as.name("/"), as.name("d"),
                      as.call(list(as.name("dt"), as.name(state)))))
  as.call(list(as.name("<-"), lhs, rhs))
}

# -- .extract_thetas ----------------------------------------------------------

test_that("extracts theta name, init, bounds", {
  ini <- theta_row("tvcl", est = 0.134, lower = 0.001, upper = 10)
  out <- .extract_thetas(ini)
  expect_length(out$thetas, 1L)
  expect_equal(out$thetas[[1]]$name,  "TVCL")
  expect_equal(out$thetas[[1]]$init,  0.134)
  expect_equal(out$thetas[[1]]$lower, 0.001)
  expect_equal(out$thetas[[1]]$upper, 10)
})

test_that("normalises theta name with dots", {
  ini <- theta_row("tv.cl", est = 0.1)
  out <- .extract_thetas(ini)
  expect_equal(out$thetas[[1]]$name, "TV_CL")
})

test_that("FIXED theta emits INFO warning", {
  ini <- theta_row("tvcl", est = 0.134, fix = TRUE)
  out <- .extract_thetas(ini)
  expect_length(out$warnings, 1L)
  expect_match(out$warnings[1], "INFO")
  expect_match(out$warnings[1], "TVCL")
  expect_match(out$warnings[1], "FIXED")
})

test_that("multiple thetas extracted in order", {
  ini <- rbind(theta_row("tvcl", 0.134), theta_row("tvv", 8.1))
  out <- .extract_thetas(ini)
  expect_equal(vapply(out$thetas, `[[`, "", "name"), c("TVCL", "TVV"))
})

# -- .extract_omegas ----------------------------------------------------------

test_that("extracts diagonal omega", {
  ini <- eta_row("eta.cl", 0.07, neta1 = 1L)
  out <- .extract_omegas(ini)
  expect_length(out$omegas, 1L)
  expect_equal(out$omegas[[1]]$type,   "diagonal")
  expect_equal(out$omegas[[1]]$names,  "ETA_CL")
  expect_equal(out$omegas[[1]]$values, 0.07)
})

test_that("extracts multiple diagonal omegas in order", {
  ini <- rbind(eta_row("eta.cl", 0.07, 1L), eta_row("eta.v", 0.02, 2L))
  out <- .extract_omegas(ini)
  expect_length(out$omegas, 2L)
  expect_equal(out$omegas[[1]]$names, "ETA_CL")
  expect_equal(out$omegas[[2]]$names, "ETA_V")
})

test_that("block omega detected from off-diagonal entry", {
  ini <- rbind(
    eta_row("eta.cl", 0.07, neta1 = 1L, neta2 = 1L),
    eta_row("eta.v",  0.02, neta1 = 2L, neta2 = 1L),
    eta_row("eta.v",  0.02, neta1 = 2L, neta2 = 2L)
  )
  out <- .extract_omegas(ini)
  expect_length(out$omegas, 1L)
  expect_equal(out$omegas[[1]]$type,   "block")
  expect_equal(out$omegas[[1]]$names,  c("ETA_CL", "ETA_V"))
  expect_equal(out$omegas[[1]]$values, c(0.07, 0.02, 0.02))
})

test_that("block + diagonal omegas: block first, diagonal after", {
  ini <- rbind(
    eta_row("eta.cl", 0.07, neta1 = 1L, neta2 = 1L),
    eta_row("eta.v",  0.02, neta1 = 2L, neta2 = 1L),
    eta_row("eta.v",  0.02, neta1 = 2L, neta2 = 2L),
    eta_row("eta.ka", 0.40, neta1 = 3L, neta2 = 3L)
  )
  out <- .extract_omegas(ini)
  expect_length(out$omegas, 2L)
  expect_equal(out$omegas[[1]]$type,  "block")
  expect_equal(out$omegas[[2]]$type,  "diagonal")
  expect_equal(out$omegas[[2]]$names, "ETA_KA")
})

test_that("empty iniDf gives empty omegas", {
  ini <- rbind(theta_row("tvcl", 0.1))
  out <- .extract_omegas(ini)
  expect_equal(out$omegas, list())
})

# -- .extract_kappas ----------------------------------------------------------

test_that("extracts IOV kappa from non-id condition", {
  ini <- eta_row("kappa.cl", 0.04, neta1 = 1L, condition = "OCC")
  out <- .extract_kappas(ini)
  expect_length(out$kappas, 1L)
  expect_equal(out$kappas[[1]]$name,  "KAPPA_CL")
  expect_equal(out$kappas[[1]]$value, 0.04)
  expect_equal(out$iov_column, "OCC")
})

test_that("block IOV emits WARN", {
  ini <- rbind(
    eta_row("kappa.cl", 0.04, neta1 = 1L, neta2 = 1L, condition = "OCC"),
    eta_row("kappa.v",  0.02, neta1 = 2L, neta2 = 1L, condition = "OCC"),
    eta_row("kappa.v",  0.02, neta1 = 2L, neta2 = 2L, condition = "OCC")
  )
  out <- .extract_kappas(ini)
  expect_match(out$warnings[1], "WARN")
})

test_that("empty IOV gives empty kappas", {
  ini <- eta_row("eta.cl", 0.07, neta1 = 1L, condition = "id")
  out <- .extract_kappas(ini)
  expect_equal(out$kappas, list())
})

# -- .extract_sigmas ----------------------------------------------------------

test_that("extracts sigma on sd scale", {
  ini <- sigma_row("err.prop", 0.01, err = "prop")
  out <- .extract_sigmas(ini)
  expect_length(out$sigmas, 1L)
  expect_equal(out$sigmas[[1]]$name,  "ERR_PROP")
  expect_equal(out$sigmas[[1]]$value, 0.01)
  expect_equal(out$sigmas[[1]]$scale, "sd")
})

# -- expression classifiers ---------------------------------------------------

test_that(".is_ddt_lhs detects d/dt LHS in assignment", {
  assign_expr <- ddt("depot", quote(-KA * depot))
  expect_true(.is_ddt_lhs(assign_expr[[2]]))
  expect_false(.is_ddt_lhs(quote(cl)))
  expect_false(.is_ddt_lhs(quote(linCmt())))
})

test_that(".is_tilde detects tilde expression", {
  expect_true(.is_tilde(quote(linCmt() ~ prop(err.prop))))
  expect_true(.is_tilde(quote(DV ~ prop(err.prop))))
  expect_false(.is_tilde(quote(cl <- tvcl)))
})

test_that(".is_lincmt_tilde detects linCmt on LHS", {
  expect_true(.is_lincmt_tilde(quote(linCmt() ~ prop(err.prop))))
  expect_false(.is_lincmt_tilde(quote(DV ~ prop(err.prop))))
})

test_that(".is_assignment detects <- and = assignment", {
  expect_true(.is_assignment(quote(cl <- tvcl)))
  expect_true(.is_assignment(as.call(list(as.name("="), as.name("cl"), as.name("tvcl")))))
  expect_false(.is_assignment(quote(linCmt() ~ prop(err))))
  # d/dt(depot) <- rhs IS an assignment (the d/dt is in the LHS)
  expect_true(.is_assignment(ddt("depot", quote(-KA))))
})

# -- .normalise_expr ----------------------------------------------------------

test_that("normalises eta name in expression", {
  map  <- c("eta.cl" = "ETA_CL", "tvcl" = "TVCL")
  expr <- quote(tvcl * exp(eta.cl))
  out  <- .normalise_expr(expr, map)
  expect_equal(deparse(out), "TVCL * exp(ETA_CL)")
})

test_that("leaves unknown names unchanged", {
  map  <- c("eta.cl" = "ETA_CL")
  expr <- quote(WT / 70)
  out  <- .normalise_expr(expr, map)
  expect_equal(deparse(out), "WT/70")
})

test_that("does not normalise function name", {
  map  <- c("exp" = "EXP")
  expr <- quote(exp(eta.cl))
  out  <- .normalise_expr(expr, map)
  expect_equal(deparse(out), "exp(eta.cl)")
})

# -- .parse_model_exprs -------------------------------------------------------

test_that("assignment parsed and normalised", {
  ini <- rbind(theta_row("tvcl", 0.134), eta_row("eta.cl", 0.07, 1L))
  map <- .norm_map_from_ini(ini)
  lst <- list(quote(cl <- tvcl * exp(eta.cl)))
  out <- .parse_model_exprs(lst, map)
  expect_length(out$indiv_params, 1L)
  expect_equal(out$indiv_params[[1]]$lhs, "CL")
  expect_match(out$indiv_params[[1]]$rhs, "TVCL")
  expect_match(out$indiv_params[[1]]$rhs, "ETA_CL")
})

test_that("linCmt tilde sets structural type lincmt", {
  map <- .norm_map_from_ini(sigma_row("err.prop", 0.01))
  lst <- list(quote(linCmt() ~ prop(err.prop)))
  out <- .parse_model_exprs(lst, map)
  expect_equal(out$structural$type, "lincmt")
})

test_that("d/dt assignment sets structural type ode", {
  map <- list()
  lst <- list(ddt("depot", quote(-KA * depot)))
  out <- .parse_model_exprs(lst, map)
  expect_equal(out$structural$type, "ode")
  expect_length(out$odes, 1L)
  expect_equal(out$odes[[1]]$state, "depot")
  expect_match(out$odes[[1]]$rhs,   "KA")
})

test_that("proportional error parsed", {
  map <- .norm_map_from_ini(sigma_row("err.prop", 0.01))
  out <- .parse_error_rhs(quote(prop(err.prop)), map)
  expect_equal(out$type,   "proportional")
  expect_equal(out$params, "ERR_PROP")
})

test_that("additive error parsed", {
  map <- .norm_map_from_ini(sigma_row("err.add", 0.5, err = "add"))
  out <- .parse_error_rhs(quote(add(err.add)), map)
  expect_equal(out$type,   "additive")
  expect_equal(out$params, "ERR_ADD")
})

test_that("combined error parsed", {
  map <- c("err.add" = "ERR_ADD", "err.prop" = "ERR_PROP")
  out <- .parse_error_rhs(quote(add(err.add) + prop(err.prop)), map)
  expect_equal(out$type,   "combined")
  expect_equal(out$params, c("ERR_ADD", "ERR_PROP"))
})

test_that("combined error parsed when prop comes first (prop + add)", {
  map <- c("err.add" = "ERR_ADD", "err.prop" = "ERR_PROP")
  out <- .parse_error_rhs(quote(prop(err.prop) + add(err.add)), map)
  expect_equal(out$type,   "combined")
  expect_equal(out$params, c("ERR_ADD", "ERR_PROP"))
})

# -- .infer_pk_macro ----------------------------------------------------------

test_that("1-cpt oral inferred from ka + v (no q)", {
  params <- list(list(lhs = "CL"), list(lhs = "V"), list(lhs = "KA"))
  out    <- .infer_pk_macro(params)
  expect_equal(out$pk_call, "one_cpt_oral")
  expect_equal(out$pk_args$cl, "CL")
  expect_equal(out$pk_args$v,  "V")
  expect_equal(out$pk_args$ka, "KA")
})

test_that("1-cpt iv inferred when no ka", {
  params <- list(list(lhs = "CL"), list(lhs = "V"))
  out    <- .infer_pk_macro(params)
  expect_equal(out$pk_call, "one_cpt_iv_bolus")
  expect_null(out$pk_args$ka)
})

test_that("2-cpt oral inferred from ka + q", {
  params <- list(list(lhs = "CL"), list(lhs = "V1"), list(lhs = "Q"),
                 list(lhs = "V2"), list(lhs = "KA"))
  out    <- .infer_pk_macro(params)
  expect_equal(out$pk_call, "two_cpt_oral")
  expect_equal(out$pk_args$v1, "V1")
  expect_equal(out$pk_args$q,  "Q")
  expect_equal(out$pk_args$ka, "KA")
})

test_that("v alias: V used when v1 expected", {
  params <- list(list(lhs = "CL"), list(lhs = "V"), list(lhs = "Q"),
                 list(lhs = "V2"), list(lhs = "KA"))
  out    <- .infer_pk_macro(params)
  expect_equal(out$pk_call, "two_cpt_oral")
  expect_equal(out$pk_args$v1, "V")
})

test_that("3-cpt oral emits ERROR and NA pk_call", {
  params <- list(list(lhs = "CL"), list(lhs = "V1"), list(lhs = "Q"),
                 list(lhs = "V2"), list(lhs = "Q2"), list(lhs = "KA"))
  out    <- .infer_pk_macro(params)
  expect_true(is.na(out$pk_call))
  expect_match(out$warnings[1], "ERROR")
  expect_length(out$unsupported, 1L)
})

test_that("3-cpt IV bolus (Q2, no KA) also unsupported -- ferx has no 3-cpt analytical", {
  params <- list(list(lhs = "CL"), list(lhs = "V1"), list(lhs = "Q2"),
                 list(lhs = "V2"), list(lhs = "V3"))
  out    <- .infer_pk_macro(params)
  expect_true(is.na(out$pk_call))
  expect_match(out$warnings[1], "ERROR")
  expect_length(out$unsupported, 1L)
  expect_match(out$unsupported[1], "three_cpt_iv_bolus")
})

test_that("bioavailability f added to pk_args when present", {
  params <- list(list(lhs = "CL"), list(lhs = "V"), list(lhs = "KA"),
                 list(lhs = "F"))
  out    <- .infer_pk_macro(params)
  expect_equal(out$pk_args$f, "F")
})

# -- rxui_to_ir integration (mock UI) ----------------------------------------

test_that("rxui_to_ir produces ferx_ir from mock 1-cpt oral", {
  ini <- rbind(
    theta_row("tvcl", 0.134, 0.001, 10),
    theta_row("tvv",  8.1,   0.1,   500),
    theta_row("tvka", 1.0,   0.01,  50),
    eta_row("eta.cl", 0.07, 1L),
    eta_row("eta.v",  0.02, 2L),
    sigma_row("err.prop", 0.01)
  )
  lst <- list(
    quote(cl <- tvcl * exp(eta.cl)),
    quote(v  <- tvv  * exp(eta.v)),
    quote(ka <- tvka),
    quote(linCmt() ~ prop(err.prop))
  )
  ir <- rxui_to_ir(mock_ui(ini, lst), source_format = "nlmixr2")

  expect_s3_class(ir, "ferx_ir")
  expect_equal(ir$source_format, "nlmixr2")
  expect_length(ir$thetas, 3L)
  expect_length(ir$omegas, 2L)
  expect_length(ir$sigmas, 1L)
  expect_equal(ir$structural$type,    "pk_macro")
  expect_equal(ir$structural$pk_call, "one_cpt_oral")
  expect_equal(ir$indiv_params[[1]]$lhs, "CL")
  expect_match(ir$indiv_params[[1]]$rhs, "ETA_CL")
  expect_equal(ir$error_model[[1]]$type, "proportional")
  expect_true(isTRUE(ir$fit_options$covariance))
})

test_that("rxui_to_ir ODE model sets structural type ode with states and obs_cmt", {
  ini <- rbind(
    theta_row("tvcl", 0.134, 0.001, 10),
    theta_row("tvv",  8.1,   0.1,   500),
    theta_row("tvka", 1.0,   0.01,  50),
    eta_row("eta.cl", 0.07, 1L),
    sigma_row("err.prop", 0.01)
  )
  lst <- list(
    quote(cl <- tvcl * exp(eta.cl)),
    quote(v  <- tvv),
    quote(ka <- tvka),
    ddt("depot",   quote(-ka * depot)),
    ddt("central", quote(ka * depot / v - cl / v * central)),
    quote(DV ~ prop(err.prop))
  )
  ir <- rxui_to_ir(mock_ui(ini, lst))
  expect_equal(ir$structural$type,   "ode")
  expect_equal(ir$structural$states, c("depot", "central"))
  expect_equal(ir$structural$obs_cmt, "central")
  expect_length(ir$odes, 2L)
  expect_equal(ir$odes[[1]]$state, "depot")
})

test_that("rxui_to_ir 3-cpt oral: structural is empty (structural_model omitted)", {
  ini <- rbind(
    theta_row("tvcl", 0.1), theta_row("tvv1", 5), theta_row("tvq",  0.5),
    theta_row("tvv2", 10),  theta_row("tvq2", 0.2), theta_row("tvv3", 20),
    theta_row("tvka", 1.0),
    sigma_row("err.prop", 0.01)
  )
  lst <- list(
    quote(cl <- tvcl), quote(v1 <- tvv1), quote(q  <- tvq),
    quote(v2 <- tvv2), quote(q2 <- tvq2), quote(v3 <- tvv3),
    quote(ka <- tvka),
    quote(linCmt() ~ prop(err.prop))
  )
  ir <- rxui_to_ir(mock_ui(ini, lst))
  expect_length(ir$structural, 0L)
  expect_length(ir$unsupported, 1L)
  expect_match(ir$warnings[1], "ERROR")
})

test_that("rxui_to_ir IOV model sets iov_column in fit_options", {
  ini <- rbind(
    theta_row("tvcl", 0.134, 0.001, 10),
    eta_row("eta.cl",   0.07, 1L),
    eta_row("kappa.cl", 0.04, 1L, condition = "OCC"),
    sigma_row("err.prop", 0.01)
  )
  lst <- list(
    quote(cl <- tvcl * exp(eta.cl + kappa.cl)),
    quote(linCmt() ~ prop(err.prop))
  )
  ir <- rxui_to_ir(mock_ui(ini, lst))
  expect_length(ir$kappas, 1L)
  expect_equal(ir$fit_options$iov_column, "OCC")
})

# -- rxode2 integration tests (require rxode2) --------------------------------

test_that("1-cpt oral nlmixr2 function converts correctly", {
  skip_if_not_installed("rxode2")
  f_1cpt <- function() {
    ini({
      tvcl <- 0.134; tvv <- 8.1; tvka <- 1.0
      eta.cl ~ 0.07; eta.v ~ 0.02
      prop.err <- 0.01
    })
    model({
      cl <- tvcl * exp(eta.cl)
      v  <- tvv  * exp(eta.v)
      ka <- tvka
      linCmt() ~ prop(prop.err)
    })
  }
  ui <- rxode2::rxode2(f_1cpt)
  ir <- rxui_to_ir(ui, source_format = "nlmixr2")
  expect_equal(ir$structural$pk_call, "one_cpt_oral")
  expect_length(ir$thetas, 3L)
})

test_that("2-cpt oral nlmixr2 function with q infers two_cpt_oral", {
  skip_if_not_installed("rxode2")
  f_2cpt <- function() {
    ini({
      tvcl <- 5; tvv1 <- 50; tvq <- 10; tvv2 <- 100; tvka <- 1.2
      eta.cl ~ 0.10
      prop.err <- 0.02
    })
    model({
      cl <- tvcl * exp(eta.cl); v1 <- tvv1; q <- tvq; v2 <- tvv2; ka <- tvka
      linCmt() ~ prop(prop.err)
    })
  }
  ui <- rxode2::rxode2(f_2cpt)
  ir <- rxui_to_ir(ui, source_format = "nlmixr2")
  expect_equal(ir$structural$pk_call, "two_cpt_oral")
})

test_that("ODE nlmixr2 model sets structural type ode", {
  skip_if_not_installed("rxode2")
  f_ode <- function() {
    ini({ tvcl <- 0.134; tvv <- 8.1; tvka <- 1.0
          eta.cl ~ 0.07; prop.err <- 0.01 })
    model({
      cl <- tvcl * exp(eta.cl); v <- tvv; ka <- tvka
      d/dt(depot)   = -ka * depot
      d/dt(central) =  ka * depot / v - cl / v * central
      central ~ prop(prop.err)
    })
  }
  ui <- rxode2::rxode2(f_ode)
  ir <- rxui_to_ir(ui, source_format = "nlmixr2")
  expect_equal(ir$structural$type, "ode")
  expect_length(ir$odes, 2L)
})
