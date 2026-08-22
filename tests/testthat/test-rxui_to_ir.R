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

test_that("FIXED theta sets fixed = TRUE in theta list", {
  ini <- theta_row("tvcl", est = 0.134, fix = TRUE)
  out <- .extract_thetas(ini)
  expect_true(out$thetas[[1]]$fixed)
  expect_length(out$warnings, 0L)
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

test_that(".iov_flattening_warnings flags a KAPPA-named IIV omega", {
  omegas <- list(
    list(type = "diagonal", names = "ETA_CL",   values = 0.10),
    list(type = "diagonal", names = "KAPPA_CL", values = 0.04)
  )
  w <- .iov_flattening_warnings(omegas)
  expect_length(w, 1L)
  expect_match(w, "KAPPA_CL", fixed = TRUE)
  expect_match(w, "inter-occasion")
  expect_match(w, "iov_column", fixed = TRUE)
})

test_that(".iov_flattening_warnings matches IOV* names and block etas", {
  omegas <- list(list(type = "block", names = c("ETA_CL", "IOV_V"),
                      values = c(0.1, 0.02, 0.05)))
  w <- .iov_flattening_warnings(omegas)
  expect_length(w, 1L)
  expect_match(w, "IOV_V", fixed = TRUE)
})

test_that(".iov_flattening_warnings is silent without IOV-named etas", {
  omegas <- list(list(type = "diagonal", names = "ETA_CL", values = 0.10))
  expect_length(.iov_flattening_warnings(omegas), 0L)
  expect_length(.iov_flattening_warnings(list()), 0L)
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

test_that("combined error parsed -- proportional first (ferx-core order)", {
  map <- c("err.add" = "ERR_ADD", "err.prop" = "ERR_PROP")
  out <- .parse_error_rhs(quote(add(err.add) + prop(err.prop)), map)
  expect_equal(out$type,   "combined")
  expect_equal(out$params, c("ERR_PROP", "ERR_ADD"))
})

test_that("combined error parsed when prop comes first (prop + add) -- proportional first", {
  map <- c("err.add" = "ERR_ADD", "err.prop" = "ERR_PROP")
  out <- .parse_error_rhs(quote(prop(err.prop) + add(err.add)), map)
  expect_equal(out$type,   "combined")
  expect_equal(out$params, c("ERR_PROP", "ERR_ADD"))
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
  expect_equal(out$pk_call, "one_cpt_iv")
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

test_that("3-cpt oral maps to three_cpt_oral pk_call", {
  # NONMEM ADVAN11 names the two inter-compartmental clearances Q2 (first,
  # paired with V2) and Q3 (second, paired with V3). Both must reach the pk
  # macro: q2 -> ferx slot Q, q3 -> ferx slot Q3.
  params <- list(list(lhs = "CL"), list(lhs = "V1"), list(lhs = "Q2"),
                 list(lhs = "V2"), list(lhs = "Q3"), list(lhs = "V3"),
                 list(lhs = "KA"))
  out    <- .infer_pk_macro(params)
  expect_equal(out$pk_call, "three_cpt_oral")
  expect_length(out$unsupported, 0L)
  expect_equal(out$pk_args$ka, "KA")
  expect_equal(out$pk_args$q2, "Q2")
  expect_equal(out$pk_args$q3, "Q3")
})

test_that("3-cpt IV maps to three_cpt_iv pk_call", {
  params <- list(list(lhs = "CL"), list(lhs = "V1"), list(lhs = "Q2"),
                 list(lhs = "V2"), list(lhs = "Q3"), list(lhs = "V3"))
  out    <- .infer_pk_macro(params)
  expect_equal(out$pk_call, "three_cpt_iv")
  expect_length(out$unsupported, 0L)
  expect_equal(out$pk_args$q2, "Q2")
  expect_equal(out$pk_args$q3, "Q3")
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
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst), source_format = "nlmixr2"))

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
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
  expect_equal(ir$structural$type,   "ode")
  expect_equal(ir$structural$states, c("depot", "central"))
  expect_equal(ir$structural$obs_cmt, "central")
  expect_length(ir$odes, 2L)
  expect_equal(ir$odes[[1]]$state, "depot")
})

test_that("rxui_to_ir 3-cpt oral: translates to three_cpt_oral pk macro", {
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
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
  expect_equal(ir$structural$type,    "pk_macro")
  expect_equal(ir$structural$pk_call, "three_cpt_oral")
  expect_length(ir$unsupported, 0L)
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
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
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
  ir <- suppressWarnings(rxui_to_ir(ui, source_format = "nlmixr2"))
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
  ir <- suppressWarnings(rxui_to_ir(ui, source_format = "nlmixr2"))
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
  ir <- suppressWarnings(rxui_to_ir(ui, source_format = "nlmixr2"))
  expect_equal(ir$structural$type, "ode")
  expect_length(ir$odes, 2L)
})

# -- random-effect name uniqueness --------------------------------------------

test_that(".uniquify_random_names leaves already-distinct names alone", {
  e <- list(list(name = "ETA_CL", raw = "eta.cl"),
            list(name = "ETA_V",  raw = "eta.v"))
  out <- .uniquify_random_names(e)
  expect_equal(vapply(out$entries, function(x) x$name, ""), c("ETA_CL", "ETA_V"))
  expect_length(out$warnings, 0L)
})

test_that("two source names normalising onto one spelling are made distinct", {
  # .norm() folds every illegal character onto `_`, which is many-to-one. The
  # theta channel has a uniqueness check and the state channel has one; this
  # channel had none, so both emitted `omega CL_IIV` and every reference
  # resolved to the first -- one IIV dropped, one double-counted, engine ok.
  e <- list(list(name = "CL_IIV", raw = "CL.IIV"),
            list(name = "CL_IIV", raw = "CL_IIV"))
  out <- .uniquify_random_names(e)
  expect_equal(vapply(out$entries, function(x) x$name, ""),
               c("CL_IIV", "CL_IIV_1"))
  expect_match(out$warnings, "CL_IIV_1", all = FALSE)
})

test_that("uniqueness is case-insensitive, as ferx compares names", {
  e <- list(list(name = "ETA_CL", raw = "eta.cl"),
            list(name = "eta_cl", raw = "ETA.CL"))
  out <- .uniquify_random_names(e)
  expect_equal(vapply(out$entries, function(x) x$name, ""),
               c("ETA_CL", "eta_cl_1"))
})

test_that("the first occurrence keeps its name", {
  # Renaming the first instead would churn a name the user already reads for the
  # benefit of a later duplicate.
  e <- list(list(name = "X", raw = "x1"), list(name = "X", raw = "x2"),
            list(name = "X", raw = "x3"))
  out <- .uniquify_random_names(e)
  expect_equal(vapply(out$entries, function(x) x$name, ""), c("X", "X_1", "X_2"))
})

test_that("a block omega uniquifies each of its names", {
  e <- list(list(name = c("A", "B"), raw = c("a", "b")),
            list(name = "A",         raw = "a2"))
  out <- .uniquify_random_names(e)
  expect_equal(out$entries[[1]]$name, c("A", "B"))
  expect_equal(out$entries[[2]]$name, "A_1")
})

test_that("colliding etas emit distinct omegas AND distinct references", {
  # The end-to-end version: the rename is worthless if name_map still points the
  # second eta's references at the name its twin took, which just moves the
  # merge from the declaration to the reference.
  ini <- rbind(theta_row("t.TCL", 1), theta_row("t.TV", 10),
               eta_row("CL.IIV", 0.09, 1L), eta_row("CL_IIV", 0.04, 2L))
  lst <- list(quote(cl <- t.TCL * exp(CL.IIV)),
              quote(v  <- t.TV  * exp(CL_IIV)),
              quote(d/dt(central) <- -cl/v * central))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  nms <- unlist(lapply(ir$omegas, function(o) o$names))
  expect_length(unique(nms), 2L)
  rhs <- vapply(ir$indiv_params, function(p) p$rhs, "")
  expect_match(rhs, "exp(CL_IIV)",   fixed = TRUE, all = FALSE)
  expect_match(rhs, "exp(CL_IIV_1)", fixed = TRUE, all = FALSE)
  # Each omega is referenced exactly once -- neither dropped nor double-counted.
  for (n in nms)
    expect_length(grep(paste0("\\b", n, "\\b"), rhs), 1L)
})

# -- theta / individual-parameter de-shadowing --------------------------------

test_that(".free_theta_name never returns the theta's own name", {
  # `base` was a candidate in .free_name's preference chain, so with both
  # prefixes taken the de-shadower could hand back the very name it was asked to
  # replace -- a non-NA answer, so the caller reports "theta 'CL' renamed to
  # 'CL'" and the shadowing survives silently. Reachable only if `taken` ever
  # stops containing the theta's own name; .deshadow_theta_names() seeds it, but
  # nothing stated that, so the guard belongs in the function.
  expect_equal(.free_theta_name("CL", c("TVCL", "THETA_CL")), "CL_1")
  for (taken in list(character(), "CL", c("TVCL", "THETA_CL"),
                     c("CL", "TVCL", "THETA_CL"),
                     c("TVCL", "THETA_CL", "CL_1", "CL_2")))
    expect_false(identical(.free_theta_name("CL", taken), "CL"))
})

test_that(".free_name still lets a non-theta caller keep its own name", {
  # The state sanitiser and the carrier-name path both WANT base as a fallback:
  # renaming `central` to CENTRAL_1 for nothing churns a name the user indexes
  # by. Only the theta path opts out.
  expect_equal(.free_name("central", c("DEPOT")), "central")
  expect_equal(.free_name("central", c("CENTRAL")), "central_1")
  expect_equal(.free_name("CL", c("X"), prefer = c("TVCL")), "TVCL")
  expect_equal(.free_name("CL", c("TVCL"), prefer = c("TVCL")), "CL")
  expect_equal(.free_name("CL", c("TVCL"), prefer = c("TVCL"),
                          allow_base = FALSE), "CL_1")
})

# -- internal invariant: states and individual parameters stay disjoint --------

test_that(".assert_state_param_disjoint passes on disjoint names", {
  expect_true(.assert_state_param_disjoint(
    list(list(state = "central"), list(state = "depot")),
    list(list(lhs = "CL"), list(lhs = "V"))))
  expect_true(.assert_state_param_disjoint(list(), list()))
})

test_that(".assert_state_param_disjoint stops, case-insensitively, on a clash", {
  # Unreachable through the current pipeline by construction -- .parse_model_exprs()
  # puts every state into `aux_vars` and pass 3 keeps aux_vars out of
  # indiv_params -- so it is called directly. It is a stop() and not a warning
  # because reaching it means the translator broke its own contract: the
  # assignment would be absorbed, dropped from [individual_parameters] and
  # inlined into itself, leaving an ODE that references a name nothing declares.
  expect_error(
    .assert_state_param_disjoint(list(list(state = "CENT")),
                                 list(list(lhs = "CENT"))),
    "internal error.*CENT.*names both an ODE state and an individual parameter")
  expect_error(
    .assert_state_param_disjoint(list(list(state = "cent")),
                                 list(list(lhs = "CENT"))),
    "internal error")
})

test_that("no rename when theta and individual parameter names are disjoint", {
  out <- .deshadow_theta_names(c("TVCL", "TVV"), c("CL", "V"))
  expect_true(all(is.na(out$map)))
  expect_length(out$warnings, 0L)
})

test_that("a theta colliding with an individual parameter is renamed to TV<name>", {
  out <- .deshadow_theta_names(c("CL", "TVV"), c("CL", "V"))
  expect_equal(out$map, c("TVCL", NA_character_))
  expect_equal(out$reasons[[1]], "shadow")
  expect_length(out$warnings, 0L)   # prose is written by the caller
})

test_that("collision detection is case-insensitive, as in ferx", {
  out <- .deshadow_theta_names("cl", "CL")
  expect_equal(out$map, "TVcl")
})

test_that("rename falls back when TV<name> is already taken", {
  out <- .deshadow_theta_names(c("CL", "TVCL"), "CL")
  expect_equal(out$map[1], "THETA_CL")
  expect_true(is.na(out$map[2]))
})

test_that("rename falls back to a numeric suffix when both prefixes are taken", {
  out <- .deshadow_theta_names(c("CL", "TVCL", "THETA_CL"), "CL")
  expect_equal(out$map[1], "CL_1")
})

test_that("rename avoids omega, kappa and sigma names", {
  out <- .deshadow_theta_names("CL", "CL", reserved = c("TVCL", "THETA_CL"))
  expect_equal(out$map, "CL_1")
})

test_that("two colliding thetas do not rename onto each other", {
  out <- .deshadow_theta_names(c("CL", "TVCL"), c("CL", "TVCL"))
  # CL cannot take TVCL (already a theta), so it falls through to THETA_CL;
  # TVCL takes TVTVCL. Ugly, but the names stay distinct, which is the point.
  expect_equal(sort(out$map), c("THETA_CL", "TVTVCL"))
  expect_length(unique(out$map), 2L)
})

test_that("duplicate theta names are made unique even without any shadowing", {
  # ferx accepts duplicate theta names and silently resolves every reference to
  # the first, leaving the second dead -- so this must be fixed, not just warned
  # about. Nothing here shadows an individual parameter.
  out <- .deshadow_theta_names(c("CL", "CL"), "NOTHING")
  expect_true(is.na(out$map[1]))
  expect_equal(out$map[2], "TVCL")
  expect_equal(out$reasons[[2]], "duplicate")
})

test_that("a duplicated AND shadowing theta gets two distinct names", {
  out <- .deshadow_theta_names(c("CL", "CL"), "CL")
  expect_length(unique(out$map), 2L)
  expect_false(anyNA(out$map))
})

test_that("rxui_to_ir de-shadows so downstream references reach the individual parameter", {
  ini <- rbind(theta_row("t.CL", 1), theta_row("t.V", 10),
               eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.CL * exp(eta1)),
              quote(v <- t.V),
              quote(k20 <- cl/v),
              ddt("central", quote(-k20 * central)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  # Only CL is renamed: `v <- t.V` is a bare theta alias, so V never becomes an
  # individual parameter and there is nothing for the theta V to shadow.
  expect_equal(vapply(ir$thetas, function(t) t$name, ""), c("TVCL", "V"))
  rhs <- setNames(vapply(ir$indiv_params, function(p) p$rhs, ""),
                  vapply(ir$indiv_params, function(p) p$lhs, ""))
  expect_equal(rhs[["CL"]], "TVCL * exp(ETA1)")
  # The reference that means the individual CL must stay bare, or ferx would
  # resolve it to the theta and the IIV would be silently dropped.
  expect_equal(rhs[["K20"]], "CL/V")
})

test_that("emitted theta names never collide with individual parameter names", {
  ini <- rbind(theta_row("t.CL", 1), theta_row("t.V", 10),
               eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.CL * exp(eta1)), quote(v <- t.V),
              ddt("central", quote(-cl/v * central)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
  expect_length(
    intersect(toupper(vapply(ir$thetas, function(t) t$name, "")),
              toupper(vapply(ir$indiv_params, function(p) p$lhs, ""))),
    0L
  )
})

test_that("every shadowed theta is renamed when all PK params carry an ETA", {
  ini <- rbind(theta_row("t.CL", 1), theta_row("t.V", 10),
               eta_row("eta1", 0.09, 1L), eta_row("eta2", 0.02, 2L))
  lst <- list(quote(cl <- t.CL * exp(eta1)),
              quote(v <- t.V * exp(eta2)),
              quote(k20 <- cl/v),
              ddt("central", quote(-k20 * central)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
  expect_equal(vapply(ir$thetas, function(t) t$name, ""), c("TVCL", "TVV"))
  rhs <- setNames(vapply(ir$indiv_params, function(p) p$rhs, ""),
                  vapply(ir$indiv_params, function(p) p$lhs, ""))
  expect_equal(rhs[["V"]], "TVV * exp(ETA2)")
  expect_equal(rhs[["K20"]], "CL/V")
})

test_that("a self-referential assignment resolves its RHS to the theta", {
  skip_if_not_installed("rxode2")
  # rxode2 lets a theta and the variable it defines share a name. The RHS `cl`
  # means the theta; installing the LHS alias before normalising the RHS turned
  # it into a self-reference that ferx rejects as a forward reference.
  f <- function() {
    ini({ cl <- 1.0; v <- 10.0; eta.cl ~ 0.09; prop.err <- 0.1 })
    model({ cl <- cl * exp(eta.cl); v <- v; linCmt() ~ prop(prop.err) })
  }
  ir <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))
  rhs <- setNames(vapply(ir$indiv_params, function(p) p$rhs, ""),
                  vapply(ir$indiv_params, function(p) p$lhs, ""))
  expect_equal(rhs[["CL"]], "TVCL * exp(ETA_CL)")
  expect_false(any(vapply(ir$thetas, function(t) t$name, "") == "CL"))
})

test_that("duplicate $THETA labels do not leave a reference dangling", {
  ini <- rbind(theta_row("theta1", 1), theta_row("theta2", 10),
               eta_row("eta1", 0.09, 1L))
  ini$label <- c("CL", "CL", NA_character_)
  lst <- list(quote(cl <- theta1 * exp(eta1)), quote(v <- theta2),
              ddt("central", quote(-cl/v * central)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  nms <- vapply(ir$thetas, function(t) t$name, "")
  expect_length(unique(nms), 2L)
  # Every symbol on an individual-parameter RHS must be a declared name, not a
  # stale reference that ferx would silently read as a covariate.
  rhs_syms <- unlist(lapply(ir$indiv_params, function(p)
    regmatches(p$rhs, gregexpr("[A-Za-z_][A-Za-z0-9_]*", p$rhs))[[1]]))
  known <- c(nms, vapply(ir$indiv_params, function(p) p$lhs, ""), "ETA1", "exp")
  expect_equal(setdiff(rhs_syms, known), character())
})

test_that("de-shadowing iterates: a rename must not create a new shadow", {
  # Pass 3 filters a theta alias by comparing names (`V <- V` is dropped as a
  # self-assignment), so renaming one theta can turn a filtered alias into an
  # individual parameter that did not exist in the previous parse -- and the
  # surviving theta then shadows it. Duplicate $THETA labels are how that
  # happens: only one of the pair is renamed. A single pass left
  # `theta V` beside `V = TVV`, with `CL = V * exp(ETA1)` reading the theta,
  # and the engine reported ok = TRUE with no diagnostic.
  ini <- rbind(theta_row("theta1", 1), theta_row("theta2", 10),
               eta_row("eta1", 0.09, 1L))
  ini$label <- c("V", "V", NA_character_)
  lst <- list(quote(v <- theta2), quote(cl <- theta1 * exp(eta1)),
              ddt("central", quote(-cl/v * central)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  th <- toupper(vapply(ir$thetas, function(t) t$name, ""))
  ip <- toupper(vapply(ir$indiv_params, function(p) p$lhs, ""))
  expect_equal(intersect(th, ip), character())
  expect_length(unique(th), length(th))
  expect_length(ir$unsupported, 0L)
})

test_that("an unresolvable name collision is reported, never emitted silently", {
  # The fixpoint is bounded, so the invariant is asserted rather than assumed.
  # Drive the real code path: with de-shadowing suppressed, pk_1cmt_oral.mod
  # emits `theta CL` beside `CL = ...`, and the assertion must catch it.
  # (The previous version of this test re-computed an intersection on a
  # hand-built IR and asserted nothing about the package at all.)
  skip_if_not_installed("nonmem2rx")
  local_mocked_bindings(
    .deshadow_theta_names = function(theta_names, indiv_names, reserved = character())
      list(map = rep(NA_character_, length(theta_names)), warnings = character())
  )
  ir <- suppressWarnings(rxui_to_ir(nonmem2rx::nonmem2rx(nm_path("pk_1cmt_oral.mod")),
                                    source_format = "nonmem"))
  expect_match(ir$unsupported, "theta/individual-parameter name collision", all = FALSE)
  expect_match(ir$warnings, "^ERROR \\| could not give theta", all = FALSE)
  # ...and it must reach the artefact, not only the result object.
  expect_match(emit_ferx(ir), "# WARNING: could not give theta")
})

test_that("covariate detection does not report iniDf names as covariates", {
  # name_map keys are raw iniDf names (`t.CL`), and .norm() does not strip the
  # prefix, so every theta reference was being reported as a covariate.
  ini <- rbind(theta_row("t.TVCL", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.TVCL * exp(eta1)),
              ddt("central", quote(-cl * central)))
  covs <- .covariate_names(lst, .norm_map_from_ini(ini))
  expect_false(any(grepl("^T_|^E_", covs)))
})

test_that("an ODE intermediate keeps the binding it had when written", {
  # Kills the mutation that reverts pass 2b to re-normalising the raw expression
  # with the FINAL name_map. De-shadowing makes that map time-varying, so `frac`
  # -- written when `cl` still meant the theta -- resolved to the individual
  # parameter. Both forms parse; ferx cannot tell them apart.
  #
  # Phase 5b emits the intermediate instead of substituting it into the d/dt
  # line, which makes the binding directly readable rather than buried in the
  # ODE right-hand side -- but the property under test is unchanged.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ cl <- 1.0; v <- 10.0; ka <- 1; prop.err <- 0.1; eta.cl ~ 0.09 })
    model({ frac <- central/cl; cl <- cl*exp(eta.cl); v <- v; ka <- ka
            d/dt(depot)   = -ka*depot
            d/dt(central) =  ka*depot - cl/v*central - frac
            central ~ prop(prop.err) })
  }
  ir  <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))
  frac <- Filter(function(o) identical(o$lhs, "FRAC"), ir$odes)
  expect_length(frac, 1L)
  rhs <- frac[[1]]$rhs

  # `frac` must carry the THETA's value, and it does so through a carrier
  # parameter rather than by naming the theta: a theta is not in scope in [odes]
  # (issue #6 defect 2), so the emitted reference used to be unresolvable. What
  # this test exists to pin is unchanged -- the term must not read the individual
  # parameter CL, whose value has the IIV applied. Asserted on the whole string
  # because "central/CL" is a substring of the correct "central/TVCL_ODE".
  #
  # The carrier is named off the THETA (`TVCL_ODE`), not the source name: `CL_ODE`
  # would read as "the CL used in the ODE", which is the individual value, and
  # telling those two apart is the whole point.
  expect_equal(rhs, "central/TVCL_ODE")
  # And the d/dt line reads the intermediate rather than re-deriving it.
  ddt_c <- Filter(function(o) identical(o$state, "central"), ir$odes)
  expect_match(ddt_c[[1]]$rhs, "FRAC", fixed = TRUE)
  ip <- vapply(ir$indiv_params, function(p) paste0(p$lhs, "=", p$rhs), "")
  expect_true("TVCL_ODE=TVCL" %in% ip)
  # No numbered fallback was needed, so no warning about one.
  expect_false(any(grepl("had to be numbered", ir$warnings)))
})

test_that("a carrier falls back to a numbered name loudly, never silently", {
  # Both preferred spellings taken: `CL` is a real individual parameter and
  # `TVCL_ODE` is a model variable in its own right. The numbered name is
  # positional -- it moves if another carrier is added before it -- so anything
  # indexing the emitted parameters by name breaks silently unless this is said
  # out loud.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ cl <- 1.0; v <- 10.0; ka <- 1; prop.err <- 0.1; eta.cl ~ 0.09 })
    model({ frac <- central/cl
            TVCL_ODE <- v * 2
            cl <- cl*exp(eta.cl); v <- v; ka <- ka
            d/dt(depot)   = -ka*depot
            d/dt(central) =  ka*depot - cl/v*central - frac + 0*TVCL_ODE
            central ~ prop(prop.err) })
  }
  ir <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))

  lhs <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_true("TVCL_ODE" %in% lhs)          # the model's own variable survives
  # Numbered off the SUFFIXED form, so the last resort is still recognisable as a
  # carrier. Numbering the source name would give `CL_1`, which is the spelling
  # this naming scheme exists to avoid.
  expect_true(any(grepl("^TVCL_ODE_[0-9]+$", lhs)))
  expect_match(ir$warnings, "had to be numbered", all = FALSE)
})

test_that("a self-reference to a plain local resolves to the emitted name", {
  # Kills the mutation that stops seeding rhs_map. .normalise_expr() leaves an
  # unmapped symbol untouched, so `k <- k * 2` emitted a bare lower-case `k`
  # that is declared nowhere and which ferx silently reads as a covariate.
  ini <- rbind(theta_row("t.CL", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- k * 2), quote(cl <- t.CL * exp(eta1)),
              ddt("central", quote(-cl * central)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
  rhs <- vapply(ir$indiv_params, function(p) p$rhs, "")
  expect_false(any(grepl("\\bk\\b", rhs)))
})

test_that("a $THETA label that is not an identifier is rejected", {
  # Kills the mutation reverting the label check to whitespace-only. `; CL/F` is
  # the standard NONMEM label for apparent clearance; emitting `theta CL/F(...)`
  # gives E_PARSE, and it diverges from the name references resolve to.
  ini <- rbind(theta_row("t.CL", 1), eta_row("eta1", 0.09, 1L))
  ini$label <- c("CL/F", NA_character_)
  lst <- list(quote(cl <- t.CL * exp(eta1)), ddt("central", quote(-cl * central)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
  nms <- vapply(ir$thetas, function(t) t$name, "")
  expect_true(all(grepl("^[A-Za-z][A-Za-z0-9_]*$", nms)))
})

test_that("a rename avoids a covariate name, not just a parameter name", {
  # ferx resolves theta before covariate too, so renaming onto a data column
  # reintroduces the shadowing on a new pair, with no diagnostic from the
  # engine. Driven end to end: asserting on .deshadow_theta_names() alone would
  # not prove rxui_to_ir() actually passes the covariates in, which is the part
  # that can break.
  expect_equal(unname(.deshadow_theta_names("CL", "CL", reserved = "TVCL")$map),
               "THETA_CL")
  expect_equal(unname(.deshadow_theta_names("CL", "CL")$map), "TVCL")

  # `TVCL` here is a data column: referenced, never assigned, not in iniDf.
  ini <- rbind(theta_row("t.CL", 1), eta_row("eta1", 0.09, 1L))
  ini$label <- c("CL", NA_character_)
  lst <- list(quote(cl <- t.CL * (TVCL/70) * exp(eta1)),
              ddt("central", quote(-cl * central)))
  expect_true("TVCL" %in% .covariate_names(lst, .norm_map_from_ini(ini)))

  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
  nms <- vapply(ir$thetas, function(t) t$name, "")
  expect_false("TVCL" %in% toupper(nms))
})

test_that("an aux var is inlined with the name it was normalised under", {
  # Kills the pass-2 mutation. A dotted local uppercases to `C.2` while aux_vars
  # holds `C_2`, so comparing raw symbols missed it: the var was not marked
  # auxiliary and got emitted referencing a name the block never declares.
  ini <- rbind(theta_row("t.CL", 1), theta_row("t.V", 10),
               eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.CL * exp(eta1)), quote(v <- t.V),
              quote(c.2 <- central/v), quote(eff <- c.2 * 3),
              ddt("central", quote(-cl * central)),
              ddt("resp",    quote(eff - resp)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  lhs <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_false("EFF" %in% lhs)            # state-dependent: an ODE intermediate
  # Phase 5b emits it instead of substituting it away, so the name now DOES
  # appear -- and the property to hold is the one the old assertion stood in for:
  # never referenced undeclared. `C_2` is the normalised spelling of `c.2`, so
  # this also still pins that the emitted reference and the emitted declaration
  # agree on which normalisation was applied.
  intermediates <- vapply(Filter(function(o) identical(o$kind, "assign"), ir$odes),
                          function(o) o$lhs, "")
  expect_true("C_2" %in% intermediates)
  expect_true(all(c("C_2", "EFF") %in%
                  c(intermediates, .ode_states(ir$odes),
                    vapply(ir$indiv_params, function(p) p$lhs, ""))))
})

test_that("an error model whose sigma arrives via a dotted name is classified", {
  # Kills the pass-3 mutation. Detection and classification must both see the
  # normalised expression; handing the raw one to the classifier made detection
  # succeed and classification find no sigma, emitting `DV ~ proportional()`.
  # The sigma must be one whose raw and normalised forms DIFFER, or the test
  # cannot tell the two apart: `eps.1` uppercases to EPS.1 but normalises to
  # EPS_1, which is the name the sigma is declared under.
  ini <- rbind(theta_row("t.CL", 1), eta_row("eta1", 0.09, 1L),
               sigma_row("eps.1", 0.04))
  lst <- list(quote(cl <- t.CL * exp(eta1)),
              ddt("central", quote(-cl * central)),
              quote(y <- central * (1 + eps.1)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
  expect_length(ir$error_model, 1L)
  expect_equal(ir$error_model[[1]]$type, "proportional")
  expect_equal(ir$error_model[[1]]$params, "EPS_1")
})

test_that("a theta renamed more than once is reported by its final name only", {
  # Renaming one theta can reveal another individual parameter, so a theta may
  # be renamed again in a later round. Reporting each hop named intermediate
  # values that appear nowhere in the output -- and in the aux-var-flip case the
  # intermediate ends up as an individual PARAMETER, so the message pointed at
  # the wrong kind of thing entirely.
  ini <- rbind(theta_row("th1", 1), theta_row("th2", 2), theta_row("th3", 3),
               eta_row("eta1", 0.09, 1L))
  ini$label <- c("CENTRAL", "CENTRAL", "CL", NA_character_)
  lst <- list(quote(tvcl <- th2 * 2), quote(cl <- th3 * exp(eta1)),
              ddt("central", quote(-cl * central)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  nms  <- vapply(ir$thetas, function(t) t$name, "")
  lhs  <- vapply(ir$indiv_params, function(p) p$lhs, "")
  renames <- grep("renamed to", ir$warnings, value = TRUE)
  # One line per theta, and every name it points at must be a theta in the
  # output -- never an individual parameter, never a value that vanished.
  expect_length(renames, length(unique(renames)))
  for (w in renames) {
    to <- sub(".*renamed to '([^']+)'.*", "\\1", w)
    expect_true(to %in% nms)
    expect_false(to %in% setdiff(lhs, nms))
  }
  expect_equal(intersect(toupper(nms), toupper(lhs)), character())
})

test_that("the invariant covers individual parameters the linCmt passthrough adds", {
  # The passthrough appends to indiv_params AFTER the de-shadow loop, so an
  # assertion placed before it could not see the `V = V` self-shadow the
  # passthrough comment says it exists to prevent.
  skip_if_not_installed("nonmem2rx")
  local_mocked_bindings(
    .deshadow_theta_names = function(theta_names, indiv_names, reserved = character())
      list(map = rep(NA_character_, length(theta_names)),
           reasons = vector("list", length(theta_names)), warnings = character())
  )
  ir <- suppressWarnings(rxui_to_ir(
    nonmem2rx::nonmem2rx(nm_path("pk_1cmt_oral_ampsim.ctl")), source_format = "nonmem"))
  clash <- intersect(toupper(vapply(ir$thetas, function(t) t$name, "")),
                     toupper(vapply(ir$indiv_params, function(p) p$lhs, "")))
  expect_gt(length(clash), 0L)                       # the corpus really does clash
  expect_match(ir$unsupported, "name collision", all = FALSE)
  for (nm in clash) expect_match(ir$unsupported, nm, all = FALSE)
})

# -- identifier sanitisation --------------------------------------------------

test_that(".ferx_ident maps any name onto the ferx grammar", {
  # Illegal characters become underscores.
  expect_equal(.ferx_ident("c.RTOT"), "c_RTOT")
  expect_equal(.ferx_ident("A-B+C"), "A_B_C")
  expect_equal(.ferx_ident("has space"), "has_space")
  # Case is preserved -- .norm() is where uppercasing happens.
  expect_equal(.ferx_ident("central"), "central")
  # A leading digit cannot be substituted away, and an empty name has nothing
  # to substitute; both need a prefix instead.
  expect_equal(.ferx_ident("2CPT"), "X_2CPT")
  expect_equal(.ferx_ident(""), "X")
  # Whatever goes in, a legal identifier comes out.
  for (nm in c("c.RTOT", "A-B+C", "2CPT", "", "_x", "9", "..", "a b c"))
    expect_true(.is_ferx_ident(.ferx_ident(nm)), info = nm)
})

test_that(".is_ferx_ident matches the ferx grammar", {
  expect_true(all(.is_ferx_ident(c("A", "_a", "A1", "a_B_9"))))
  expect_false(any(.is_ferx_ident(c("c.RTOT", "1A", "", "a-b", "a b"))))
})

test_that(".norm is .ferx_ident plus the uppercase convention", {
  # The old .norm() substituted only the dot. Anything it used to do it must
  # still do, or every emitted name changes at once.
  expect_equal(.norm("t.CL"), "T_CL")
  expect_equal(.norm("central"), "CENTRAL")
  expect_equal(.norm("c.RTOT"), "C_RTOT")
  # ...and it now also covers the characters the dot rule missed.
  expect_equal(.norm("2CPT"), "X_2CPT")
})

test_that("a dotted state name is renamed at every reference site", {
  # nonmem2rx prefixes `c.` onto a compartment whose name collides with a
  # variable, and the result appeared in four places: obs_cmt=, states=[...],
  # the d/dt target, and inlined into other ODE right-hand sides. Renaming the
  # declaration alone leaves the references pointing at nothing.
  ini <- rbind(theta_row("t.KEL", 0.05), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(kel <- t.KEL * exp(eta1)),
              ddt("CENT", quote(-kel * CENT - `c.RTOT`)),
              ddt("c.RTOT", quote(-kel * `c.RTOT`)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_equal(ir$structural$states, c("CENT", "c_RTOT"))
  states <- vapply(ir$odes, function(o) o$state, "")
  expect_equal(states, c("CENT", "c_RTOT"))
  # The reference inlined into the OTHER state's RHS is the whole point.
  expect_match(ir$odes[[1]]$rhs, "c_RTOT", fixed = TRUE)
  # No CODE line may still carry the dot. Comment lines may and now do -- the
  # `# renamed:` provenance line names the source spelling on purpose -- so the
  # assertion is scoped to the lines ferx actually parses.
  code <- grep("^\\s*#", strsplit(emit_ferx(ir), "\n")[[1]], invert = TRUE, value = TRUE)
  expect_false(any(grepl("c.RTOT", code, fixed = TRUE)))
  expect_match(ir$warnings, "^INFO  \\| state 'c\\.RTOT' is not a legal ferx",
               all = FALSE)
})

test_that("obs_cmt taken from ui$central is renamed with the state", {
  # `ui$central` is a separate, RAW channel for the observed compartment -- it
  # bypasses the odes entirely, so it needs its own translation. Every other
  # test reaches obs_cmt through the tail(state_names) guess, which is already
  # sanitised, leaving this branch unexercised: a mutation that dropped the
  # lookup emitted `obs_cmt=c.RTOT` beside `states=[c_RTOT]` and the whole suite
  # still passed.
  ini <- rbind(theta_row("t.KEL", 0.05), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(kel <- t.KEL * exp(eta1)),
              ddt("c.RTOT", quote(-kel * `c.RTOT`)))
  ui  <- c(mock_ui(ini, lst), list(central = "c.RTOT"))
  ir  <- suppressWarnings(rxui_to_ir(ui))

  expect_equal(ir$structural$obs_cmt, "c_RTOT")
  # And no guess was needed, so the guess warning must be absent -- otherwise
  # this could pass by falling through to the fallback path.
  expect_length(grep("obs_cmt could not be inferred", ir$warnings), 0L)
  # Code lines only -- the `# renamed:` provenance comment names the source
  # spelling deliberately.
  code <- grep("^\\s*#", strsplit(emit_ferx(ir), "\n")[[1]], invert = TRUE, value = TRUE)
  expect_false(any(grepl("c.RTOT", code, fixed = TRUE)))
})

test_that("obs_cmt is resolved through the state renames, not the parameter map", {
  # The previous test cannot separate the two maps: `c.RTOT` is a key in both
  # and they agree on it. Here they disagree. `ui$central` names a state whose
  # raw spelling is also an iniDf key, so resolving through `name_map` rewrites
  # obs_cmt to the THETA's emitted name -- producing an obs_cmt that names no
  # state at all, which ferx rejects with
  # `E_PARSE: Observable compartment 'VC' not in states [...]`.
  ini <- rbind(theta_row("CENT", 1), eta_row("eta1", 0.09, 1L))
  ini$label <- c("VC", NA_character_)
  lst <- list(quote(k <- CENT * exp(eta1)), ddt("CENT", quote(-k * CENT)))
  ui  <- c(mock_ui(ini, lst), list(central = "CENT"))
  ir  <- suppressWarnings(rxui_to_ir(ui))

  expect_equal(ir$structural$obs_cmt, "CENT")
  # The invariant that actually matters: obs_cmt must name a declared state.
  expect_true(ir$structural$obs_cmt %in% ir$structural$states)
})

test_that("a state rename is recorded in the emitted file as provenance", {
  # The .ferx is the artefact that gets shared. Without this line a reader
  # holding only that file cannot map `c_RTOT` back to the $MODEL compartment or
  # the A(n) index it came from -- the rename lived only in result$warnings,
  # which does not travel with the file.
  ini <- rbind(theta_row("t.KEL", 0.05), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(kel <- t.KEL * exp(eta1)),
              ddt("c.RTOT", quote(-kel * `c.RTOT`)))
  # `central` is supplied so obs_cmt needs no guess: the rename is then the ONLY
  # note the translation produces, which is what makes the warning-count
  # assertion below meaningful rather than incidental.
  ir  <- suppressWarnings(rxui_to_ir(c(mock_ui(ini, lst), list(central = "c.RTOT"))))
  txt <- emit_ferx(ir)

  expect_equal(ir$state_renames, c("c.RTOT" = "c_RTOT"))
  expect_match(txt, "# renamed: state c.RTOT -> c_RTOT", fixed = TRUE)
  # Provenance, not a diagnostic: it must not inflate the warning count, and a
  # model whose only note is a rename must not advertise warnings at all.
  expect_false(grepl("# Warnings:", txt, fixed = TRUE))
})

test_that("a model with no renames gets no provenance line", {
  ini <- rbind(theta_row("t.KA", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(ka <- t.KA * exp(eta1)), ddt("depot", quote(-ka * depot)))
  txt <- emit_ferx(suppressWarnings(rxui_to_ir(mock_ui(ini, lst))))
  expect_false(grepl("# renamed:", txt, fixed = TRUE))
})

test_that("a state that is already legal is left alone", {
  # A rename is a user-visible change to a name they index by. Only names that
  # must change may change -- `depot`/`central` are legal ferx identifiers and
  # uppercasing them buys nothing.
  ini <- rbind(theta_row("t.KA", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(ka <- t.KA * exp(eta1)),
              ddt("depot", quote(-ka * depot)),
              ddt("central", quote(ka * depot)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_equal(vapply(ir$odes, function(o) o$state, ""), c("depot", "central"))
  expect_length(grep("state '", ir$warnings), 0L)
})

test_that("two states that sanitise to the same name are disambiguated", {
  # `.ferx_ident()` is not injective: distinct illegal names can collapse onto
  # one legal one, which would silently merge two compartments.
  ini <- rbind(theta_row("t.K", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- t.K * exp(eta1)),
              ddt("A.B", quote(-k * `A.B`)),
              ddt("A-B", quote(-k * `A-B`)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  states <- vapply(ir$odes, function(o) o$state, "")
  expect_length(unique(states), 2L)
  expect_true(all(.is_ferx_ident(states)))
  # Each RHS must reference its OWN state, not the one it collided with.
  for (o in ir$odes) expect_match(o$rhs, o$state, fixed = TRUE)
})

test_that("an already-legal state keeps its name against a sanitised collider", {
  # Processing in source order let an ILLEGAL name claim the spelling of a
  # legal one that appeared later: `A.B` sanitised to `A_B` and displaced the
  # existing `A_B` to `A_B_1`, renaming the one name in the pair that was
  # already fine. The legal name has first claim.
  ini <- rbind(theta_row("t.K", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- t.K * exp(eta1)),
              ddt("A.B", quote(-k * `A.B`)),
              ddt("A_B", quote(-k * A_B)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  states <- vapply(ir$odes, function(o) o$state, "")
  # Source order is preserved in the odes; the legal name is untouched.
  expect_equal(states, c("A_B_1", "A_B"))
  expect_length(grep("state 'A_B' renamed", ir$warnings), 0L)
  # Each RHS still references its own state.
  for (o in ir$odes) expect_match(o$rhs, o$state, fixed = TRUE)
})

test_that("a state colliding with an eta name is renamed, not the eta", {
  ini <- rbind(theta_row("t.K", 1), eta_row("ETA1", 0.09, 1L))
  lst <- list(quote(k <- t.K * exp(ETA1)), ddt("eta1", quote(-k * eta1)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  # The eta keeps its name; the state gives way. Critically, the eta reference
  # in [individual_parameters] must NOT have been rewritten along with it.
  expect_equal(vapply(ir$omegas, function(o) o$names, ""), "ETA1")
  expect_match(ir$indiv_params[[1]]$rhs, "exp(ETA1)", fixed = TRUE)
  expect_false(toupper(ir$odes[[1]]$state) == "ETA1")
})

test_that("state renaming does not rename thetas -- one owner per collision", {
  # `.deshadow_theta_names()` is the single owner of theta naming (CLAUDE.md).
  # Reserving theta names in the state sanitiser too made two owners for one
  # collision: a state `central` beside a theta labelled CENTRAL got renamed to
  # `central_1` while the theta was independently renamed, churning both names.
  ini <- rbind(theta_row("th1", 1), eta_row("eta1", 0.09, 1L))
  ini$label <- c("CENTRAL", NA_character_)
  lst <- list(quote(cl <- th1 * exp(eta1)), ddt("central", quote(-cl * central)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_equal(.ode_states(ir$odes), "central")
  expect_length(grep("state 'central'", ir$warnings), 0L)
})

test_that("a covariate whose name is illegal is reported, never renamed", {
  # ferx matches covariates to data columns by exact name, case included, so a
  # rename cannot fix an illegal covariate -- it only moves the failure from
  # E_PARSE to E_MISSING_COVARIATE. It has to be said out loud instead.
  ini <- rbind(theta_row("t.CL", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.CL * `WT.BASE` * exp(eta1)),
              ddt("central", quote(-cl * central)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_match(ir$warnings, "^ERROR \\| covariate reference", all = FALSE)
  expect_match(ir$unsupported, "covariate name is not a legal ferx identifier",
               all = FALSE)
  expect_match(ir$unsupported, "WT.BASE", all = FALSE, fixed = TRUE)
  # It must reach the artefact, not only the result object.
  expect_match(emit_ferx(ir), "# WARNING: covariate reference")
})

test_that("covariate case is preserved exactly as written", {
  # ferx matches data columns case-SENSITIVELY (ferx-core datareader.rs uses
  # case-insensitive matching only for the standard columns). nonmem2rx keeps
  # the $INPUT case for data items while lowercasing assigned variables, and a
  # covariate is absent from name_map so .normalise_expr() leaves it alone --
  # which is how this works today. Nothing tested it, so a blanket uppercase
  # pass over emitted symbols would have broken every covariate silently.
  ini <- rbind(theta_row("t.CL", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.CL * (Wt/70) * exp(eta1)),
              ddt("central", quote(-cl * central)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_match(ir$indiv_params[[1]]$rhs, "Wt/70", fixed = TRUE)
  expect_false(grepl("WT/70", ir$indiv_params[[1]]$rhs, fixed = TRUE))
  # A legal covariate is not reported as a problem.
  expect_length(grep("covariate reference", ir$warnings), 0L)
})

test_that("the emitted state name does not depend on statement order", {
  # The d/dt target used to be resolved through `name_map`, which
  # .parse_model_exprs() keeps extending as it walks -- every ordinary
  # assignment installs an alias. So an unrelated `central <- 0` written ABOVE
  # d/dt(central) renamed the state to CENTRAL, and the same model with that
  # line BELOW did not. The renames are decided once, before parsing, and the
  # declaration reads only those.
  #
  # The RHS is asserted, not only the state. Comparing `o$state` alone passed
  # while the equations still differed: the declaration was already read from
  # the pre-parse decision, but the references were not, so `central <- 0`
  # standing above rebound the state's key mid-walk and emitted
  # `d/dt(central_1) = -K * CENTRAL` -- a derivative that reads the constant
  # individual parameter and never its own compartment, which ferx validates
  # clean.
  ini <- rbind(theta_row("t.K", 1), eta_row("eta1", 0.09, 1L))
  ode <- function(l) {
    o <- suppressWarnings(rxui_to_ir(mock_ui(ini, l)))$odes
    vapply(o, function(x) paste0("d/dt(", x$state, ") = ", x$rhs), "")
  }
  above <- ode(list(quote(k <- t.K * exp(eta1)), quote(central <- 0),
                    ddt("central", quote(-k * central))))
  below <- ode(list(quote(k <- t.K * exp(eta1)),
                    ddt("central", quote(-k * central)), quote(central <- 0)))
  expect_equal(above, below)
  # and the surviving equation is the right one -- the state, not the constant.
  expect_equal(above, "d/dt(central_1) = -K * central_1")
})

test_that("a state whose raw name is a parameter key keeps its own name", {
  # Same root cause, different symptom: `name_map` holds every iniDf key, so a
  # state that merely happened to share one was renamed to that PARAMETER's
  # emitted name -- a theta keyed CENT and labelled VC turned d/dt(CENT) into
  # d/dt(VC), silently making the compartment the theta. The sanitiser had
  # decided no such rename.
  ini <- rbind(theta_row("CENT", 1), eta_row("eta1", 0.09, 1L))
  ini$label <- c("VC", NA_character_)
  lst <- list(quote(k <- CENT * exp(eta1)), ddt("CENT", quote(-k * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_equal(vapply(ir$odes, function(o) o$state, ""), "CENT")
  expect_equal(vapply(ir$thetas, function(t) t$name, ""), "VC")

  # The declaration was the only half this test originally checked, and the
  # REFERENCES are where the same map does the real damage: the RHS `-k * CENT`
  # resolved CENT through `name_map` to the theta's emitted name and emitted
  # `-K * VC`. That parses -- VC is a declared theta -- so the model runs, with
  # the elimination term reading a fixed effect instead of the amount in the
  # compartment. Nothing downstream distinguishes the two.
  #
  # Inside [odes] the state wins; in $PK the same symbol still means the theta.
  expect_equal(ir$odes[[1]]$rhs, "-K * CENT")
  expect_false(grepl("VC", ir$odes[[1]]$rhs, fixed = TRUE))
  expect_equal(ir$indiv_params[[1]]$rhs, "VC * exp(ETA1)")
  # ...and the ambiguity itself is reported, since the source cannot say which
  # reading it meant.
  expect_match(ir$unsupported, "state/parameter source-name collision: CENT",
               all = FALSE, fixed = TRUE)
})

test_that("a sanitised state never swallows an individual parameter", {
  # The worst failure this review found, because ferx accepts the result. When a
  # sanitised state lands on an assignment target, that assignment is absorbed
  # into aux_vars, dropped from [individual_parameters], and its references
  # resolve to the state: `A_B = K * exp(ETA1)` beside `d/dt(A.B) = -A_B * A.B`
  # emitted `d/dt(A_B) = -A_B * A_B` -- the amount squared, rate constant and
  # IIV gone -- and validated clean.
  ini <- rbind(theta_row("t.K", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(A_B <- t.K * exp(eta1)), ddt("A.B", quote(-A_B * `A.B`)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  lhs <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_true("A_B" %in% lhs)                       # the parameter survives
  expect_false(identical(ir$odes[[1]]$state, "A_B")) # the state gave way
  # and the IIV is still wired to something
  expect_match(paste(vapply(ir$indiv_params, function(p) p$rhs, ""), collapse = " "),
               "ETA1", fixed = TRUE)
})

# NOTE: the state/individual-parameter invariant in rxui_to_ir() has no test
# because I could not construct an input that reaches it. Once the sanitiser
# reserves every assignment target in `lst`, a state can no longer be renamed
# onto an individual parameter, and the parameters that appear later (the linCmt
# passthrough) only exist for models with no ODEs. It is kept as defence in
# depth -- ferx is silent about this collision and the symptom is a deleted
# parameter -- but it is unasserted, and a test that merely built the colliding
# IR by hand would assert the constructor, not the guard.

test_that("a symbol referenced before its assignment gets the emitted name", {
  # The alias used to be installed only once the walk had passed the assignment,
  # so a forward reference was emitted EXACTLY as written -- two spellings of one
  # variable, `f.rac` illegal at every reference site and `F_RAC` declared beside
  # it. Nothing reported it either: .unmapped_symbols() subtracts assignment
  # targets regardless of position, so the leaked name was filtered out of the
  # legality check as "assigned, therefore legal by construction".
  ini <- rbind(theta_row("t.CL", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.CL * exp(eta1)),
              ddt("central", quote(-cl * central * `f.rac`)),
              quote(f.rac <- 0.5))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_equal(ir$odes[[1]]$rhs, "-CL * central * F_RAC")
  expect_true("F_RAC" %in% vapply(ir$indiv_params, function(p) p$lhs, ""))
  code <- grep("^\\s*#", strsplit(emit_ferx(ir), "\n")[[1]], invert = TRUE, value = TRUE)
  expect_false(any(grepl("f.rac", code, fixed = TRUE)))
})

test_that("the order-independent alias does not override an iniDf binding", {
  # Seeding every assignment target up front must not disturb de-shadowing: in
  # `cl <- cl * exp(eta.cl)` the RHS `cl` is the THETA and the following
  # `k20 <- cl/v` is the individual parameter, and only the mid-walk rebinding
  # tells them apart. The seed therefore skips names the map already holds.
  ini <- rbind(theta_row("cl", 1), theta_row("v", 10), eta_row("eta.cl", 0.09, 1L))
  lst <- list(quote(cl <- cl * exp(eta.cl)), quote(v <- v), quote(k20 <- cl / v),
              ddt("central", quote(-k20 * central)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  lhs <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_match(ir$indiv_params[[match("CL", lhs)]]$rhs, "exp(ETA_CL)", fixed = TRUE)
  # the theta was renamed away, and K20 reads the individual parameter
  expect_equal(ir$indiv_params[[match("K20", lhs)]]$rhs, "CL/V")
})

test_that("a theta named after a ferx solver builtin is renamed", {
  # The quietest of the three theta failures. ferx does not reject `theta TIME`
  # -- it resolves the bare name to the value the solver injects, so KA reads the
  # integrator clock, the theta is declared and estimated and never referenced,
  # and the only diagnostic is a W_UNUSED_PARAM that explains none of it.
  ini <- rbind(theta_row("t.KA", 0.1), eta_row("eta1", 0.09, 1L))
  ini$label <- c("TIME", NA_character_)
  lst <- list(quote(ka <- t.KA * exp(eta1)), ddt("depot", quote(-ka * depot)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_equal(vapply(ir$thetas, function(t) t$name, ""), "TVTIME")
  expect_equal(ir$indiv_params[[1]]$rhs, "TVTIME * exp(ETA1)")
  expect_match(ir$warnings, "collides with a ferx solver builtin", all = FALSE)
  # Every name in the list, not just TIME -- T and TAD are ordinary $DES and
  # $PK variable names, so the whole set has to be covered.
  for (nm in .RESERVED_ODE_NAMES) {
    i2 <- rbind(theta_row("t.KA", 0.1), eta_row("eta1", 0.09, 1L))
    i2$label <- c(nm, NA_character_)
    r2 <- suppressWarnings(rxui_to_ir(mock_ui(i2, lst)))
    expect_false(toupper(r2$thetas[[1]]$name) %in% .RESERVED_ODE_NAMES, info = nm)
  }
})

test_that("an individual parameter named after a builtin is reported", {
  # It cannot be renamed the way a theta can -- the name is the source's, and
  # [scaling] and [error_model] reference it by that name -- so it is reported.
  # ferx-core rejects the model outright, but names no source variable.
  ini <- rbind(theta_row("t.K", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(TAD <- t.K * exp(eta1)), ddt("depot", quote(-TAD * depot)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_match(ir$warnings, "^ERROR \\| TAD names both an individual parameter",
               all = FALSE)
  expect_match(ir$unsupported, "individual parameter collides with a ferx builtin",
               all = FALSE)
  expect_match(emit_ferx(ir), "# WARNING: TAD names both an individual parameter")
})

test_that("a sigma whose source name needs sanitising still reaches the error model", {
  # nonmem2rx keeps sigma out of iniDf, so nothing bound its source spelling to
  # its emitted name: the declaration and the eps reference in $ERROR agreed only
  # while both were plain uppercase. Once the declaration was sanitised the two
  # diverged, pass 3 found no sigma in the error assignment, and the file came
  # out with `sigma EPS_1` declared and NO [error_model] block at all -- a model
  # with no residual error and not one word about it.
  ini <- rbind(theta_row("t.K", 1), eta_row("eta1", 0.09, 1L))
  sig <- matrix(0.04, 1, 1, dimnames = list("eps.1", "eps.1"))
  ui  <- c(mock_ui(ini, list(quote(k <- t.K * exp(eta1)),
                             ddt("cent", quote(-k * cent)),
                             quote(y <- cent * (1 + `eps.1`)))),
           list(sigma = sig))
  ir  <- suppressWarnings(rxui_to_ir(ui))

  expect_equal(vapply(ir$sigmas, function(s) s$name, ""), "EPS_1")
  expect_length(ir$error_model, 1L)
  expect_equal(ir$error_model[[1]]$params, "EPS_1")
  # and the sigma is not misreported as an illegal covariate on the way through
  expect_length(grep("covariate reference", ir$warnings), 0L)
})

test_that("a non-scalar ui$central does not abort the translation", {
  # is.character() is TRUE for character(0) and for length 2, so both reached an
  # `if` that errors ("argument is of length zero" / "the condition has length
  # > 1") on a model that previously translated.
  ini <- rbind(theta_row("t.K", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- t.K * exp(eta1)), ddt("cent", quote(-k * cent)))
  for (v in list(character(0), c("a", "b"))) {
    ui <- c(mock_ui(ini, lst), list(central = v))
    expect_no_error(suppressWarnings(rxui_to_ir(ui)))
  }
})

test_that(".ferx_ident maps NA to a legal identifier", {
  # nzchar(NA) is TRUE and gsub/sub pass NA through, so NA slipped past both
  # guards and out of a function documented to return a legal identifier for
  # any input.
  expect_true(.is_ferx_ident(.ferx_ident(NA_character_)))
  expect_false(is.na(.ferx_ident(NA_character_)))
})

test_that("an illegal free symbol is reported even when it normalises to a known name", {
  # The check used to run on the covariate set, which is classified by
  # NORMALISED name: a raw `c.RTOT` beside an eta named `c_RTOT` normalises to
  # a known `C_RTOT` and was filtered out as "not a covariate" -- while
  # .normalise_expr(), which matches raw keys, still emitted `c.RTOT` verbatim.
  ini <- rbind(theta_row("t.K", 1), eta_row("c_RTOT", 0.09, 1L))
  lst <- list(quote(k <- t.K * `c.RTOT` * exp(c_RTOT)),
              ddt("central", quote(-k * central)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_match(ir$unsupported, "not a legal ferx identifier", all = FALSE)
  expect_match(ir$unsupported, "c.RTOT", all = FALSE, fixed = TRUE)
})

test_that("a state whose source name is also a parameter name is scope-resolved", {
  # Both deparse to the same symbol, so the source cannot say which is meant and
  # the collision has to be reported. What it must NOT do is corrupt the model:
  # refusing the rename left the state sharing a name with the eta, so the
  # assignment referencing it was absorbed into aux_vars, dropped from
  # [individual_parameters] -- the block came out EMPTY -- and self-inlined to
  # the depth cap, emitting exp(ETA_X) fifteen times.
  ini <- rbind(theta_row("t.K", 1), eta_row("ETA_X", 0.09, 1L))
  lst <- list(quote(k <- t.K * exp(ETA_X)), ddt("ETA_X", quote(-k * ETA_X)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_match(ir$warnings, "^ERROR \\| 'ETA_X' names both an ODE state and a",
               all = FALSE)
  expect_match(ir$unsupported, "state/parameter source-name collision",
               all = FALSE)

  # Resolved by scope, which is how ferx itself reads them: inside [odes] the
  # bare name is the state (etas are out of scope there), outside it the eta.
  expect_equal(ir$odes[[1]]$state, "ETA_X_1")
  expect_equal(ir$odes[[1]]$rhs, "-K * ETA_X_1")
  lhs <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_true("K" %in% lhs)
  expect_match(ir$indiv_params[[match("K", lhs)]]$rhs, "exp(ETA_X)", fixed = TRUE)
  # One exp(), not the depth-cap cascade.
  expect_equal(lengths(regmatches(ir$odes[[1]]$rhs,
                                  gregexpr("exp(", ir$odes[[1]]$rhs, fixed = TRUE))),
               0L)
})

# -- $MODEL DEFOBS ------------------------------------------------------------

test_that(".extract_nm_defobs reads DEFOBS wherever it is declared", {
  f <- file.path(tmp_ctl_dir(), "defobs.ctl")
  writeLines(c("$PROBLEM x", "$MODEL",
               "  COMP=(CENT, DEFDOSE, DEFOBS)", "  COMP=(PERIPH)", "$PK"), f)
  out <- .extract_nm_defobs(f)
  expect_equal(out$index, 1L)
  expect_equal(out$name, "CENT")
  expect_equal(out$n_comp, 2L)
})

test_that(".extract_nm_defobs tolerates the legal $MODEL spellings", {
  n  <- 0L
  mk <- function(...) {
    n <<- n + 1L
    f <- file.path(tmp_ctl_dir(), paste0("model", n, ".ctl"))
    writeLines(c("$PROBLEM x", "$MODEL", ..., "$PK"), f); f
  }
  # `COMP (...)` without `=`, space-separated attributes.
  expect_equal(.extract_nm_defobs(mk("  COMP (DEPOT DEFDOSE)",
                                     "  COMP (CENTRAL DEFOBS)"))$index, 2L)
  # DEFOBSERVATION, the unabbreviated spelling.
  expect_equal(.extract_nm_defobs(mk("  COMP=(A)",
                                     "  COMP=(B, DEFOBSERVATION)"))$index, 2L)
  # DEFDOSE must not be mistaken for DEFOBS -- they share the DEF prefix.
  # `index` NA rather than a NULL return: the COMP list is worth having on its
  # own, since it is the only way to check a NONMEM compartment NUMBER against
  # the d/dt order it indexes into. Callers already require a non-NA index.
  no_defdose <- .extract_nm_defobs(mk("  COMP=(A, DEFDOSE)", "  COMP=(B)"))
  expect_true(is.na(no_defdose$index))
  expect_true(is.na(no_defdose$name))
  expect_equal(no_defdose$comps, c("A", "B"))
  # A DEFOBS that only appears in a comment is not a declaration.
  commented <- .extract_nm_defobs(mk("  COMP=(A) ; DEFOBS goes here one day",
                                     "  COMP=(B)"))
  expect_true(is.na(commented$index))
  expect_equal(commented$comps, c("A", "B"))
  # A bare `COMP=NAME` occupies an ordinal, so it must appear in the list or
  # every later ordinal is off by one.
  bare <- .extract_nm_defobs(mk("  COMP=DEPOT", "  COMP=(CENTRAL, DEFOBS)"))
  expect_equal(bare$index, 2L)
  expect_equal(bare$comps, c("DEPOT", "CENTRAL"))
  # No $MODEL at all -- still NULL, since there is no COMP list either.
  f <- file.path(tmp_ctl_dir(), "nomodel.ctl")
  writeLines(c("$PROBLEM x", "$PK"), f)
  expect_null(.extract_nm_defobs(f))
})

# -- $MODEL DEFDOSE and the $INPUT CMT item (issue #27) -----------------------

test_that(".extract_nm_defobs reads DEFDOSE in every legal abbreviation", {
  n  <- 0L
  mk <- function(...) {
    n <<- n + 1L
    f <- file.path(tmp_ctl_dir(), paste0("dd", n, ".ctl"))
    writeLines(c("$PROBLEM x", "$MODEL", ..., "$PK"), f); f
  }
  # DEFDOSE and DEFOBSERVATION diverge at the fourth character, so DEFD is the
  # shortest unambiguous prefix. `^DEFOBS` once rejected the legal `DEFO`, which
  # is the same bug one attribute over -- cover the short spellings here so it
  # cannot repeat.
  for (sp in c("DEFD", "DEFDOS", "DEFDOSE"))
    expect_equal(.extract_nm_defobs(mk("  COMP=(A)",
                                       paste0("  COMP=(B, ", sp, ")")))$defdose,
                 2L, info = sp)
  # DEFOBS must not answer the DEFDOSE question.
  only_obs <- .extract_nm_defobs(mk("  COMP=(A)", "  COMP=(B, DEFOBS)"))
  expect_true(is.na(only_obs$defdose))
  expect_equal(only_obs$index, 2L)
  # No attributes anywhere. NA, not 1: NONMEM's own default IS compartment 1 and
  # so is ferx's, so the two agree and there is nothing to report -- filling the
  # NA in here would turn that agreement into a claim the source never made.
  expect_true(is.na(.extract_nm_defobs(mk("  COMP=(A)", "  COMP=(B)"))$defdose))
  # Both attributes on different compartments, which is the shape #27 needs.
  both <- .extract_nm_defobs(mk("  COMP=(PERIPH)",
                                "  COMP=(CENTRAL, DEFDOSE, DEFOBS)"))
  expect_equal(both$defdose, 2L)
  expect_equal(both$index, 2L)
})

test_that(".nm_input_has_cmt reads the $INPUT data items", {
  n  <- 0L
  mk <- function(...) {
    n <<- n + 1L
    f <- file.path(tmp_ctl_dir(), paste0("inp", n, ".ctl"))
    writeLines(c("$PROBLEM x", ..., "$PK"), f); f
  }
  expect_true(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV CMT")))
  expect_false(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV")))
  # Synonym pairs, both directions: NONMEM reads the column as CMT either way.
  expect_true(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV CMT=COMPT")))
  expect_true(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV COMPT=CMT")))
  # A DROPped item is not read by NONMEM, so NM-TRAN falls back to DEFDOSE
  # exactly as if the column were absent -- which is the question being asked.
  expect_false(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV CMT=DROP")))
  expect_false(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV CMT=SKIP")))
  # Exact match on the item name. `CMTX` is a different column.
  expect_false(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV CMTX")))
  # Comma-separated items are legal.
  expect_true(.nm_input_has_cmt(mk("$INPUT ID,TIME,AMT,DV,MDV,CMT")))
  # Items may continue onto following lines -- pk_1cmt_oral.mod is written that
  # way, so reading only the $INPUT line would answer FALSE for a real model.
  expect_true(.nm_input_has_cmt(mk("$INPUT", "ID TIME AMT DV MDV CMT")))
  # A CMT that appears only in a comment is not a data item.
  expect_false(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV ; no CMT here")))
  # Whitespace around a synonym pair's `=` must not change the answer. Splitting
  # on whitespace before parsing the `=` handed `CMT = DROP` back as three
  # tokens, of which the bare `CMT` answered TRUE -- so the DROP rule held for
  # `CMT=DROP` and not for `CMT = DROP`, and NEWS.md documents a rule the code
  # applied to only one spelling.
  expect_false(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV CMT = DROP")))
  expect_false(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV CMT =DROP")))
  expect_true(.nm_input_has_cmt(mk("$INPUT ID TIME AMT DV MDV CMT = COMPT")))
  # Lower case throughout, which NM-TRAN accepts.
  expect_true(.nm_input_has_cmt(mk("$input id time amt dv mdv cmt")))
  # `$INP` is the shortest unambiguous abbreviation ($IN collides with $INFN).
  expect_true(.nm_input_has_cmt(mk("$INP ID TIME AMT DV MDV CMT")))
  # Unknown, not FALSE: with no $INPUT there is no evidence either way, and
  # reporting a divergence on no evidence is worse than staying quiet.
  expect_true(is.na(.nm_input_has_cmt(mk("$DATA x.csv"))))
  expect_true(is.na(.nm_input_has_cmt(mk("$INPUT"))))
})

test_that("the DEFDOSE/CMT divergence needs BOTH halves to fire", {
  skip_if_not_installed("nonmem2rx")
  # Neither half alone discriminates: with a CMT column the data decides and
  # DEFDOSE never applies, and with DEFDOSE on compartment 1 the wrong code
  # gives the right answer. So the guard is exercised through all four
  # combinations, not just the one that fires.
  p  <- nm_path("defdose_no_cmt.ctl")
  ui <- nonmem2rx::nonmem2rx(p)
  mk <- function(has_cmt) suppressWarnings(rxui_to_ir(
    ui, source_format = "nonmem", scaling_hint = .extract_nm_scaling(p),
    obs_hint = .extract_nm_defobs(p), has_cmt_col = has_cmt))

  # Targeted on DEFDOSE rather than on `length(unsupported) == 0`: dropping
  # obs_hint below removes the DEFOBS evidence too, and the pre-existing
  # obs_cmt-guess ERROR then fires for reasons that have nothing to do with
  # this guard. A silence assertion that counts unrelated entries fails for
  # unrelated reasons.
  quiet <- function(ir) expect_false(any(grepl("DEFDOSE", ir$unsupported)))

  fires <- mk(FALSE)
  expect_true(any(grepl("DEFDOSE", fires$unsupported)))
  expect_true(any(grepl("^ERROR.*DEFDOSE", fires$warnings)))

  # CMT column present -> the data decides -> silence.
  quiet(mk(TRUE))
  # Unknown -> silence. NA must not be read as FALSE.
  quiet(mk(NA))

  # DEFDOSE on compartment 1 -> the two rules agree -> silence, even with no
  # CMT column.
  hint <- .extract_nm_defobs(p)
  hint$defdose <- 1L
  quiet(suppressWarnings(rxui_to_ir(
    ui, source_format = "nonmem", scaling_hint = .extract_nm_scaling(p),
    obs_hint = hint, has_cmt_col = FALSE)))
  # No $MODEL evidence at all -> silence.
  quiet(suppressWarnings(rxui_to_ir(
    ui, source_format = "nonmem", scaling_hint = .extract_nm_scaling(p),
    obs_hint = NULL, has_cmt_col = FALSE)))

  # A double `defdose` must fire exactly as an integer one does. rxui_to_ir()
  # is exported and obs_hint is a documented argument, so `list(defdose = 2)`
  # -- 2, not 2L, which is how anyone writes it -- is a caller spelling, and
  # is.integer() is FALSE for it. That skipped the entire guard and restored
  # the silence it exists to break, which no assertion above could see.
  dbl <- .extract_nm_defobs(p)
  dbl$defdose <- 2
  fires_dbl <- suppressWarnings(rxui_to_ir(
    ui, source_format = "nonmem", scaling_hint = .extract_nm_scaling(p),
    obs_hint = dbl, has_cmt_col = FALSE))
  expect_true(any(grepl("DEFDOSE", fires_dbl$unsupported)))
})

test_that("the DEFDOSE test is by compartment NAME, not by $MODEL ordinal", {
  skip_if_not_installed("nonmem2rx")
  # `defdose != 1` and "DEFDOSE is not the compartment ferx puts first" are the
  # same question only while .nm_cmt_order() reconciled, because states=[...] is
  # then $MODEL COMP order. When it did not -- the case the #25 ERROR twenty
  # lines above reports -- `defdose` stays a $MODEL ordinal while states stays
  # in $DES order, and the two come apart in BOTH directions. None of this is
  # reachable through nm_to_ferx() on a bundled model, so it is driven through
  # obs_hint, which is a documented argument of an exported function.
  p  <- nm_path("defdose_no_cmt.ctl")
  ui <- nonmem2rx::nonmem2rx(p)
  fires <- function(comps, dd) {
    h <- .extract_nm_defobs(p); h$comps <- comps; h$defdose <- dd
    ir <- suppressWarnings(rxui_to_ir(
      ui, source_format = "nonmem", scaling_hint = .extract_nm_scaling(p),
      obs_hint = h, has_cmt_col = FALSE))
    any(grepl("DEFDOSE", ir$unsupported))
  }

  # Ordinal 2, and the name it picks IS the compartment ferx puts first. No
  # divergence -- and `defdose != 1L` called it one, on a message that would
  # have named PERIPH as both the wrong compartment and the right one.
  expect_false(fires(c("YYY", "PERIPH"), 2L))

  # Ordinal 1, and the name it picks is NOT the compartment ferx puts first.
  # A real divergence that `defdose != 1L` stayed silent on -- the direction
  # that matters, because silence here is the SILENT-WRONG the guard exists to
  # break.
  expect_true(fires(c("CENTRAL", "YYY"), 1L))

  # Reconciled and DEFDOSE on compartment 1: both readings agree, still quiet.
  expect_false(fires(c("CENTRAL", "PERIPH"), 1L))

  # When it does fire on unreconciled names, both names in the message must be
  # sourced honestly: the $MODEL name for what NONMEM dosed, the emitted state
  # name for where the dose lands. They are two naming universes precisely
  # because reconciliation failed, and collapsing them would hide that.
  h <- .extract_nm_defobs(p); h$comps <- c("XXX", "YYY"); h$defdose <- 2L
  ir <- suppressWarnings(rxui_to_ir(
    ui, source_format = "nonmem", scaling_hint = .extract_nm_scaling(p),
    obs_hint = h, has_cmt_col = FALSE))
  expect_match(ir$structural$note, "doses land in 'PERIPH'", fixed = TRUE)
  expect_match(ir$structural$note, "NONMEM dosed 'YYY'",     fixed = TRUE)
})

test_that("emit_ferx renders structural$note above the line it explains", {
  # Tier 1 on the emitter itself. The integration test covers this only through
  # nonmem2rx and a .ctl, so it skips entirely where that package is absent --
  # and it cannot reach the pk_macro branch or the empty-note branch at all.
  mk <- function(structural) new_ferx_ir(
    source_format = "nonmem",
    thetas = list(list(name = "TVCL", init = 1, lower = 0.001, upper = 10)),
    structural = structural,
    odes = list(list(kind = "ddt", state = "A", rhs = "-A"),
                list(kind = "ddt", state = "C", rhs = "A")))
  ode <- function(note) {
    st <- list(type = "ode", obs_cmt = "C", states = c("A", "C"))
    if (!missing(note)) st$note <- note
    strsplit(emit_ferx(mk(st)), "\n")[[1]]
  }

  ln  <- ode("WARNING: doses land in 'A' here")
  hdr <- which(ln == "[structural_model]")
  expect_length(hdr, 1L)
  expect_equal(ln[hdr + 1L], "  # WARNING: doses land in 'A' here")
  expect_equal(ln[hdr + 2L], "  ode(obs_cmt=C, states=[A, C])")

  # No note, and an empty one, both leave the block exactly as it was.
  for (ln2 in list(ode(), ode("")))
    expect_equal(ln2[which(ln2 == "[structural_model]") + 1L],
                 "  ode(obs_cmt=C, states=[A, C])")

  # The emitter is type-agnostic about where a note goes, so a pk_macro note
  # lands above its `pk` line. Nothing produces one today; this is what makes
  # the first producer's output predictable rather than discovered.
  pk <- strsplit(emit_ferx(new_ferx_ir(
    source_format = "nonmem",
    thetas = list(list(name = "TVCL", init = 1, lower = 0.001, upper = 10)),
    structural = list(type = "pk_macro", pk_call = "one_cpt_iv",
                      pk_args = list(cl = "CL", v = "V"),
                      note = "WARNING: check the dose compartment"))), "\n")[[1]]
  h2 <- which(pk == "[structural_model]")
  expect_equal(pk[h2 + 1L], "  # WARNING: check the dose compartment")
  expect_equal(pk[h2 + 2L], "  pk one_cpt_iv(cl=CL, v=V)")
})

test_that("a malformed structural$note is rejected, not passed to nzchar()", {
  # nzchar() aborts on character(0) with "missing value where TRUE/FALSE
  # needed" and on a length-2 vector with "'length = 2' in coercion to
  # 'logical(1)'" -- base-R errors naming neither the field nor the function.
  # R/ir.R records exactly this trap for obs_cmt and says it must not be
  # reintroduced; new_ferx_ir() is exported and `structural` is a user-supplied
  # list, so the note field reopened it.
  mk <- function(note) new_ferx_ir(
    source_format = "nonmem",
    thetas = list(list(name = "TVCL", init = 1, lower = 0.001, upper = 10)),
    structural = list(type = "ode", obs_cmt = "C", states = c("A", "C"),
                      note = note),
    odes = list(list(kind = "ddt", state = "A", rhs = "-A"),
                list(kind = "ddt", state = "C", rhs = "A")))
  for (bad in list(c("x", "y"), character(0), NA_character_, 1L))
    expect_error(emit_ferx(mk(bad)), "structural\\$note")
})

test_that("the observation expression outranks DEFOBS when they disagree", {
  skip_if_not_installed("nonmem2rx")
  # NONMEM's DEFOBS is the default observation compartment for records with no
  # CMT -- it is not a statement about what $ERROR reads. When $ERROR names a
  # compartment outright (`Y = A(2)/S2`), that is the model saying it, and
  # preferring DEFOBS there regressed models that were previously correct.
  p  <- nm_path("defobs_expression_wins.ctl")
  ui <- nonmem2rx::nonmem2rx(p)
  ir <- suppressWarnings(rxui_to_ir(ui, source_format = "nonmem",
                                    scaling_hint = .extract_nm_scaling(p),
                                    obs_hint     = .extract_nm_defobs(p)))
  expect_equal(.extract_nm_defobs(p)$name, "DEPOT")   # $MODEL says DEPOT
  expect_equal(ir$structural$obs_cmt, "CENT")         # $ERROR says A(2)
  expect_equal(ir$scaling$obs_scale, "V")             # and S2 = V follows it
  expect_match(ir$warnings, "\\$MODEL declares 'DEPOT' as DEFOBS", all = FALSE)
})

test_that("DEFOBS decides when the DV expression goes through F", {
  skip_if_not_installed("nonmem2rx")
  # `$ERROR IPRE = F` carries no compartment of its own -- in NONMEM `F` IS the
  # DEFOBS compartment. nonmem2rx nonetheless resolves `f <- CENTRAL` ignoring
  # DEFOBS, so following `f` would launder its guess into an "explicit" answer.
  # pkpd_ir.mod declares COMP=(EFFECT,DEFOBS) and is observed on EFFECT.
  p  <- nm_path("pkpd_ir.mod")
  ui <- nonmem2rx::nonmem2rx(p)
  ir <- suppressWarnings(rxui_to_ir(ui, source_format = "nonmem",
                                    scaling_hint = .extract_nm_scaling(p),
                                    obs_hint     = .extract_nm_defobs(p)))
  expect_equal(ir$structural$obs_cmt, "EFFECT")
  expect_length(grep("could not be inferred", ir$warnings), 0L)
})

test_that("DEFOBS decides obs_cmt AND the scaling compartment", {
  skip_if_not_installed("nonmem2rx")
  # The regression this exists for: with DEFOBS declared first and no explicit
  # compartment in $ERROR, the old tail(state_names, 1) guess picked PERIPH,
  # and because the scaling lookup derives its compartment number from the same
  # guess, the correctly parsed `S1 = V` was discarded with no diagnostic.
  p  <- nm_path("defobs_not_last.ctl")
  ui <- nonmem2rx::nonmem2rx(p)
  sc <- .extract_nm_scaling(p)
  expect_equal(sc[["1"]], "V")          # parsing was never the problem

  fixed <- suppressWarnings(rxui_to_ir(ui, source_format = "nonmem",
                                       scaling_hint = sc,
                                       obs_hint = .extract_nm_defobs(p)))
  expect_equal(fixed$structural$obs_cmt, "CENT")
  expect_equal(fixed$scaling$obs_scale, "V")
  expect_length(grep("obs_cmt could not be inferred", fixed$warnings), 0L)
})

test_that("a DEFOBS that disagrees with d/dt order is refused, not trusted", {
  skip_if_not_installed("nonmem2rx")
  # $MODEL says compartment 1 is NOPE; the first d/dt is for CENT. Believing the
  # index anyway would silently observe the wrong compartment. Uses a model
  # whose DV expression names no compartment, so the hint is what is on trial.
  p  <- nm_path("pkpd_ir.mod")
  ui <- nonmem2rx::nonmem2rx(p)
  ir <- suppressWarnings(rxui_to_ir(ui, source_format = "nonmem",
                                    obs_hint = list(index = 1L, name = "NOPE",
                                                    n_comp = 4L)))
  expect_match(ir$warnings, "orderings disagree", all = FALSE)
  expect_equal(ir$structural$obs_cmt, "EFFECT")   # fell back to the guess
})

test_that("a malformed obs_hint index is refused, not evaluated", {
  skip_if_not_installed("nonmem2rx")
  # rxui_to_ir() is exported and obs_hint is a documented @param, so a length-0,
  # NA or length-2 index must not reach `&&` and abort with a base-R condition
  # error naming neither field.
  ui <- nonmem2rx::nonmem2rx(nm_path("pkpd_ir.mod"))
  for (idx in list(NA_integer_, integer(0), c(1L, 2L), "1")) {
    ir <- expect_no_error(suppressWarnings(
      rxui_to_ir(ui, source_format = "nonmem",
                 obs_hint = list(index = idx, name = "EFFECT", n_comp = 4L))))
    expect_true(nzchar(ir$structural$obs_cmt))
  }
})

test_that(".same_cmt_name sees through nonmem2rx's c. prefix and case", {
  expect_true(.same_cmt_name("c.RTOT", "RTOT"))
  expect_true(.same_cmt_name("central", "CENTRAL"))
  expect_true(.same_cmt_name("A.B", "A_B"))
  expect_false(.same_cmt_name("CENT", "PERIPH"))
  expect_false(.same_cmt_name(NA_character_, "CENT"))
})

# -- statements with a call-shaped assignment target --------------------------

test_that("a dropped f()/alag() statement is reported, not silently discarded", {
  ini <- rbind(theta_row("t.TCL", 1), theta_row("t.TV", 10), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.TCL * exp(eta1)), quote(v <- t.TV),
              quote(f(depot) <- BIO.AV), quote(alag(depot) <- 0.5),
              quote(d/dt(depot) <- -cl * depot),
              quote(d/dt(central) <- cl * depot - cl/v * central))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_match(ir$warnings, "bioavailability for compartment 'depot'", all = FALSE)
  expect_match(ir$warnings, "dose lag time for compartment 'depot'",   all = FALSE)
})

test_that("a feature ferx supports is NOT filed as a ferx-core gap", {
  # ferx-core maps F{cmt} -> `f=`, ALAG{cmt} -> `lagtime=` and D/R -> duration
  # and rate, and ferx-r ships bioavailability.ferx and warfarin_ode_lagtime.ferx
  # as worked examples. So the translator not emitting these is a ferxtranslate
  # limitation, and $unsupported -- which CLAUDE.md defines as the ferx-core
  # prioritisation signal -- must stay empty for them. Filing them there asked
  # the engine team to build what they had already shipped.
  ini <- rbind(theta_row("t.TCL", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.TCL * exp(eta1)), quote(f(depot) <- 0.7),
              quote(alag(depot) <- 0.5), quote(dur(depot) <- 2),
              quote(d/dt(depot) <- -cl * depot))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_length(ir$unsupported, 0L)
  # ... but the user is still told, at ERROR level, that it was dropped.
  expect_match(ir$warnings, "supported by ferx but is not yet emitted", all = FALSE)
})

test_that("a state initial condition is translated, not reported as a gap", {
  # ferx parses `init(STATE) = <expr>` inside [odes]. The expression may
  # reference individual parameters, other states and literals.
  ini <- rbind(theta_row("t.TBL", 3), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(bl <- t.TBL * exp(eta1)),
              quote(EFFECT(0) <- bl),
              quote(d/dt(EFFECT) <- -bl * EFFECT))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_length(ir$initial_conditions, 1L)
  expect_equal(ir$initial_conditions[[1]]$state, "EFFECT")
  expect_equal(ir$initial_conditions[[1]]$rhs, "BL")
  expect_match(emit_ferx(ir), "init(EFFECT) = BL", fixed = TRUE)
  expect_length(ir$unsupported, 0L)
})

test_that("an init referencing a solver builtin ferx rejects is dropped", {
  # An init expression has a NARROWER scope than a d/dt right-hand side. Measured
  # against ferx 0.3.0: MACHEPS and TIME resolve there, T / TAFD / TAD do not,
  # and the engine says so -- "may only reference declared states (0 at init
  # time), individual parameters, or the MACHEPS constant". Treating the whole
  # of .RESERVED_ODE_NAMES as free let this through the scope guard and emitted
  # a file the engine then refused, which is the outcome the guard exists to
  # prevent.
  for (nm in c("TAD", "TAFD")) {
    ini <- rbind(theta_row("t.TK", 1), eta_row("eta1", 0.09, 1L))
    lst <- list(quote(k <- t.TK * exp(eta1)),
                bquote(EFFECT(0) <- .(as.name(nm))),
                quote(d/dt(EFFECT) <- -k * EFFECT))
    ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                      source_format = "nonmem"))
    expect_length(ir$initial_conditions, 0L)
    expect_match(ir$warnings, "does not resolve inside an init expression",
                 all = FALSE, info = nm)
  }
})

test_that("an init referencing MACHEPS is emitted", {
  # The other half: MACHEPS does resolve there, so narrowing the allowlist must
  # not have narrowed it past what ferx accepts.
  ini <- rbind(theta_row("t.TK", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- t.TK * exp(eta1)), quote(EFFECT(0) <- MACHEPS),
              quote(d/dt(EFFECT) <- -k * EFFECT))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_length(ir$initial_conditions, 1L)
})

test_that("TIME in a COMPOUND init is substituted, not grounds for dropping", {
  # The case the bare-TIME test cannot see, and the reason substitution beats
  # dropping. Model time at init is zero and the engine computes exactly that
  # (`init = TIME + 50` and `init = 0 + 50` agree to max |diff| = 0.0), so the
  # rest of the expression is still the value the source asked for. Dropping
  # started the compartment at 0 when the correct value was in hand.
  mk <- function(rhs) {
    ini <- rbind(theta_row("t.TK", 1), eta_row("eta1", 0.09, 1L))
    lst <- list(quote(k <- t.TK * exp(eta1)),
                bquote(EFFECT(0) <- .(rhs)),
                quote(d/dt(EFFECT) <- -k * EFFECT))
    suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                source_format = "nonmem"))
  }
  a <- mk(quote(TIME + 50))
  expect_length(a$initial_conditions, 1L)
  expect_equal(a$initial_conditions[[1]]$rhs, "0 + 50")
  expect_match(a$warnings, "replaced with 0", all = FALSE)

  b <- mk(quote(k + TIME))
  expect_length(b$initial_conditions, 1L)
  expect_equal(b$initial_conditions[[1]]$rhs, "K + 0")

  # An init with no TIME at all is untouched by any of this.
  d <- mk(quote(50))
  expect_equal(d$initial_conditions[[1]]$rhs, "50")
  expect_length(grep("replaced with 0", d$warnings), 0L)

  # The substitution is annotated AT THE LINE, not only in the header block.
  # `K + 0` reads as a translator bug to anyone who does not know why, and the
  # header warning is twenty lines away. Verified separately that ferx accepts a
  # comment inside [odes], so this cannot make a valid file invalid.
  # Anchored on the STATEMENT, not on the text: the header WARNING block quotes
  # `init(EFFECT) = K + 0` too, so a plain fixed match finds two lines and the
  # assertion silently tests the wrong one.
  odes <- strsplit(emit_ferx(b), "\n")[[1]]
  i <- grep("^\\s*init\\(EFFECT\\)", odes)
  expect_length(i, 1L)
  expect_match(odes[i - 1], "TIME replaced with 0")

  # ... and an untouched init carries no note.
  odes_d <- strsplit(emit_ferx(d), "\n")[[1]]
  j <- grep("^\\s*init\\(EFFECT\\)", odes_d)
  expect_length(j, 1L)
  expect_false(grepl("^\\s*#", odes_d[j - 1]))
})

test_that("an init referencing TIME is dropped even though ferx accepts it", {
  # The ONE place this deliberately does not mirror the engine. ferx accepts
  # `init(X) = TIME` -- by accident, since a bare TIME is an AST node rather
  # than a variable and the undefined-name check never sees it -- and then reads
  # it as 0, because model time at init is zero by definition. Measured by
  # forward simulation: `init = TIME + 50` and `init = 0 + 50` give identical
  # PRED, and `init = TIME` gives PRED 0 throughout. Emitting it would produce a
  # model that differs from the source with nothing to show for it.
  ini <- rbind(theta_row("t.TK", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- t.TK * exp(eta1)), quote(EFFECT(0) <- TIME),
              quote(d/dt(EFFECT) <- -k * EFFECT))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_length(ir$initial_conditions, 0L)
  # Dropped because substituting TIME := 0 reduces it to a bare 0, which is the
  # value every compartment already has -- the same no-op arm as `F1 = 1`. It is
  # NOT a special case for TIME; it falls out of the substitution.
  expect_match(ir$warnings, "dropping it does not change the model", all = FALSE)
})

test_that("an init referencing a theta is dropped with a scope explanation", {
  # Measured against the engine: `init(X) = TVBL` is an E_PARSE naming TVBL
  # undefined, because a theta is not in scope inside an init expression. Emit
  # it and the file does not parse, so it is dropped and said out loud.
  ini <- rbind(theta_row("t.TBL", 3), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- t.TBL * exp(eta1)),
              quote(EFFECT(0) <- t.TBL),
              quote(d/dt(EFFECT) <- -k * EFFECT))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_length(ir$initial_conditions, 0L)
  expect_match(ir$warnings, "does not resolve inside an init expression",
               all = FALSE)
})

test_that("an initial condition of 0 is a no-op, not a dropped feature", {
  ini <- rbind(theta_row("t.TK", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- t.TK * exp(eta1)), quote(EFFECT(0) <- 0),
              quote(d/dt(EFFECT) <- -k * EFFECT))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_length(ir$initial_conditions, 0L)
  expect_length(ir$unsupported, 0L)
  expect_match(ir$warnings, "already every compartment's initial value",
               all = FALSE)
})

test_that("an alias bound twice is not read as a constant", {
  # The relative-bioavailability idiom `F1 = 1` / `IF (FORM.EQ.2) F1 = THETA(4)`.
  # Taking the FIRST binding made const(1) true, so the translator said "sets
  # the value ferx already uses, so dropping it does not change the model" --
  # asserting nothing was lost while the formulation-dependent F disappeared.
  ini <- rbind(theta_row("t.TCL", 1), theta_row("t.T4", 0.7), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.TCL * exp(eta1)),
              quote(rxf.rxddta1. <- 1),
              quote(f(depot) <- rxf.rxddta1.),
              quote(if (FORM == 2) { rxf.rxddta1. <- t.T4; f(depot) <- rxf.rxddta1. }),
              quote(d/dt(depot) <- -cl * depot))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_length(grep("does not change the model", ir$warnings), 0L)
  expect_match(ir$warnings, "binds it more than once", all = FALSE)
})

test_that("a nonmem2rx alias temporary never reaches [individual_parameters]", {
  # `.resolve_alias()` knew rxdur/rxrate/rxalag and pass 3's skip list did not,
  # so the companion binding of a dropped statement was emitted as a model
  # parameter named after a nonmem2rx internal.
  ini <- rbind(theta_row("t.TK", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- t.TK * exp(eta1)),
              quote(rxdur.rxddta1. <- 2), quote(dur(depot) <- rxdur.rxddta1.),
              quote(rxalag.rxddta2. <- 0.5), quote(alag(depot) <- rxalag.rxddta2.),
              quote(d/dt(depot) <- -k * depot))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  lhs <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_length(grep("^RX", lhs), 0L)
})

test_that("a name only reachable through a dropped statement is not a covariate", {
  # BIO.AV appears solely on the RHS of `f(depot) <- BIO.AV`, which is never
  # emitted. Reporting it as an illegal covariate prescribed renaming a data
  # column that fixes nothing, and put a phantom entry in $unsupported -- which
  # per CLAUDE.md is the ferx-core prioritisation signal.
  ini <- rbind(theta_row("t.TCL", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.TCL * exp(eta1)), quote(f(depot) <- BIO.AV),
              quote(d/dt(depot) <- -cl * depot))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_length(grep("covariate reference", ir$warnings), 0L)
  expect_length(grep("BIO.AV", ir$unsupported, fixed = TRUE), 0L)
})

test_that("a genuinely illegal covariate is still reported", {
  # The other half: this one IS emitted, so the diagnostic is real and must
  # survive the narrowing above.
  ini <- rbind(theta_row("t.TCL", 1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(cl <- t.TCL * exp(eta1) * WT.KG),
              quote(d/dt(depot) <- -cl * depot))
  ir <- suppressWarnings(rxui_to_ir(list(iniDf = ini, lstExpr = lst),
                                    source_format = "nonmem"))
  expect_match(ir$warnings, "covariate reference", all = FALSE)
  expect_match(ir$unsupported, "WT.KG", fixed = TRUE, all = FALSE)
})

# -- theta pass-through into [odes] -------------------------------------------
#
# Issue #6 defect 2. A theta is not in scope in [odes]: ferx resolves a d/dt
# right-hand side against states, individual parameters and covariates only.
# Measured against ferx 0.2.0 and 0.3.0, a theta IS readable from
# [individual_parameters], `[scaling] y` and `obs_scale`, and is NOT readable
# from a d/dt right-hand side, an ODE-block intermediate, an `init()` expression
# or a pk macro argument. Only the second group needs a pass-through.

test_that("a theta referenced from an ODE gets a pass-through parameter", {
  # The shape with nothing to convert: the source names the theta straight in
  # $DES (`DADT(1) = -THETA(2)*A(1)`) and never assigns it in $PK, so there is no
  # alias for de-shadowing to turn into an individual parameter. This emitted a
  # bare theta into [odes] -- ferx cannot resolve it -- until the pass-through was
  # appended explicitly.
  ini <- rbind(theta_row("KTP", 0.5), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), ddt("CENT", quote(-k * CENT + KTP * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_equal(ir$odes[[1]]$rhs, "-K * CENT + KTP * CENT")
  ip <- vapply(ir$indiv_params, function(p) paste0(p$lhs, "=", p$rhs), "")
  expect_true("KTP=TVKTP" %in% ip)
  # The theta itself must be renamed, or the pass-through is the `KTP = KTP`
  # self-shadow that defect 14 is about.
  expect_true("TVKTP" %in% vapply(ir$thetas, function(t) t$name, ""))
  expect_match(ir$warnings, "is referenced from an ODE", all = FALSE)
})

test_that("a source-supplied alias is the pass-through, not a second one", {
  # `KTP = THETA(2)` in $PK arrives as `KTP <- KTP` once both sides normalise.
  # Pass 3 drops that as a self-assignment, but de-shadowed to `KTP <- TVKTP` it
  # survives and IS the pass-through -- the form ferx's own examples use. The
  # appended entry must stand down, or the parameter is defined twice.
  ini <- rbind(theta_row("KTP", 0.5), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(KTP <- KTP), quote(k <- K * exp(eta1)),
              ddt("CENT", quote(-k * CENT + KTP * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  lhs <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_equal(sum(lhs == "KTP"), 1L)
  expect_equal(ir$indiv_params[[match("KTP", lhs)]]$rhs, "TVKTP")
})

test_that("the ODE reference does not depend on where the alias is written", {
  # The same order dependency as the state renames, and the same cause: the ODE
  # right-hand side was resolved through `name_map`, which grows as the statements
  # are walked. `KTP <- KTP` above the d/dt line installed the alias in time and
  # the reference read the parameter; the identical line below it did not, and the
  # reference read the theta while the pass-through sat there unused.
  ini <- rbind(theta_row("KTP", 0.5), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  alias <- quote(KTP <- KTP)
  ka    <- quote(k <- K * exp(eta1))
  dd    <- ddt("CENT", quote(-k * CENT + KTP * CENT))

  before <- suppressWarnings(rxui_to_ir(mock_ui(ini, list(alias, ka, dd))))
  after  <- suppressWarnings(rxui_to_ir(mock_ui(ini, list(ka, dd, alias))))

  expect_equal(before$odes[[1]]$rhs, "-K * CENT + KTP * CENT")
  expect_equal(after$odes[[1]]$rhs,  before$odes[[1]]$rhs)
})

test_that("a theta read only from an individual parameter stays a theta", {
  # The negative case, and the one that keeps the pass-through from firing on
  # every theta in the model. KSS in the TMDD model is read from
  # [individual_parameters], where thetas ARE in scope, so renaming it and adding
  # a pass-through would be churn -- an extra parameter and a misleading shadow
  # INFO for a model that was already correct.
  ini <- rbind(theta_row("KSS", 2), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1) * KSS), ddt("CENT", quote(-k * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_true("KSS" %in% vapply(ir$thetas, function(t) t$name, ""))
  expect_false("KSS" %in% vapply(ir$indiv_params, function(p) p$lhs, ""))
  expect_false(any(grepl("KSS", ir$warnings, fixed = TRUE)))
})

test_that("a state is not mistaken for a theta reference in an ODE", {
  # A theta NAMED CENT alongside a state CENT. `-k * CENT` in [odes] is the
  # compartment, not the theta, so nothing here needs a carrier. The visible
  # symptom of getting it wrong is the rename: the theta was moved to TVCENT and
  # reported as shadowing an individual parameter that does not exist. (The ODE
  # itself survives, because discovery is recomputed each round and the renamed
  # theta no longer matches -- so this asserts the rename, which is what breaks.)
  ini <- rbind(theta_row("CENT", 1), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), ddt("CENT", quote(-k * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_true("CENT" %in% vapply(ir$thetas, function(t) t$name, ""))
  expect_equal(ir$odes[[1]]$rhs, "-K * CENT")
  expect_false("CENT" %in% vapply(ir$indiv_params, function(p) p$lhs, ""))
  expect_false(any(grepl("referenced from an ODE", ir$warnings)))
  expect_false(any(grepl("theta 'CENT' shares a name", ir$warnings)))
})

test_that("a theta reaching the ODE through an intermediate is scoped too", {
  # An intermediate that touches a state cannot be an individual parameter, so it
  # is emitted into [odes] -- and the text that gets emitted was normalised for a
  # different context, never against the [odes] scope. `ki <- KTP*CENT` therefore
  # arrived as `TVKTP * CENT`: a bare theta in [odes], with the pass-through
  # parameter defined and unreferenced beside it. This is why the scope is
  # applied to the emitted right-hand side rather than during the walk.
  #
  # Phase 5b changed HOW the intermediate reaches the block -- declared in source
  # order rather than substituted into the d/dt line -- but not what must hold:
  # no theta name may appear anywhere in [odes].
  ini <- rbind(theta_row("KTP", 0.5), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), quote(ki <- KTP * CENT),
              ddt("CENT", quote(-k * CENT + ki)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  ki <- Filter(function(o) identical(o$lhs, "KI"), ir$odes)
  expect_length(ki, 1L)
  expect_equal(ki[[1]]$rhs, "KTP * CENT")
  # The carrier, never the theta it carries.
  expect_false(any(grepl("TVKTP", .emitted_ode_symbols(ir$odes), fixed = TRUE)))
  expect_true("KTP" %in% vapply(ir$indiv_params, function(p) p$lhs, ""))
  # Declared before the line that reads it: [odes] has no use-before-def check,
  # so the order is the whole correctness property.
  kinds <- vapply(ir$odes, function(o) if (is.null(o$kind)) "ddt" else o$kind, "")
  expect_lt(which(kinds == "assign")[1], which(kinds == "ddt")[1])
})

test_that("an intermediate that survives as a parameter carries the theta", {
  # The other half of the same case: `ki <- KTP` has no state reference, so it
  # becomes the individual parameter `KI = KTP` and no pass-through is needed --
  # thetas are in scope in [individual_parameters]. What must hold either way is
  # that no theta name reaches the ODE, and here that is achieved by the model's
  # own intermediate rather than by anything this package adds.
  ini <- rbind(theta_row("KTP", 0.5), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), quote(ki <- KTP),
              ddt("CENT", quote(-k * CENT + ki * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_equal(ir$odes[[1]]$rhs, "-K * CENT + KI * CENT")
  expect_true("KI=KTP" %in%
              vapply(ir$indiv_params, function(p) paste0(p$lhs, "=", p$rhs), ""))
  # KTP is read from [individual_parameters] only, so it keeps its name.
  expect_true("KTP" %in% vapply(ir$thetas, function(t) t$name, ""))
})


# -- [odes] declaredness ------------------------------------------------------
#
# Every name the emitted [odes] block references must resolve. ferx reports this
# itself, but only where the engine runs -- and the phase-2 legality check beside
# it tests the GRAMMAR, not whether anything declares the name, so `KTP` passes it.
# Two defects reached the corpus through that gap (issue #6 defects 2 and 4).
#
# The set of names a ferx ODE RHS may reference is CLOSED: declared states,
# individual parameters, ODE-block intermediates and the reserved time variables.
# Not thetas, etas or sigmas -- and not covariates either, which is what makes the
# check possible without consulting a covariate list. That matters, because
# `.covariate_names()` defines a covariate as a symbol nothing binds and rxode2's
# `ui$allCovs` does much the same, so checking against either whitelists exactly the
# names worth reporting. An earlier version did, passed the whole suite, and could
# not fire.

test_that("an eta referenced from an ODE is reported as untranslatable", {
  # Not just a check test -- this is a real gap. An eta is out of scope in [odes]
  # (measured on ferx 0.2.0 and 0.3.0, same as thetas and sigmas), so ferx rejects
  # it as E_PARSE, and we used to emit it without comment. There is no carrier for
  # it either: an eta cannot be read from [individual_parameters] into an ODE
  # without restructuring the model, so this is reported rather than fixed.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ cl <- 1.0; v <- 10.0; prop.err <- 0.1; eta.cl ~ 0.09 })
    model({ k <- cl/v
            d/dt(central) = -k*central*exp(eta.cl)
            central ~ prop(prop.err) })
  }
  ir <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))

  expect_match(ir$warnings, "^ERROR \\| \\[odes\\] references ETA_CL", all = FALSE)
  expect_match(ir$unsupported, "theta/eta/sigma referenced from \\[odes\\]", all = FALSE)
})

test_that("a covariate referenced from an ODE is reported, not assumed valid", {
  # This started out as the opposite test -- a covariate in [odes] was assumed
  # legitimate and only warned about when no data-column list was available to
  # confirm it. ferx does not allow one at all: "An ODE RHS may only reference
  # declared states, individual parameters, ODE-block intermediates, or the reserved
  # TIME/TAFD/TAD/MACHEPS variables ... pre-compute the covariate-dependent term in
  # [individual_parameters]". So the set is closed, every leftover is an error, and
  # no column list is needed to decide.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ cl <- 1.0; v <- 10.0; prop.err <- 0.1; eta.cl ~ 0.09 })
    model({ cl <- cl*exp(eta.cl); k <- cl/v
            d/dt(central) = -k*central*(WT/70)
            central ~ prop(prop.err) })
  }
  ir <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))

  expect_match(ir$warnings, "^ERROR \\| \\[odes\\] references WT", all = FALSE)
  # The message has to carry the remedy: the term is expressible, just one block
  # earlier, and a user told only "undeclared" would not know that.
  expect_match(ir$warnings, "pre-computed in \\[individual_parameters\\]", all = FALSE)
  expect_match(ir$unsupported, "undeclared name referenced from \\[odes\\]", all = FALSE)
})

test_that("TIME and the reserved ODE variables are not reported", {
  # The closed set includes ferx's own time variables. Omitting them would make the
  # check fire on any model with a time-dependent rate -- the kind of false positive
  # that gets a check deleted rather than fixed.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ cl <- 1.0; v <- 10.0; prop.err <- 0.1; eta.cl ~ 0.09 })
    model({ cl <- cl*exp(eta.cl); k <- cl/v
            d/dt(central) = -k*central*exp(-TIME/24)
            central ~ prop(prop.err) })
  }
  ir <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))

  expect_false(any(grepl("\\[odes\\] references", ir$warnings)))
  expect_length(ir$unsupported, 0L)
})

test_that("a carrier leaves the [odes] block fully declared", {
  # The positive control for the check: the model that needs a carrier must come
  # out clean, or the check and the fix disagree and one of them is wrong.
  ini <- rbind(theta_row("KTP", 0.5), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), ddt("CENT", quote(-k * CENT + KTP * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_false(any(grepl("\\[odes\\] references", ir$warnings)))
  expect_length(ir$unsupported, 0L)
})

test_that("an undeclared name spelled like a math function is still reported", {
  # The check used to declare a list of function names on the theory that
  # `.emitted_ode_symbols()` returns call heads. It does not -- `.collect_symbols()`
  # recurses over `expr[-1]` -- so the list could only ever whitelist ordinary
  # identifiers that happen to share a spelling with a function, and it did:
  # this model reported nothing while the identical one with WT reported WT, and
  # ferx rejected the file with "RHS references undefined name(s): MAX".
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ cl <- 1.0; v <- 10.0; prop.err <- 0.1; eta.cl ~ 0.09 })
    model({ cl <- cl*exp(eta.cl); k <- cl/v
            d/dt(central) = -k*central*(MAX/70)
            central ~ prop(prop.err) })
  }
  ir <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))

  expect_match(ir$warnings, "^ERROR \\| \\[odes\\] references MAX", all = FALSE)
  expect_match(ir$unsupported, "undeclared name referenced from \\[odes\\]: MAX",
               all = FALSE)
})

test_that("function calls in an ODE are not reported as undeclared names", {
  # The other half of removing that whitelist: dropping it must not start
  # reporting the functions themselves. It cannot, because a call head is never
  # collected -- but the check is worth pinning, since the whole list was added
  # to prevent a false positive that could not happen.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ cl <- 1.0; v <- 10.0; prop.err <- 0.1; eta.cl ~ 0.09 })
    model({ cl <- cl*exp(eta.cl); k <- cl/v
            d/dt(central) = -k*central*exp(-TIME/24) - sqrt(abs(k))*log(1 + k)
            central ~ prop(prop.err) })
  }
  ir <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))

  expect_false(any(grepl("\\[odes\\] references", ir$warnings)))
  expect_length(ir$unsupported, 0L)
})


# -- carrier / [scaling] ordering ---------------------------------------------

test_that("[scaling] still resolves when the scaled variable needs a carrier", {
  # S2 = VC, with VC referenced from the ODE and never assigned in $PK. The
  # carrier moves the name that answers to `VC` from the theta list to the
  # individual parameters, so [scaling] must be resolved AFTER the carrier is
  # appended. Resolved before, the lookup found the theta renamed to TVVC and no
  # parameter named VC, `matched` came back NULL, and [scaling] was dropped with
  # no diagnostic -- an emitted model that predicts amounts against concentration
  # data and validates clean, which is exactly the S2=V failure CLAUDE.md warns
  # about with the loud half removed.
  ini <- rbind(theta_row("VC", 10), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), ddt("CENT", quote(-k * CENT/VC)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst),
                                     scaling_hint = list("1" = "VC")))

  expect_equal(ir$scaling$obs_scale, "VC")
  expect_match(ir$warnings, "S1 = VC detected", all = FALSE)
  # ...and it must be the carrier that answers, not the theta: a theta IS in
  # scope in obs_scale, so pointing [scaling] at TVVC would also validate, and
  # would silently drop the IIV-free/IIV distinction the carrier exists to keep.
  expect_true("VC=TVVC" %in%
              vapply(ir$indiv_params, function(p) paste0(p$lhs, "=", p$rhs), ""))
})


# -- carrier naming and the state exclusion, case-correctly -------------------

test_that("a carrier never takes the name of a ferx solver builtin", {
  # `.deshadow_theta_names()` renames a theta off a builtin name (theta TIME ->
  # TVTIME) because ferx resolves the bare name to the integrator's clock and the
  # theta would be estimated and never read. The carrier then took the theta's
  # SOURCE name as its first choice and put the collision straight back, one block
  # lower: `TIME = TVTIME` in [individual_parameters], with the ODE term reading
  # the clock. `builtin_params` reported it, but as a source defect -- "rename the
  # variable in the source model" -- for a variable the translator had invented.
  ini <- rbind(theta_row("TIME", 0.5), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), ddt("CENT", quote(-k * CENT + TIME * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  lhs <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_false(any(toupper(lhs) %in% .RESERVED_ODE_NAMES))
  # The carrier still exists and still carries the theta -- the fix is the name,
  # not the mechanism.
  expect_true("TVTIME_ODE=TVTIME" %in%
              vapply(ir$indiv_params, function(p) paste0(p$lhs, "=", p$rhs), ""))
  expect_equal(ir$odes[[1]]$rhs, "-K * CENT + TVTIME_ODE * CENT")
  # ...and the translator no longer reports a collision it created itself.
  expect_false(any(grepl("collides with a ferx solver builtin \\(", ir$warnings) &
                   grepl("^ERROR", ir$warnings)))
  expect_length(ir$unsupported, 0L)
})

test_that("a theta sharing a lowercase state's name is left alone", {
  # The state exclusion compares `theta_orig`, which is always `.norm()`ed to
  # uppercase, against the state names. Keyed on the RAW d/dt target it never
  # matched for an nlmixr2 source, where states are lowercase: theta CENTRAL was
  # renamed to TVCENTRAL and reported as shadowing an individual parameter that
  # does not exist. Nothing here references the theta from the ODE -- `-central *
  # KK` is the compartment -- so no rename and no carrier are warranted.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ CENTRAL <- 1.0; v <- 10.0; prop.err <- 0.1; eta.v ~ 0.09 })
    model({ v  <- v*exp(eta.v)
            kk <- CENTRAL/v
            d/dt(central) = -central*kk
            central ~ prop(prop.err) })
  }
  ir <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))

  expect_true("CENTRAL" %in% vapply(ir$thetas, function(t) t$name, ""))
  expect_false("CENTRAL" %in% vapply(ir$indiv_params, function(p) p$lhs, ""))
  expect_false(any(grepl("theta 'CENTRAL' shares a name", ir$warnings)))
  expect_false(any(grepl("theta 'CENTRAL' is referenced from an ODE", ir$warnings)))

  # What the spurious rename was accidentally masking, and why the exclusion
  # cannot ship on its own. `kk <- CENTRAL/v` is an honest read of the theta in
  # [individual_parameters] scope, but its RHS "references a state" under the
  # case-folded comparison, so the inliner drags the text into [odes] -- the one
  # block where CENTRAL means the compartment. Measured against ferx 0.3.0 the
  # result validates clean and reports `W_UNUSED_PARAM: theta 'CENTRAL' ... not
  # referenced in any model expression`: the engine read it as the state, the
  # theta went dead, and the term became the amount squared over V. Renaming the
  # theta hid that by moving it off the colliding spelling; nothing detected it.
  # It is now reported, so `strict = TRUE` aborts instead of shipping the model.
  expect_match(ir$warnings,
               "^ERROR \\| 'central' names both an ODE state and a model parameter",
               all = FALSE)
  expect_match(ir$unsupported, "state/parameter source-name collision: central",
               all = FALSE)
})

test_that("a duplicate $THETA label produces one carrier, not two", {
  # The discovery predicate is per-theta, but duplicate labels give two thetas the
  # same source name, and collapsing the result to names threw away which of them
  # matched: the carrier loop tested `theta_orig[i] %in% ode_theta` and let both
  # through. One ODE reference then defined two parameters -- the one the
  # reference resolves to, and a dead `KTP = TVKTP` that ferx reports as
  # `computed but never used` -- while the run simultaneously warned that ferx
  # resolves every reference to the first of the pair.
  ini <- rbind(theta_row("KTP", 0.5), theta_row("KTP", 0.7),
               theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), ddt("CENT", quote(-k * CENT + KTP * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  ip  <- vapply(ir$indiv_params, function(p) paste0(p$lhs, "=", p$rhs), "")
  expect_length(grep("is referenced from an ODE", ir$warnings), 1L)
  # Every emitted parameter is read by something -- no orphan beside the one the
  # ODE actually names.
  carriers <- grep("^[^=]+=(TV|THETA_)?KTP", ip, value = TRUE)
  expect_length(carriers, 1L)
  expect_true(grepl(sub("=.*$", "", carriers), ir$odes[[1]]$rhs, fixed = TRUE))
  # The duplicate itself is still reported; this changes the carrier, not that.
  expect_match(ir$warnings, "duplicate \\$THETA label", all = FALSE)
})

test_that("scoping the ODEs does not resurrect the dropped rhs_expr field", {
  # The inlining pass rebuilds each ode as `list(state, rhs)` and its only
  # consumer runs earlier, so `rhs_expr` is gone by then. Writing it back in
  # `.scope_odes_to_params()` made the field present on carrier models and absent
  # everywhere else -- a phase-5 reader would get a correct expression from some
  # models and NULL from the rest.
  ini <- rbind(theta_row("KTP", 0.5), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), ddt("CENT", quote(-k * CENT + KTP * CENT)))
  carried <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))
  plain   <- suppressWarnings(rxui_to_ir(mock_ui(
    rbind(theta_row("K", 0.1), eta_row("eta1", 0.09, 1L)),
    list(quote(k <- K * exp(eta1)), ddt("CENT", quote(-k * CENT))))))

  # Same shape either way -- that is the assertion, not which fields there are.
  expect_equal(names(carried$odes[[1]]), names(plain$odes[[1]]))
  expect_false("rhs_expr" %in% names(carried$odes[[1]]))
})

# -- accidental dose attributes (issue #17) -----------------------------------

test_that(".is_dose_attr_name matches ferx's grammar, including its edges", {
  # Mirrors DoseAttr::from_indexed_name (ferx-core src/types.rs): one of the
  # prefixes LAGTIME/ALAG/F/D/R followed by a pure digit string denoting a
  # compartment >= 1, plus the bare F/LAGTIME/ALAG that occupy RESERVED_PK_SLOTS.
  expect_true(all(.is_dose_attr_name(
    c("F1", "D1", "R1", "ALAG1", "LAGTIME1", "F10", "D23"))))
  expect_true(all(.is_dose_attr_name(c("F", "ALAG", "LAGTIME"))))
  # Case-insensitive, as the engine's own lowercase comparison is.
  expect_true(all(.is_dose_attr_name(c("f1", "alag2", "lagtime3", "f"))))
  # `cmt >= 1`, so an all-zero suffix is NOT an attribute but a leading zero is.
  expect_false(.is_dose_attr_name("F0"))
  expect_false(.is_dose_attr_name("F00"))
  expect_true(.is_dose_attr_name("F01"))
  # A non-numeric suffix is not an attribute -- this is what keeps every
  # ordinary pharmacometric name (and our own generated ones) out of the net.
  expect_false(any(.is_dose_attr_name(
    c("CL", "V", "KA", "F_BIO", "FRAC", "DUR", "RATE", "D_1", "F1_PAR",
      "TVF1", "TVF1_ODE", "RBASE", "DSC"))))
})

test_that("a dose-attribute parameter name is renamed and every reference follows", {
  # The reporter's model (translator#17): an ordinary elimination rate constant
  # that happens to be called F1. Emitted verbatim, ferx reads it as
  # bioavailability AND as the rate constant, and every prediction is off by
  # exactly its value -- with no error and no warning from the engine.
  ini <- rbind(theta_row("F1", 0.1), theta_row("V", 50), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(F1 <- F1 * exp(eta1)), ddt("CENT", quote(-F1 * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  ip <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_false(any(.is_dose_attr_name(ip)))
  # The rename is only correct if the ODE follows it. A renamed declaration
  # beside an unrewritten reference is an undefined name, not a fix.
  expect_true(any(grepl("F1_PAR", ip, fixed = TRUE)))
  expect_match(ir$odes[[1]]$rhs, "F1_PAR", fixed = TRUE)
  expect_false(grepl("\\bF1\\b", ir$odes[[1]]$rhs))
  expect_match(ir$warnings, "shape of a ferx dose attribute", all = FALSE)
})

test_that("the carrier declines a source name ferx reads as a dose attribute", {
  # Reachable without the source ever assigning the name: a theta labelled F1 and
  # referenced only from $DES is de-shadowed to TVF1 -- because the carrier about
  # to be created predicts an individual parameter called F1 -- which frees F1 for
  # the carrier to take. Measured before the guard: `theta TVF1` beside
  # `F1 = TVF1`, a bioavailability nothing in the source asked for.
  ini <- rbind(theta_row("F1", 0.6), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), ddt("CENT", quote(-k * CENT - F1 * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  ip <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_false(any(.is_dose_attr_name(ip)))
  expect_match(ir$odes[[1]]$rhs, "TVF1_ODE", fixed = TRUE)
  # One message, naming the parameter that is actually in the file. Renaming the
  # carrier after the fact instead would leave the carrier's own INFO line
  # describing `F1 = TVF1`, which no longer exists.
  expect_length(grep("dose attribute", ir$warnings), 0L)
})

test_that("a name that only looks like a dose attribute is left alone", {
  # The opposite failure: over-matching would churn ordinary names and, for
  # F0 specifically, contradict the engine -- `cmt >= 1` makes F0 an ordinary
  # parameter, so renaming it would be wrong rather than merely noisy.
  ini <- rbind(theta_row("F0", 0.1), theta_row("FRAC", 0.5),
               eta_row("eta1", 0.09, 1L))
  lst <- list(quote(F0 <- F0 * exp(eta1)), quote(FRAC <- FRAC),
              ddt("CENT", quote(-F0 * FRAC * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  ip <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_true("F0" %in% ip)
  expect_false(any(grepl("_PAR$", ip)))
  expect_length(grep("dose attribute", ir$warnings), 0L)
})

test_that("the replacement dodges a name the model already uses", {
  # `.free_name()`'s fallback, reached when F1_PAR is itself taken. The point is
  # that the two names stay distinct -- collapsing them would merge two
  # parameters into one and lose an estimate silently.
  taken <- c("F1_PAR", "CL")
  out   <- .deconflict_dose_attr_names(c("F1", "F1_PAR", "CL"), taken = taken)
  expect_equal(names(out$map), "F1")
  expect_false(out$map[["F1"]] %in% taken)
  expect_false(.is_dose_attr_name(out$map[["F1"]]))
})

test_that("another individual parameter reading the renamed one follows it", {
  # The declaration is only half of it: [individual_parameters] right-hand sides
  # are a reference site too, and one left unrewritten names a parameter that no
  # longer exists. Rewriting goes through the parse tree rather than gsub, so
  # `R1` does not also hit `R10` or `XR1`.
  ini <- rbind(theta_row("R1", 2), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(R1 <- R1 * exp(eta1)), quote(KTOT <- K * R1),
              ddt("CENT", quote(-KTOT * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  ip  <- vapply(ir$indiv_params, function(p) p$lhs, "")
  rhs <- vapply(ir$indiv_params, function(p) p$rhs, "")
  expect_false(any(.is_dose_attr_name(ip)))
  expect_match(rhs[ip == "KTOT"], "R1_PAR", fixed = TRUE)
  # Every name a right-hand side references must still be declared somewhere.
  expect_false(any(grepl("\\bR1\\b", rhs)))
})

test_that(".is_dose_attr_name answers for NA rather than aborting", {
  # `startsWith(NA, p)` is NA, which made `any(hit)` NA and the guard an abort --
  # a translation dying mid-run instead of reporting anything.
  expect_false(.is_dose_attr_name(NA_character_))
  expect_equal(.is_dose_attr_name(c("F1", NA, "CL")), c(TRUE, FALSE, FALSE))
})

test_that(".dose_attr_kind names what ferx would do for every prefix", {
  # The switch() falls through to the lag default for any prefix it does not
  # name, so a typo in the D or R case would describe a modelled-infusion
  # parameter as a lag time -- telling the user the wrong reason their model
  # was changed. Only the F{n} branch was covered before.
  expect_match(.dose_attr_kind("F1"),  "bioavailability for doses into compartment 1")
  expect_match(.dose_attr_kind("D1"),  "duration for compartment 1")
  expect_match(.dose_attr_kind("R1"),  "rate for compartment 1")
  expect_match(.dose_attr_kind("ALAG2"),    "lag time for compartment 2")
  expect_match(.dose_attr_kind("LAGTIME3"), "lag time for compartment 3")
  expect_match(.dose_attr_kind("F"),        "bioavailability applied to every dose")
  expect_match(.dose_attr_kind("ALAG"),     "lag time applied to every dose")
  expect_match(.dose_attr_kind("LAGTIME"),  "lag time applied to every dose")
  # ferx parses the suffix as an integer, so F01 is compartment 1. Echoing "01"
  # points the reader at a numbering that appears nowhere in their model.
  expect_match(.dose_attr_kind("F01"), "compartment 1$")
})

test_that("a bare F individual parameter is renamed too", {
  # The highest-consequence name in the class: ferx applies a bare `F` as
  # bioavailability to every dose, with no compartment index and no RATE column
  # needed. Every other end-to-end fixture here uses an indexed name.
  ini <- rbind(theta_row("F", 0.1), theta_row("K", 0.5), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(F <- F * exp(eta1)), ddt("CENT", quote(-K * CENT * F)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  ip <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_false(any(.is_dose_attr_name(ip)))
  expect_true("F_PAR" %in% ip)
  expect_match(ir$odes[[1]]$rhs, "F_PAR", fixed = TRUE)
  expect_match(ir$warnings, "bioavailability applied to every dose", all = FALSE)
})

test_that("an init() expression follows the rename", {
  # phase 4 added init() to [odes], and its expressions name individual
  # parameters. An unrewritten one is an init referencing a name nothing
  # declares -- the undefined-name failure this whole change exists to prevent,
  # in the newest block.
  ini <- rbind(theta_row("F1", 5), theta_row("K", 0.1), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(F1 <- F1 * exp(eta1)),
              as.call(list(as.name("<-"), as.call(list(as.name("CENT"), 0)),
                           as.name("F1"))),
              ddt("CENT", quote(-K * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_length(ir$initial_conditions, 1L)
  expect_equal(ir$initial_conditions[[1]]$rhs, "F1_PAR")
  expect_true("F1_PAR" %in% vapply(ir$indiv_params, function(p) p$lhs, ""))
})

test_that("a rename does not land on a covariate name", {
  # This is the third naming authority in the file; the other two reserve
  # foreign names. A rename onto a data column is silent -- the declared
  # parameter wins and the covariate becomes unreferenceable.
  ini <- rbind(theta_row("F1", 0.1), theta_row("K", 0.5), eta_row("eta1", 0.09, 1L))
  # F1_PAR is bound by nothing, so .covariate_names() classifies it as a covariate.
  lst <- list(quote(F1 <- F1 * exp(eta1)),
              ddt("CENT", quote(-K * CENT * F1 * F1_PAR)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  ip <- vapply(ir$indiv_params, function(p) p$lhs, "")
  expect_false(any(.is_dose_attr_name(ip)))
  expect_false("F1_PAR" %in% ip)
  # The covariate reference must survive untouched.
  expect_match(ir$odes[[1]]$rhs, "F1_PAR", fixed = TRUE)
})

# -- phase 5b: statements, conditionals and the block partition ---------------

# Statement helpers, so the expectations read as what the emitted block says.
ip_of  <- function(ir) vapply(ir$indiv_params, function(p)
  if (identical(p$kind, "if")) paste0("if (", p$cond, ")") else paste0(p$lhs, " = ", p$rhs), "")
ode_of <- function(ir) vapply(ir$odes, function(o)
  switch(if (is.null(o$kind)) "ddt" else o$kind,
         `if`   = paste0("if (", o$cond, ")"),
         assign = paste0(o$lhs, " = ", o$rhs),
         paste0("d/dt(", o$state, ")")), "")

test_that("a $DES conditional is emitted into [odes] with both arms", {
  # Defect 4. The conditional used to match no branch in the parse loop and be
  # discarded, leaving the name it defines undefined in the output.
  ini <- rbind(theta_row("K", 0.1), theta_row("KS", 2), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)),
              quote(ct <- CENT / KS),
              quote(if (ct < 0) { cf <- 0 } else { cf <- ct * 2 }),
              ddt("CENT", quote(-k * CENT * cf)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  o <- ode_of(ir)
  expect_true(any(grepl("^if \\(CT < 0\\)$", o)))
  cond <- Filter(function(x) identical(x$kind, "if"), ir$odes)[[1]]
  expect_equal(vapply(cond$then,  function(x) paste0(x$lhs, " = ", x$rhs), ""), "CF = 0")
  expect_equal(vapply(cond$else_, function(x) paste0(x$lhs, " = ", x$rhs), ""), "CF = CT * 2")
  # And the name it defines is no longer undeclared.
  expect_false("CF" %in% setdiff(.emitted_ode_symbols(ir$odes),
                                 toupper(c(.stmt_declared(ir$odes, "ddt", "assign"),
                                           .ode_states(ir$odes),
                                           .ip_names(ir$indiv_params)))))
})

test_that("a state-dependent variable lands in [odes], never [individual_parameters]", {
  # The half of defect 4 that produces a file which PARSES and is numerically
  # wrong: [individual_parameters] is evaluated once per subject, but a variable
  # reading a compartment amount has to be evaluated at every integration step.
  ini <- rbind(theta_row("K", 0.1), theta_row("KS", 2), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), quote(ct <- CENT / KS),
              quote(fb <- ct / (KS + ct)),
              ddt("CENT", quote(-k * CENT * fb)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_false("FB" %in% .ip_names(ir$indiv_params))
  expect_true("FB" %in% .stmt_declared(ir$odes, "ddt", "assign"))
})

test_that("an ODE intermediate is declared before the line that reads it", {
  # [odes] has no use-before-def check. An intermediate below its consumer stays
  # valid, reads a stale slot, and collapses the prediction to a constant with no
  # diagnostic -- so ordering is the correctness property, not presentation.
  ini <- rbind(theta_row("K", 0.1), theta_row("KS", 2), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), quote(ct <- CENT / KS),
              ddt("CENT", quote(-k * CENT * ct)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  o <- ode_of(ir)
  expect_lt(grep("^CT = ", o), grep("^d/dt\\(CENT\\)", o))
})

test_that("a $PK conditional is emitted into [individual_parameters]", {
  # Defect 8. Capturing conditionals without routing them here would be WORSE
  # than the original defect -- captured and then dropped, in neither block and
  # with no diagnostic, so the code looks like it handles them.
  ini <- rbind(theta_row("T1", 1), theta_row("T3", 2), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(tvcl <- T1),
              quote(if (SEX == 1) tvcl <- T1 * T3),
              quote(cl <- tvcl * exp(eta1)),
              ddt("CENT", quote(-cl * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  ip <- ip_of(ir)
  expect_true(any(grepl("^if \\(SEX == 1\\)$", ip)))
  # Source order: defined, conditionally overridden, then read.
  expect_lt(grep("^TVCL = ", ip)[1], grep("^if \\(SEX", ip))
  expect_lt(grep("^if \\(SEX", ip), grep("^CL = ", ip))
})

test_that("$ERROR indicator variables are dropped, not emitted dead", {
  # Defect 6. W1/W2 reference nothing and nothing reads them, so no reachability
  # rule evicts them -- they arrive through pass 3's default. ferx would report
  # them as computed but never used.
  ini <- rbind(theta_row("K", 0.1), eta_row("eta1", 0.09, 1L),
               sigma_row("eps1", 0.1), sigma_row("eps2", 0.2))
  lst <- list(quote(k <- K * exp(eta1)),
              quote(w1 <- 0), quote(w2 <- 0),
              ddt("CENT", quote(-k * CENT)),
              quote(y <- CENT * (1 + w1 * eps1 + w2 * eps2)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_false(any(c("W1", "W2") %in% .ip_names(ir$indiv_params)))
  expect_match(ir$warnings, "\\$ERROR scaffolding", all = FALSE)
})

test_that("a theta is de-shadowed against a name assigned only in a branch", {
  # The defect CLAUDE.md opens with, one nesting level down. `.deshadow_theta_names()`
  # is fed the individual-parameter names; a name assigned only inside an `if`
  # never reached that list, so a theta of the same name silently shadowed it --
  # the branch assignment would be dead, with no diagnostic from ferx.
  ini <- rbind(theta_row("CL", 1), theta_row("F2", 2), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(if (SEX == 1) CL <- CL * F2),
              quote(k <- CL * exp(eta1)),
              ddt("CENT", quote(-k * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  thetas <- vapply(ir$thetas, function(t) t$name, "")
  expect_false("CL" %in% thetas)          # renamed away from the branch target
  expect_true("CL" %in% .ip_names(ir$indiv_params))
})

test_that("an unused ODE intermediate is dropped rather than emitted", {
  # `pk_1cmt_oral.mod` has an unused `CP = A(2)/V` in $DES. Inlining discarded it
  # for free; emitting it would produce a `computed but never used` warning from
  # the engine and put a name in the file the model does not use.
  ini <- rbind(theta_row("K", 0.1), theta_row("V", 10), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(k <- K * exp(eta1)), quote(cp <- CENT / V),
              ddt("CENT", quote(-k * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_false("CP" %in% .stmt_declared(ir$odes, "ddt", "assign"))
  expect_false("CP" %in% .ip_names(ir$indiv_params))
})

test_that("the scaffolding drop stays aligned when a conditional is present", {
  # `.ip_names()` walks into branches, so a conditional contributes as many names
  # as it assigns. A name-indexed logical then misaligns with the STATEMENT list
  # and deletes the wrong entry. This fixture forces the mismatch: one
  # conditional assigning two names, so five names span four statements.
  ini <- rbind(theta_row("K", 0.1), theta_row("T2", 2), eta_row("eta1", 0.09, 1L),
               sigma_row("eps1", 0.1))
  lst <- list(quote(k <- K * exp(eta1)),
              quote(if (SEX == 1) { aa <- T2; bb <- T2 * 2 }),
              quote(w1 <- 0),
              ddt("CENT", quote(-k * CENT * aa * bb)),
              quote(y <- CENT * (1 + w1 * eps1)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  # The scaffolding goes, and nothing else does.
  expect_false("W1" %in% .ip_names(ir$indiv_params))
  expect_true(all(c("K", "AA", "BB") %in% .ip_names(ir$indiv_params)))
  # The conditional itself must survive intact -- both arms, both names.
  conds <- Filter(function(p) identical(p$kind, "if"), ir$indiv_params)
  expect_length(conds, 1L)
  expect_equal(sort(vapply(conds[[1]]$then, function(x) x$lhs, "")), c("AA", "BB"))
})

test_that("a parameter read only inside a conditional counts as used", {
  # The `used` set was built from `p$rhs` over the statement list, and a
  # conditional has no `rhs` -- so every name read inside a branch looked unused.
  # Combined with the scaffolding rule, which DELETES on the strength of that
  # set, a parameter referenced only in a branch could be removed and leave the
  # branch naming something undeclared.
  ini <- rbind(theta_row("K", 0.1), theta_row("T2", 2), eta_row("eta1", 0.09, 1L),
               sigma_row("eps1", 0.1))
  lst <- list(quote(k <- K * exp(eta1)),
              quote(w1 <- T2),                       # read ONLY by the branch...
              quote(if (SEX == 1) k <- k * w1),      # ...here
              ddt("CENT", quote(-k * CENT)),
              quote(y <- CENT * (1 + w1 * eps1)))    # ...and by $ERROR
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  # W1 is referenced by an emitted conditional, so it is NOT scaffolding.
  expect_true("W1" %in% .ip_names(ir$indiv_params))
  # And nothing the emitted block names is left undeclared.
  declared <- c(.ip_names(ir$indiv_params), .ode_states(ir$odes),
                .stmt_declared(ir$odes, "ddt", "assign"),
                vapply(ir$thetas, function(t) t$name, ""))
  expect_true("W1" %in% declared)
})

test_that("a dose-attribute rename reaches inside a conditional", {
  # Two review findings in one fixture. The #17 rename pass assumed every
  # [individual_parameters] entry has `lhs`/`rhs`, so a model carrying BOTH a
  # $PK conditional and a dose-attribute-shaped name aborted in
  # `dose_out$map[[NULL]]` -- measured, not hypothetical.
  #
  # And `.ip_names()` reports a name assigned both at the top level and in a
  # branch twice, so the deconflicter renamed it twice, consuming a second
  # candidate (`F1_PAR` then `F1_PAR_1`) and reporting it twice for one
  # parameter. One decision per name, not per occurrence.
  ini <- rbind(theta_row("T1", 0.1), theta_row("T2", 2), eta_row("eta1", 0.09, 1L))
  lst <- list(quote(f1 <- T1),
              quote(if (SEX == 1) f1 <- T1 * T2),
              quote(k <- f1 * exp(eta1)),
              ddt("CENT", quote(-k * CENT)))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_false(any(.is_dose_attr_name(.ip_names(ir$indiv_params))))
  expect_true("F1_PAR" %in% .ip_names(ir$indiv_params))   # not F1_PAR_1
  # The branch assignment and the downstream reference follow the same rename.
  cond <- Filter(function(p) identical(p$kind, "if"), ir$indiv_params)[[1]]
  expect_equal(vapply(cond$then, function(x) x$lhs, ""), "F1_PAR")
  k <- Filter(function(p) identical(p$lhs, "K"), ir$indiv_params)[[1]]
  expect_match(k$rhs, "F1_PAR", fixed = TRUE)
  expect_length(grep("shape of a ferx dose attribute", ir$warnings), 1L)
})

# -- Phase 6a: error model structure (issue #6 defects 10, 11, 6f) ------------

# One two-sigma model, built so the SOURCE order of the two epsilons is the
# opposite of ferx's argument order. That is the whole point: `combined(a, b)` is
# (proportional, additive) in ferx-core, so emitting traversal order transposes
# the two SDs and the fit converges on the wrong answer with no diagnostic.
err_ui <- function(y_rhs) {
  ini <- rbind(theta_row("t.CL", 1), eta_row("eta1", 0.09, 1L),
               sigma_row("eps1", 0.04), sigma_row("eps2", 0.09))
  mock_ui(ini, list(quote(cl <- t.CL * exp(eta1)),
                    ddt("central", quote(-cl * central)),
                    as.call(list(as.name("<-"), as.name("y"), y_rhs))))
}
err_of <- function(y_rhs) {
  ir <- suppressWarnings(rxui_to_ir(err_ui(y_rhs)))
  if (!length(ir$error_model)) return(NA_character_)
  paste0(ir$error_model[[1]]$type, "(",
         paste(ir$error_model[[1]]$params, collapse = ", "), ")")
}

test_that("combined() emits ferx's (proportional, additive) order, not source order", {
  # Defect 10. The additive term is written FIRST here, so traversal order gives
  # combined(EPS1, EPS2) and the correct answer is combined(EPS2, EPS1). A source
  # that happened to write the proportional term first could not tell the two
  # apart -- the fixture only discriminates because the order is reversed.
  expect_equal(err_of(quote(central + eps1 + central * eps2)),
               "combined(EPS2, EPS1)")
  # ... and the mirror image, so the rule is not "always swap".
  expect_equal(err_of(quote(central * (1 + eps1) + eps2)),
               "combined(EPS1, EPS2)")
})

test_that("Y = F + F*EPS is proportional, not additive", {
  # Defect 11, and the common case: this is how NONMEM models usually spell
  # proportional error. The old rule keyed on the top-level call being `+`.
  expect_equal(err_of(quote(central + central * eps1)), "proportional(EPS1)")
  # The two spellings of the same model must agree.
  expect_equal(err_of(quote(central * (1 + eps1))), "proportional(EPS1)")
  # ... while a genuinely additive error is still additive.
  expect_equal(err_of(quote(central + eps1)), "additive(EPS1)")
})

test_that("a prediction that is an expression is still recognised as proportional", {
  # `CENT/VC` is the prediction, but `VC` is an ordinary individual parameter
  # that references no state. Any rule that decides "is this the prediction?"
  # from a whitelist of names fails here while passing every test above it.
  ini <- rbind(theta_row("t.VC", 3), eta_row("eta1", 0.09, 1L),
               sigma_row("eps1", 0.04))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, list(
    quote(vc <- t.VC * exp(eta1)),
    ddt("central", quote(-central/vc)),
    quote(y <- central/vc * (1 + eps1))))))
  expect_length(ir$error_model, 1L)
  expect_equal(ir$error_model[[1]]$type, "proportional")
})

test_that("an error expression ferx cannot express emits no [error_model]", {
  # Each of these is a real model that the counting classifier turned into a
  # confident, wrong answer. `NA` here means no [error_model] is emitted, which
  # makes the engine reject the file rather than fit it.
  expect_true(is.na(err_of(quote(central * (1 + w1 * eps1 + w2 * eps2)))))  # multi-endpoint
  expect_true(is.na(err_of(quote(central * exp(eps1)))))                    # non-linear
  expect_true(is.na(err_of(quote(central + 2 * eps1))))                     # scaled sigma
  expect_true(is.na(err_of(quote(central + eps1 + eps2))))                  # two additive
})

test_that("an untranslatable error expression is reported and carries a suggestion", {
  ir <- suppressWarnings(rxui_to_ir(err_ui(quote(central * (1 + w1 * eps1 + w2 * eps2)))))
  expect_length(ir$error_model, 0L)
  expect_match(ir$warnings, "^ERROR \\| could not determine the error model",
               all = FALSE)
  # The reason names the offending coefficient, not just "complex $ERROR".
  expect_match(ir$warnings, "weighted by", all = FALSE)

  txt <- emit_ferx(ir)
  expect_no_match(txt, "\n\\[error_model\\]")
  expect_match(txt, "# ferxtranslate could not translate this \\$ERROR expression")
  expect_match(txt, "# A plausible reading, NOT a translation")
  # Every suggestion line must be a comment. One uncommented line would make the
  # engine parse a guess as though it were the translation.
  sug <- grep("^#", strsplit(txt, "\n")[[1]], value = TRUE, invert = TRUE)
  expect_length(grep("DV ~", sug), 0L)
})

test_that("a translatable error expression produces no suggestion", {
  ir <- suppressWarnings(rxui_to_ir(err_ui(quote(central * (1 + eps1)))))
  expect_length(ir$error_model, 1L)
  expect_length(ir$error_suggestion, 0L)
  expect_no_match(emit_ferx(ir), "could not translate")
})

test_that("an unrecognised tilde error form is reported, not guessed (defect 6f)", {
  # Was `WARN | complex $ERROR -- classified as proportional, verify`, which is
  # the failure mode this issue is about: a guess that reads as a translation.
  # `lnorm()` is not one of the three forms .parse_error_rhs() understands.
  map <- .norm_map_from_ini(sigma_row("err.x", 0.01))
  out <- .parse_error_rhs(quote(lnorm(err.x)), map)
  expect_true(is.na(out$type))
  expect_length(out$params, 0L)
  expect_match(out$warnings, "^ERROR \\| could not determine the error model",
               all = FALSE)
})

test_that("an unrecognised tilde error form emits no [error_model] block", {
  # The whole point of returning NA: the file must NOT carry a guessed model.
  # Reached through the pipeline, since the abort lever is the engine rejecting
  # a file with no [error_model].
  ini <- rbind(theta_row("t.CL", 1), eta_row("eta1", 0.09, 1L),
               sigma_row("err.x", 0.01))
  ir  <- suppressWarnings(rxui_to_ir(mock_ui(
    ini, list(quote(cl <- t.CL * exp(eta1)),
              as.call(list(as.name("~"), as.name("DV"),
                           quote(lnorm(err.x))))))))
  expect_length(ir$error_model, 0L)
  txt <- emit_ferx(ir)
  expect_no_match(txt, "\n\\[error_model\\]")
  expect_match(txt, "# ferxtranslate could not translate")
})

test_that("sigma declarations are reordered to the roles the error model needs", {
  sig <- function(...) lapply(c(...), function(n) list(name = n, value = 1, scale = "sd"))
  nm  <- function(x) vapply(x, function(s) s$name, "")

  # ferx binds single-endpoint sigmas positionally from the declaration order and
  # discards the names, so the declaration list IS the role assignment.
  em <- list(list(dv = "DV", type = "combined", params = c("EPS2", "EPS1")))
  expect_equal(nm(.order_sigmas_for_error(sig("EPS1", "EPS2"), em)),
               c("EPS2", "EPS1"))

  # A sigma the error model does not reference keeps its place behind the ones it
  # does -- dropping it would change the model's parameter count.
  em1 <- list(list(dv = "DV", type = "proportional", params = "EPS2"))
  expect_equal(nm(.order_sigmas_for_error(sig("EPS1", "EPS2", "EPS3"), em1)),
               c("EPS2", "EPS1", "EPS3"))

  # Already in role order: unchanged, so an ordinary model is byte-identical.
  em2 <- list(list(dv = "DV", type = "combined", params = c("EPS1", "EPS2")))
  expect_equal(nm(.order_sigmas_for_error(sig("EPS1", "EPS2"), em2)),
               c("EPS1", "EPS2"))

  # Nothing to do, and nothing dropped, when there is no error model at all.
  expect_equal(nm(.order_sigmas_for_error(sig("EPS1", "EPS2"), list())),
               c("EPS1", "EPS2"))
  expect_length(.order_sigmas_for_error(list(), list()), 0L)

  # A multi-entry error model is left alone. Unreachable from rxui_to_ir() today
  # -- at most one entry is ever built -- so only a hand-built list can test it,
  # and it goes untested otherwise: relaxing the guard to `length(...) == 0`
  # leaves the whole suite green.
  #
  # It matters for phase 6b, which introduces per-CMT and covariate-selected
  # error models. ferx resolves THOSE by name, so reordering is unnecessary, and
  # the entries disagree about roles, so reordering would silently follow
  # whichever endpoint happened to be first.
  per_cmt <- list(list(dv = "DV", cmt = 1L, type = "proportional", params = "EPS2"),
                  list(dv = "DV", cmt = 3L, type = "additive",     params = "EPS1"))
  expect_equal(nm(.order_sigmas_for_error(sig("EPS1", "EPS2"), per_cmt)),
               c("EPS1", "EPS2"))
})

# -- Phase 6b: endpoint dispatch (issue #6 defects 5, 12, 15) -----------------

test_that("constant folding removes exactly the identities substitution creates", {
  f <- function(s) deparse1(.fold_consts(str2lang(s)))
  # The shape that matters: an indicator substituted to 1 and 0 leaves
  # `X * (1 + 1*0 + 0*0)`, and only folding turns that back into the readout.
  expect_equal(f("X * (1 + 1*0 + 0*0)"), "X")
  expect_equal(f("A/B * (1 + 1*EPS1 + 0*EPS2)"), "A/B * (1 + EPS1)")
  expect_equal(f("x + 0 * z"), "x")
  expect_equal(f("x^1"), "x")
  expect_equal(f("2*3 + x"), "6 + x")
  # Precedence survives dropping the `(` node, because deparse() re-parenthesises
  # from the tree. If it did not, `(a + b) * c` would emit as `a + b * c`.
  expect_equal(f("(a + b) * c"), "(a + b) * c")
  # NOT an algebra system: nothing here may rewrite terms that substitution did
  # not create. Widening a $ERROR heuristic is how this issue happened.
  expect_equal(f("x - x"), "x - x")
  expect_equal(f("2 * x + 3 * x"), "2 * x + 3 * x")
})

test_that("inlining resolves a chain and reports a cycle instead of truncating it", {
  defs <- list(IPRED = quote(CTOT), CTOT = quote(CENT/VC))
  expect_equal(deparse1(.inline_defs(quote(IPRED * (1 + EPS1)), defs)),
               "CENT/VC * (1 + EPS1)")
  # A name with no definition is left alone -- that is how a covariate survives.
  expect_equal(deparse1(.inline_defs(quote(FLAG + 1), defs)), "FLAG + 1")
  # The predecessor, .inline_aux_vars(), gave up at depth 30 and returned the
  # PARTLY inlined expression, so a cycle reached the file as an undefined name
  # with no diagnostic. NULL makes the caller report it.
  expect_null(.inline_defs(quote(A), list(A = quote(B), B = quote(A))))
})

# A captured conditional as pass 2d builds it: the condition already normalised.
cnd <- function(txt, ...) {
  arms <- list(...)
  list(kind = "if", cond = str2lang(txt), pos = 0L,
       then = lapply(arms, function(a) list(lhs = a[[1]], rhs = str2lang(a[[2]]))),
       else_ = NULL)
}

test_that("only a single-column equality dispatch is read as one", {
  ok <- .dispatch_conditions(list(cnd("FLAG == 2", c("IPRED", "RTOT")),
                                  cnd("FLAG == 1", c("W1", "1"))))
  expect_equal(ok$col, "FLAG")
  expect_equal(unlist(ok$values), c(2, 1))

  # No conditionals at all is single-endpoint, not a failure: `why` is NULL and
  # the caller falls through to the pre-6b path unchanged.
  expect_null(.dispatch_conditions(list())$why)

  # Two columns: the cases are no longer mutually exclusive, so the enumeration
  # would emit a model for a combination it never evaluated.
  expect_match(.dispatch_conditions(list(cnd("FLAG == 2", c("A", "1")),
                                         cnd("CMT == 2",  c("B", "1"))))$why,
               "more than one column")
  # An inequality has the same problem: `FLAG >= 1` and `FLAG >= 2` overlap.
  expect_match(.dispatch_conditions(list(cnd("FLAG >= 1", c("A", "1"))))$why,
               "not a `<column> == <number>` test")
  nested <- cnd("FLAG == 2", c("A", "1"))
  nested$then <- list(cnd("FLAG == 1", c("B", "1")))
  expect_match(.dispatch_conditions(list(nested))$why, "nested")
})

# The reporter's $ERROR block, as pass 2d hands it to .build_endpoints().
tmdd_chain <- function() list(
  list(kind = "assign", lhs = "CTOT",  rhs = quote(CENT/VC), pos = 1L),
  list(kind = "assign", lhs = "RTOT",  rhs = quote(c_RTOT),  pos = 2L),
  list(kind = "assign", lhs = "IPRED", rhs = quote(CTOT),    pos = 3L),
  cnd("FLAG == 2", c("IPRED", "RTOT")),
  list(kind = "assign", lhs = "W1", rhs = 0, pos = 5L),
  list(kind = "assign", lhs = "W2", rhs = 0, pos = 6L),
  cnd("FLAG == 1", c("W1", "1")),
  cnd("FLAG == 2", c("W2", "1")))

tmdd_aux <- c("CENT", "TISS", "C_RTOT", "CTOT", "RTOT", "IPRED", "EPS1", "EPS2")

test_that("an indicator-weighted Y resolves to one error model per endpoint", {
  ep <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1 + W2 * EPS2)),
                         tmdd_chain(), c("EPS1", "EPS2"),
                         tmdd_aux, c("CENT", "TISS", "c_RTOT"))
  expect_equal(ep$col, "FLAG")
  expect_length(ep$cases, 3L)          # FLAG == 2, FLAG == 1, fall-through

  expect_equal(ep$cases[[1]]$value, 2)
  expect_equal(deparse1(ep$cases[[1]]$pred), "c_RTOT")
  expect_equal(ep$cases[[1]]$type, "proportional")
  expect_equal(ep$cases[[1]]$params, "EPS2")

  expect_equal(ep$cases[[2]]$value, 1)
  expect_equal(deparse1(ep$cases[[2]]$pred), "CENT/VC")
  expect_equal(ep$cases[[2]]$params, "EPS1")

  # With every indicator 0 the source computes Y = IPRED: an observation with no
  # residual error. That is the source declaring its own dispatch exhaustive,
  # and it is why the last listed case has to become the `else`.
  expect_true(is.na(ep$cases[[3]]$type))
  expect_equal(deparse1(ep$cases[[3]]$pred), "CENT/VC")
})

test_that("an [odes] intermediate in the readout is rejected, not emitted", {
  # Measured on ferx 0.3.0: `y = CT` where CT is an [odes] intermediate
  # validates as VALID with no data and dies with E_MISSING_COVARIATE the moment
  # a data file is present -- ferx resolves the name inside [odes] and reads it
  # as a covariate everywhere else. So the check has to be ours.
  chain <- list(list(kind = "assign", lhs = "IPRED", rhs = quote(CT), pos = 1L),
                cnd("FLAG == 1", c("W1", "1")))
  ep <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1)), chain,
                         "EPS1", c(tmdd_aux, "CT"), c("CENT", "TISS", "c_RTOT"))
  expect_match(ep$why, "CT")
  expect_match(ep$why, "resolves inside [odes] but not in [scaling] y", fixed = TRUE)
})

test_that("the last dispatched case becomes the else when the fall-through has none", {
  ep <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1 + W2 * EPS2)),
                         tmdd_chain(), c("EPS1", "EPS2"),
                         tmdd_aux, c("CENT", "TISS", "c_RTOT"))
  out <- .assemble_endpoints(ep, c("CENT", "TISS", "c_RTOT"))
  expect_equal(out$readout$kind, "expr")
  expect_equal(out$readout$y, "if (FLAG == 2) c_RTOT else CENT/VC")
  expect_length(out$error_model, 2L)
  expect_equal(out$error_model[[1]]$cond, "FLAG == 2")
  expect_equal(out$error_model[[1]]$params, "EPS2")
  # ferx REQUIRES a terminating bare `else`; a NULL cond is how the emitter
  # knows which entry that is.
  expect_null(out$error_model[[2]]$cond)
  expect_equal(out$error_model[[2]]$params, "EPS1")
  expect_length(grep("exhaustive", out$warnings), 1L)
})

test_that("a fall-through that IS an error model is kept as the else", {
  # `W1 = 1` by default, cleared for FLAG == 2. Nothing has to be assumed
  # exhaustive here, so no assumption is announced.
  chain <- list(
    list(kind = "assign", lhs = "IPRED", rhs = quote(CENT/VC), pos = 1L),
    list(kind = "assign", lhs = "W1", rhs = 1, pos = 2L),
    list(kind = "assign", lhs = "W2", rhs = 0, pos = 3L),
    cnd("FLAG == 2", c("W1", "0"), c("W2", "1")))
  ep <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1) + W2 * EPS2), chain,
                         c("EPS1", "EPS2"), tmdd_aux, c("CENT", "TISS", "c_RTOT"))
  out <- .assemble_endpoints(ep, c("CENT", "TISS", "c_RTOT"))
  expect_equal(out$error_model[[1]]$cond, "FLAG == 2")
  expect_equal(out$error_model[[1]]$type, "additive")
  expect_null(out$error_model[[2]]$cond)
  expect_equal(out$error_model[[2]]$type, "proportional")
  expect_length(out$warnings, 0L)
  # One readout for both endpoints, so no branching `y`.
  expect_equal(out$readout$y, "CENT/VC")
})

test_that("a CMT dispatch emits y[CMT=N] and CMT=N: bodies, ordered by compartment", {
  # ferx does not expose CMT as a covariate -- measured, `Model references
  # covariate(s) not found in data (case-sensitive): CMT` -- so Form C cannot
  # dispatch on it and this path is required, not a stylistic alternative.
  chain <- tmdd_chain()
  chain[[4]]$cond <- quote(CMT == 3)
  chain[[7]]$cond <- quote(CMT == 1)
  chain[[8]]$cond <- quote(CMT == 3)
  ep <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1 + W2 * EPS2)), chain,
                         c("EPS1", "EPS2"), tmdd_aux, c("CENT", "TISS", "c_RTOT"))
  out <- .assemble_endpoints(ep, c("CENT", "TISS", "c_RTOT"))
  expect_equal(out$readout$kind, "per_cmt")
  expect_equal(vapply(out$readout$entries, function(e) e$cmt, 0L), c(1L, 3L))
  expect_equal(vapply(out$readout$entries, function(e) e$expr, ""),
               c("CENT/VC", "c_RTOT"))
  expect_equal(vapply(out$error_model, function(e) e$cmt, 0L), c(1L, 3L))
  expect_equal(vapply(out$error_model, function(e) e$params, ""), c("EPS1", "EPS2"))
  # No entry is invented for compartment 2, and the omission is announced --
  # ferx rejects the fit naming it if the data observes it.
  expect_length(grep("compartment\\(s\\) 2", out$warnings), 1L)
})

test_that("a CMT fall-through error model is expanded over the compartments no condition named", {
  # `IPRED = CONC` by default, `RESP` for CMT 2, one error model throughout.
  # The default really does apply to every other compartment, so naming them is
  # exact rather than a guess -- and it has to be done, because y[CMT=N] has no
  # `else` and ferx requires an entry for every observed CMT.
  chain <- list(
    list(kind = "assign", lhs = "IPRED", rhs = quote(CENT/VC), pos = 1L),
    cnd("CMT == 2", c("IPRED", "PD")))
  ep <- .build_endpoints(quote(IPRED * (1 + EPS1)), chain, "EPS1",
                         c("CENT", "PD", "IPRED", "EPS1"), c("CENT", "PD"))
  out <- .assemble_endpoints(ep, c("CENT", "PD"))
  expect_equal(vapply(out$readout$entries, function(e) e$cmt, 0L), c(1L, 2L))
  expect_equal(vapply(out$readout$entries, function(e) e$expr, ""),
               c("CENT/VC", "PD"))
  # One error model for both, so it stays single-endpoint rather than being
  # split into identical CMT=N: lines.
  expect_length(out$error_model, 1L)
  expect_null(out$error_model[[1]]$cmt)
  expect_length(out$warnings, 0L)
})

test_that("a single dispatched value with no fall-through model has no else to emit", {
  # ferx requires a terminating bare `else`, and there is only one branch to
  # give it. Reported rather than emitted as an unconditional model, which would
  # apply EPS1 to observations the source excluded from it.
  chain <- list(
    list(kind = "assign", lhs = "IPRED", rhs = quote(CENT/VC), pos = 1L),
    list(kind = "assign", lhs = "W1", rhs = 0, pos = 2L),
    cnd("FLAG == 1", c("W1", "1")))
  ep <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1)), chain, "EPS1",
                         tmdd_aux, c("CENT", "TISS", "c_RTOT"))
  out <- .assemble_endpoints(ep, c("CENT", "TISS", "c_RTOT"))
  expect_match(out$why, "no branch could serve as the `else`", fixed = TRUE)
})

test_that("an epsilon with a zero coefficient is dropped from that endpoint", {
  # Normally folding removes these first: substituting `W2 = 0` leaves
  # `0 * EPS2`, which collapses to `0` and takes EPS2 out of the expression
  # entirely. This is the shape that survives folding -- `.fold_consts()`
  # deliberately has no `x - x` rule -- so it is the only way to reach the
  # guard. Without it EPS2 reaches the classifier with coefficient 0, which is
  # neither 1 nor the prediction, and the whole endpoint is reported
  # untranslatable over a term the source had already cancelled.
  chain <- list(
    list(kind = "assign", lhs = "IPRED", rhs = quote(CENT/VC), pos = 1L),
    list(kind = "assign", lhs = "W1", rhs = 0, pos = 2L),
    cnd("FLAG == 1", c("W1", "1")))
  ep <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1) + EPS2 - EPS2), chain,
                         c("EPS1", "EPS2"), tmdd_aux, c("CENT", "TISS", "c_RTOT"))
  expect_null(ep$why)
  expect_equal(ep$cases[[1]]$type, "proportional")
  expect_equal(ep$cases[[1]]$params, "EPS1")
  expect_equal(deparse1(ep$cases[[1]]$pred), "CENT/VC")
})

test_that("the classifier reports its reason as a field, and names an eps-free prediction", {
  # `reason` exists so phase 6b can place the explanation in a message naming the
  # ENDPOINT. Deriving it by stripping the prefix off `warnings` would be string
  # surgery on a sentence.
  r <- .classify_error_assignment(quote(F + EPS1 * THETA4), "EPS1")
  expect_true(is.na(r$type))
  expect_match(r$reason, "`EPS1` is weighted by `THETA4`", fixed = TRUE)
  expect_match(r$warnings, r$reason, fixed = TRUE)
  # The prediction SHOWN is eps-free. `pred` is the raw expression -- zeroing
  # happens inside the evaluator -- so deparsing it printed the epsilon terms
  # back and the message read "neither 1 nor the prediction `F + EPS1 * THETA4`".
  expect_match(r$reason, "the prediction `F`", fixed = TRUE)
  expect_no_match(r$reason, "the prediction `F + EPS1", fixed = TRUE)
})

test_that("one unexpressible endpoint yields a readout plus a marked suggestion", {
  # A clean FLAG dispatch where only the FLAG == 2 endpoint carries a scaled
  # sigma. The readout is derived per case and independently of the error model,
  # so it is complete here -- emitting it leaves only the branch a human has to
  # write. Bailing on the whole dispatch instead sent the model down the
  # single-endpoint path, which re-diagnosed the UN-substituted `Y` and blamed
  # EPS1, the epsilon that was fine.
  chain <- list(
    list(kind = "assign", lhs = "CONC",  rhs = quote(CENT/VC), pos = 1L),
    list(kind = "assign", lhs = "PERIF", rhs = quote(PERI),    pos = 2L),
    list(kind = "assign", lhs = "IPRED", rhs = quote(CONC),    pos = 3L),
    cnd("FLAG == 2", c("IPRED", "PERIF")),
    list(kind = "assign", lhs = "W1", rhs = 0, pos = 5L),
    list(kind = "assign", lhs = "W2", rhs = 0, pos = 6L),
    cnd("FLAG == 1", c("W1", "1")),
    cnd("FLAG == 2", c("W2", "1")))
  aux <- c("CENT", "PERI", "CONC", "PERIF", "IPRED", "EPS1", "EPS2")
  ep <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1) + W2 * EPS2 * THETA3),
                         chain, c("EPS1", "EPS2"), aux, c("CENT", "PERI"))
  # FLAG == 2 first in source order, so it is cases[[1]].
  expect_true(is.na(ep$cases[[1]]$type))
  expect_match(ep$cases[[1]]$why, "the endpoint for FLAG == 2", fixed = TRUE)
  expect_match(ep$cases[[1]]$why, "`EPS2` is weighted by `THETA3`", fixed = TRUE)
  expect_equal(ep$cases[[2]]$type, "proportional")

  out <- .assemble_endpoints(ep, c("CENT", "PERI"))
  # The readout is emitted for real, and it is the complete one.
  expect_equal(out$readout$y, "if (FLAG == 2) PERI else CENT/VC")
  expect_null(out$error_model)
  expect_match(out$why, "FLAG == 2", fixed = TRUE)
  expect_true(all(grepl("^\\s*#", out$suggestion)))
  expect_true(any(grepl("if (FLAG == 2) { DV ~ ??? }", out$suggestion, fixed = TRUE)))
  expect_true(any(grepl("# <- Y = PERI + EPS2 * THETA3", out$suggestion, fixed = TRUE)))
  expect_true(any(grepl("else { DV ~ proportional(EPS1) }", out$suggestion, fixed = TRUE)))
})

test_that("a failed endpoint is never folded into the else ferx requires", {
  # The `else` takes no condition, so a gap placed there could not say which
  # endpoint it belonged to. The last listed case is the one that gets folded in
  # when the fall-through has no model, so it is the one that must have worked.
  chain <- list(
    list(kind = "assign", lhs = "IPRED", rhs = quote(CENT/VC), pos = 1L),
    list(kind = "assign", lhs = "W1", rhs = 0, pos = 2L),
    list(kind = "assign", lhs = "W2", rhs = 0, pos = 3L),
    cnd("FLAG == 1", c("W1", "1")),
    cnd("FLAG == 2", c("W2", "1")))
  aux <- c("CENT", "IPRED", "EPS1", "EPS2")
  # EPS2 (FLAG == 2, the LAST listed case) is the scaled one.
  ep <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1) + W2 * EPS2 * THETA3),
                         chain, c("EPS1", "EPS2"), aux, "CENT")
  out <- .assemble_endpoints(ep, "CENT")
  expect_null(out$suggestion)
  expect_null(out$readout)
  expect_match(out$why, "the endpoint for FLAG == 2", fixed = TRUE)
})

test_that("a per-CMT dispatch with a gap keeps its CMT=N keys", {
  chain <- list(
    list(kind = "assign", lhs = "CONC",  rhs = quote(CENT/VC), pos = 1L),
    list(kind = "assign", lhs = "PERIF", rhs = quote(PERI),    pos = 2L),
    list(kind = "assign", lhs = "IPRED", rhs = quote(CONC),    pos = 3L),
    cnd("CMT == 2", c("IPRED", "PERIF")),
    list(kind = "assign", lhs = "W1", rhs = 0, pos = 5L),
    list(kind = "assign", lhs = "W2", rhs = 0, pos = 6L),
    cnd("CMT == 1", c("W1", "1")),
    cnd("CMT == 2", c("W2", "1")))
  aux <- c("CENT", "PERI", "CONC", "PERIF", "IPRED", "EPS1", "EPS2")
  ep  <- .build_endpoints(quote(IPRED * (1 + W1 * EPS1) + W2 * EPS2 * THETA3),
                          chain, c("EPS1", "EPS2"), aux, c("CENT", "PERI"))
  out <- .assemble_endpoints(ep, c("CENT", "PERI"))
  expect_equal(out$readout$kind, "per_cmt")
  expect_true(any(grepl("CMT=1: DV ~ proportional(EPS1)", out$suggestion, fixed = TRUE)))
  expect_true(any(grepl("CMT=2: DV ~ ???", out$suggestion, fixed = TRUE)))
  expect_true(any(grepl("# <- Y = PERI + EPS2 * THETA3", out$suggestion, fixed = TRUE)))
})

test_that("a per-CMT gap in the LAST listed compartment is still expressible", {
  # The "never fold a failed case into the `else`" guard is a Form C concern:
  # per-CMT has no `else`, every entry carries its own key, so a gap in the last
  # listed compartment is as expressible as one anywhere else. Applying the guard
  # to both paths sent exactly this model back to the single-endpoint path, which
  # is where the two-ERROR report and the wrong-shaped `combined()` suggestion
  # came from.
  #
  # Both fixtures written for the gap put the failure FIRST, so neither could
  # show this. The position is the whole point of this one.
  chain <- list(
    list(kind = "assign", lhs = "CONC",  rhs = quote(CENT/VC), pos = 1L),
    list(kind = "assign", lhs = "RESP",  rhs = quote(PD),      pos = 2L),
    list(kind = "assign", lhs = "IPRED", rhs = quote(CONC),    pos = 3L),
    cnd("CMT == 2", c("IPRED", "RESP")),
    list(kind = "assign", lhs = "W1", rhs = 0, pos = 5L),
    list(kind = "assign", lhs = "W2", rhs = 0, pos = 6L),
    cnd("CMT == 1", c("W1", "1")),
    cnd("CMT == 2", c("W2", "1")))
  aux <- c("CENT", "PD", "CONC", "RESP", "IPRED", "EPS1", "EPS2")
  # CMT == 2 is seen first, so CMT == 1 is the LAST listed value -- and it is the
  # one carrying the scaled sigma.
  ep <- .build_endpoints(quote(IPRED * (1 + W2 * EPS2) + W1 * EPS1 * THETA3),
                         chain, c("EPS1", "EPS2"), aux, c("CENT", "PD"))
  expect_equal(unlist(lapply(ep$cases[1:2], function(c) c$value)), c(2, 1))
  expect_true(is.na(ep$cases[[2]]$type))

  out <- .assemble_endpoints(ep, c("CENT", "PD"))
  expect_false(is.null(out$suggestion))
  expect_equal(out$readout$kind, "per_cmt")
  expect_true(any(grepl("CMT=1: DV ~ ???", out$suggestion, fixed = TRUE)))
  expect_true(any(grepl("CMT=2: DV ~ proportional(EPS2)", out$suggestion, fixed = TRUE)))

  # Form C keeps the guard: there the last branch becomes the bare `else`, which
  # takes no condition, so a gap placed there could not name its endpoint.
  ep$col <- "FLAG"
  for (i in seq_along(ep$cases)) ep$cases[[i]]$why <-
    sub("CMT", "FLAG", ep$cases[[i]]$why %||% NA_character_)
  ep$cases[[1]]$why <- NULL
  ep$cases[[3]]$why <- NULL
  out_c <- .assemble_endpoints(ep, c("CENT", "PD"))
  expect_null(out_c$suggestion)
  expect_match(out_c$why, "FLAG == 1", fixed = TRUE)
})

test_that("every failing endpoint is reported, not only the first", {
  # Two gaps in the file and one named in the console is a reader filling in the
  # branch they were told about and being surprised by the other.
  chain <- list(
    list(kind = "assign", lhs = "IPRED", rhs = quote(CENT/VC), pos = 1L),
    list(kind = "assign", lhs = "W1", rhs = 0, pos = 2L),
    list(kind = "assign", lhs = "W2", rhs = 0, pos = 3L),
    list(kind = "assign", lhs = "W3", rhs = 1, pos = 4L),
    cnd("FLAG == 1", c("W1", "1"), c("W3", "0")),
    cnd("FLAG == 2", c("W2", "1"), c("W3", "0")))
  aux <- c("CENT", "IPRED", "EPS1", "EPS2")
  # FLAG 1 and 2 both scale their sigma; the fall-through keeps a clean EPS1.
  ep <- .build_endpoints(
    quote(IPRED + W1 * EPS1 * THETA3 + W2 * EPS2 * THETA3 + W3 * EPS1),
    chain, c("EPS1", "EPS2"), aux, "CENT")
  out <- .assemble_endpoints(ep, "CENT")
  expect_length(out$why, 2L)
  expect_match(out$why[[1]], "FLAG == 1", fixed = TRUE)
  expect_match(out$why[[2]], "FLAG == 2", fixed = TRUE)
  expect_equal(sum(grepl("DV ~ ???", out$suggestion, fixed = TRUE)), 2L)
})

test_that("every failing endpoint is reported on the per-CMT path too", {
  # The per-CMT and Form C returns are separate call sites with the same
  # expression in them, so one test cannot guard both -- reverting only the
  # per-CMT one left the suite green.
  chain <- list(
    list(kind = "assign", lhs = "CONC",  rhs = quote(CENT/VC), pos = 1L),
    list(kind = "assign", lhs = "RESP",  rhs = quote(PD),      pos = 2L),
    list(kind = "assign", lhs = "IPRED", rhs = quote(CONC),    pos = 3L),
    cnd("CMT == 2", c("IPRED", "RESP")),
    list(kind = "assign", lhs = "W1", rhs = 0, pos = 5L),
    list(kind = "assign", lhs = "W2", rhs = 0, pos = 6L),
    cnd("CMT == 1", c("W1", "1")),
    cnd("CMT == 2", c("W2", "1")))
  aux <- c("CENT", "PD", "CONC", "RESP", "IPRED", "EPS1", "EPS2")
  # BOTH endpoints scale their sigma.
  ep <- .build_endpoints(
    quote(IPRED + W1 * EPS1 * THETA3 + W2 * EPS2 * THETA3),
    chain, c("EPS1", "EPS2"), aux, c("CENT", "PD"))
  out <- .assemble_endpoints(ep, c("CENT", "PD"))
  expect_length(out$why, 2L)
  expect_true(any(grepl("CMT == 2", out$why, fixed = TRUE)))
  expect_true(any(grepl("CMT == 1", out$why, fixed = TRUE)))
  expect_equal(sum(grepl("DV ~ ???", out$suggestion, fixed = TRUE)), 2L)
})

# -- issue #25: NONMEM compartment numbering --------------------------------

test_that(".nm_cmt_order permutes the state list into $MODEL COMP order", {
  # The whole point: ferx numbers compartments by position in `states=[...]`,
  # NONMEM by $MODEL COMP order, and nonmem2rx hands back $DES order.
  expect_equal(
    .nm_cmt_order(state_raw   = c("CENTRAL", "DEPOT"),
                  state_names = c("CENTRAL", "DEPOT"),
                  comps       = c("DEPOT", "CENTRAL")),
    c("DEPOT", "CENTRAL"))
})

test_that(".nm_cmt_order leaves an order that already agrees alone", {
  # Must be identical(), not merely equal: the caller decides whether to announce
  # a renumbering by comparing against the input, so a spurious attribute or a
  # name carried through here would produce an INFO about a change that did not
  # happen. unname() in the helper is what makes this hold.
  st <- c("DEPOT", "CENTRAL")
  expect_identical(.nm_cmt_order(st, st, c("DEPOT", "CENTRAL")), st)
})

test_that(".nm_cmt_order matches raw names but returns sanitised ones", {
  # The trap this guards: matching is done on RAW d/dt names because that is what
  # .same_cmt_name() compares against a $MODEL COMP name -- it strips a leading
  # `c.` from one side only, so `c.RTOT` vs its own sanitised form `c_RTOT`
  # compares FALSE. Match on raw, apply by index. Mixing the two lists up
  # reorders garbage and still returns a character vector of the right length.
  expect_equal(
    .nm_cmt_order(state_raw   = c("CENT", "c.RTOT"),
                  state_names = c("CENT", "c_RTOT"),
                  comps       = c("RTOT", "CENT")),
    c("c_RTOT", "CENT"))
})

test_that(".nm_cmt_order declines rather than guessing a partial permutation", {
  st <- c("DEPOT", "CENTRAL")
  # No COMP list at all -- a non-NONMEM source. d/dt order is all there is, and
  # it is correct there, because those languages have no separate numbering.
  expect_null(.nm_cmt_order(st, st, NULL))
  expect_null(.nm_cmt_order(st, st, character()))
  # A $MODEL compartment nonmem2rx dropped: 3 declared, 2 states (issue #26).
  # Note this one is ALSO caught by the per-COMP check below (DUMMY matches no
  # state), so it does not on its own show the length guard doing anything.
  expect_null(.nm_cmt_order(st, st, c("DEPOT", "CENTRAL", "DUMMY")))
  # This one does, and nothing else catches it. A repeated COMP name gives
  # perm = c(1, 2, 1): every COMP finds exactly one state, so the per-COMP check
  # passes, and setequal() ignores duplicates so the bijection check passes too.
  # Without the length guard the helper returns a THREE-element state list for a
  # two-state model -- a longer states=[...] than the model has compartments.
  # Verified by mutation: removing the guard makes only this line fail.
  expect_null(.nm_cmt_order(st, st, c("DEPOT", "CENTRAL", "DEPOT")))
  # A COMP that matches no state.
  expect_null(.nm_cmt_order(st, st, c("DEPOT", "PERIPH")))
  # A state no COMP claims -- caught by the bijection check, not the per-COMP
  # one: every COMP below finds exactly one state, and DEPOT is still orphaned.
  expect_null(.nm_cmt_order(st, st, c("CENTRAL", "CENTRAL")))
})

test_that(".nm_cmt_order declines when two COMP names normalise alike", {
  # `.same_cmt_name()` folds case and illegal characters, so `C-ENT` and `C.ENT`
  # are the same string to it. The ordinal they would share is ambiguous, and
  # picking one is exactly the guessing this exists to remove.
  expect_null(.nm_cmt_order(state_raw   = c("C.ENT", "C-ENT"),
                            state_names = c("C_ENT", "C_ENT2"),
                            comps       = c("C.ENT", "C-ENT")))
})

test_that(".cmt_index resolves a compartment name case-insensitively", {
  # Kept as its own function so neither call site drifts into a bare match(),
  # which would be silently case-sensitive against a lowercased nonmem2rx name.
  expect_equal(.cmt_index("central", c("DEPOT", "CENTRAL")), 2L)
  expect_equal(.cmt_index("DEPOT",   c("DEPOT", "CENTRAL")), 1L)
  expect_true(is.na(.cmt_index("PERIPH", c("DEPOT", "CENTRAL"))))
  expect_true(is.na(.cmt_index("DEPOT", NULL)))
  expect_true(is.na(.cmt_index(NA_character_, c("DEPOT"))))
})


# -- issue #32: does the source's DV expression already carry S<n>? -----------

test_that(".dv_is_scaled reports NA when there is no DV expression to read", {
  # NA means "leave the existing behaviour alone", NOT "unscaled". Both guesses
  # are wrong by a factor of S<n> and neither is safe, so the caller must be able
  # to tell "the source does not scale" from "I could not tell". Returning FALSE
  # here would silently suppress scaling for any model whose DV we fail to find.
  expect_true(is.na(.dv_is_scaled(list())))
  expect_true(is.na(.dv_is_scaled(list(quote(x <- 1)))))
})

test_that(".dv_is_scaled follows the DV chain through f, unlike .explicit_obs_states", {
  # `y <- f * (1 + eps1)` with `f <- CENTRAL/scale2` is NONMEM's `Y = F`, which
  # IS scaled. .explicit_obs_states() deliberately stops at `f` -- it asks which
  # compartment is named outright, and nonmem2rx's synthetic `f` names one in
  # every model. This asks whether the value is scaled, and `f` is exactly where
  # the scaling lives. Two questions, two walks.
  lst <- list(quote(f <- CENTRAL/scale2), quote(y <- f * (1 + eps1)))
  expect_true(.dv_is_scaled(lst))
})

test_that(".dv_is_scaled distinguishes the three DV shapes", {
  # Anchored against NONMEM 7.6.0, tests/nonmem-anchor/: S<n> applies to F and
  # not to a bare A(n), a ratio of exactly V at every timepoint.
  scaled_f    <- list(quote(f <- CENTRAL/scale2), quote(y <- f * (1 + eps1)))
  scaled_expr <- list(quote(f <- CENTRAL/scale2), quote(y <- CENTRAL/scale2 * (1 + eps1)))
  bare        <- list(quote(f <- CENTRAL/scale2), quote(y <- CENTRAL * (1 + eps1)))
  expect_true(.dv_is_scaled(scaled_f))
  expect_true(.dv_is_scaled(scaled_expr))
  expect_false(.dv_is_scaled(bare))
})

test_that(".dv_is_scaled does not mistake a name containing 'scale' for scale<n>", {
  # `scalef` is an ordinary parameter. A substring match would read this model as
  # already scaled; the suffix must be digits and nothing else.
  lst <- list(quote(f <- CENTRAL/scale2),
              quote(y <- CENTRAL * scalef * (1 + eps1)))
  expect_false(.dv_is_scaled(lst))
})

# -- FIX on variance parameters (#31) -----------------------------------------
#
# Three extractors read the flag from THREE DIFFERENT channels, measured:
#   omega / kappa  -> iniDf$fix
#   sigma (nlmixr2) -> iniDf$fix on the err != NA rows
#   sigma (NONMEM)  -> attr(ui$sigma, "lotriFix"), because iniDf carries no
#                      sigma row at all for a nonmem2rx source
# A fixture covering one of them leaves the other two unguarded, so each gets
# its own test against its own channel.

# Minimal iniDf. Only the columns the extractors read; building the real 12-column
# frame would assert nothing extra and would drift with upstream.
.mk_ini <- function(...) {
  rows <- list(...)
  do.call(rbind, lapply(rows, function(r) data.frame(
    ntheta = r$ntheta %||% NA_integer_, neta1 = r$neta1 %||% NA_integer_,
    neta2  = r$neta2  %||% NA_integer_, name = r$name, est = r$est,
    fix = r$fix, err = r$err %||% NA_character_,
    condition = r$condition %||% NA_character_,
    stringsAsFactors = FALSE)))
}

test_that(".extract_omegas carries iniDf$fix onto diagonal entries", {
  ini <- .mk_ini(
    list(neta1 = 1, neta2 = 1, name = "eta.cl", est = 0.09, fix = TRUE,  condition = "id"),
    list(neta1 = 2, neta2 = 2, name = "eta.v",  est = 0.04, fix = FALSE, condition = "id")
  )
  om <- .extract_omegas(ini)$omegas
  expect_true(om[[1]]$fixed)
  expect_false(om[[2]]$fixed)
})

test_that(".extract_omegas fixes a block only when every element is fixed", {
  blk <- function(f1, f12, f2) .extract_omegas(.mk_ini(
    list(neta1 = 1, neta2 = 1, name = "eta.cl",         est = 0.09, fix = f1,  condition = "id"),
    list(neta1 = 2, neta2 = 1, name = "(eta.cl,eta.v)", est = 0.01, fix = f12, condition = "id"),
    list(neta1 = 2, neta2 = 2, name = "eta.v",          est = 0.04, fix = f2,  condition = "id")
  ))
  all_fixed <- blk(TRUE, TRUE, TRUE)
  expect_true(all_fixed$omegas[[1]]$fixed)
  expect_length(all_fixed$warnings, 0L)

  none_fixed <- blk(FALSE, FALSE, FALSE)
  expect_false(none_fixed$omegas[[1]]$fixed)
  expect_length(none_fixed$warnings, 0L)
})

test_that("a partially fixed block is emitted free and reported, not silently coerced", {
  # NOT reachable from a control stream or a model function: NONMEM's FIX is a
  # $OMEGA RECORD attribute and nlmixr2's fix() wraps the whole `a + b ~ ...`
  # line, so both sources set the flag uniformly (measured). The iniDf is built
  # by hand here for that reason -- inventing a source spelling for it would be
  # a fixture that cannot fail.
  mixed <- .extract_omegas(.mk_ini(
    list(neta1 = 1, neta2 = 1, name = "eta.cl",         est = 0.09, fix = TRUE,  condition = "id"),
    list(neta1 = 2, neta2 = 1, name = "(eta.cl,eta.v)", est = 0.01, fix = FALSE, condition = "id"),
    list(neta1 = 2, neta2 = 2, name = "eta.v",          est = 0.04, fix = FALSE, condition = "id")
  ))
  # Emitted FREE: ferx fixes a block_omega all or nothing, so the flag cannot
  # describe this block, and free is the direction the user can see in the fit.
  expect_false(mixed$omegas[[1]]$fixed)
  expect_match(mixed$warnings, "^ERROR \\| ")
  expect_match(mixed$warnings, "mixes fixed", fixed = TRUE)
  expect_match(mixed$warnings, "ETA_CL", fixed = TRUE)
  expect_length(mixed$unsupported, 1L)
  expect_match(mixed$unsupported, "all or nothing", fixed = TRUE)
})

test_that(".extract_kappas carries fix from iniDf, not from the omega attribute", {
  # Measured on rxode2 5.1.2: an IOV model has the flag on its iniDf row while
  # attr(ui$omega, "lotriFix") is NULL for that same model. Reading the
  # attribute here would answer FALSE for every fixed kappa.
  ka <- .extract_kappas(.mk_ini(
    list(neta1 = 1, neta2 = 1, name = "iov.cl", est = 0.05, fix = TRUE,  condition = "occ"),
    list(neta1 = 2, neta2 = 2, name = "iov.v",  est = 0.03, fix = FALSE, condition = "occ")
  ))$kappas
  expect_true(ka[[1]]$fixed)
  expect_false(ka[[2]]$fixed)
})

test_that(".extract_sigmas reads fix from iniDf on the nlmixr2 branch", {
  sg <- .extract_sigmas(.mk_ini(
    list(ntheta = 1, name = "prop.sd", est = 0.2, fix = TRUE,  err = "prop", condition = "cp"),
    list(ntheta = 2, name = "add.sd",  est = 0.5, fix = FALSE, err = "add",  condition = "cp")
  ))$sigmas
  expect_true(sg[[1]]$fixed)
  expect_false(sg[[2]]$fixed)
})

test_that(".extract_sigmas reads lotriFix on the nonmem2rx branch", {
  # iniDf is byte-identical for `$SIGMA 0.04` and `$SIGMA 0.04 FIX` -- it has no
  # sigma row at all for a NONMEM source -- so the matrix attribute is the only
  # channel. Empty ini, exactly as the nonmem2rx path presents it.
  empty <- .mk_ini(list(neta1 = 1, neta2 = 1, name = "eta.cl", est = 0.09,
                        fix = FALSE, condition = "id"))
  m <- matrix(c(0.04, 0, 0, 0.01), 2, 2,
              dimnames = list(c("eps1", "eps2"), c("eps1", "eps2")))
  attr(m, "lotriFix") <- matrix(c(TRUE, FALSE, FALSE, FALSE), 2, 2,
                                dimnames = dimnames(m))
  sg <- .extract_sigmas(empty, m)$sigmas
  expect_true(sg[[1]]$fixed)
  expect_false(sg[[2]]$fixed)

  # A missing or mis-shaped attribute must answer FALSE, not error. It is an
  # upstream attribute this package does not control.
  for (bad in list(NULL, "nonsense", matrix(TRUE, 1, 1))) {
    m2 <- m
    attr(m2, "lotriFix") <- bad
    sg2 <- expect_no_error(.extract_sigmas(empty, m2)$sigmas)
    expect_false(sg2[[1]]$fixed)
  }
})

test_that(".ini_fixed answers FALSE for an iniDf with no fix column", {
  # NULL[i] is NULL and isTRUE(NULL) is FALSE, so an absent column already
  # yields the right answer -- but by accident, at whichever call site touches
  # it first. This pins the decision.
  no_fix <- data.frame(name = "eta.cl", est = 0.09, stringsAsFactors = FALSE)
  expect_false(.ini_fixed(no_fix, 1L))
})

test_that("$SIGMA BLOCK residual covariances are reported rather than dropped in silence", {
  # Only the diagonal of ui$sigma is emitted, so the off-diagonals are lost.
  # That drop predates #31 and is not fixed here; the silence is. Sigma was the
  # only one of the three channels that said nothing -- omega emits its block in
  # full and .extract_kappas() already warns when it drops an off-diagonal.
  empty <- .mk_ini(list(neta1 = 1, neta2 = 1, name = "eta.cl", est = 0.09,
                        fix = FALSE, condition = "id"))
  m <- matrix(c(0.04, 0.001, 0.001, 0.01), 2, 2,
              dimnames = list(c("eps1", "eps2"), c("eps1", "eps2")))
  out <- .extract_sigmas(empty, m)
  expect_match(out$warnings, "^ERROR \\| ")
  expect_match(out$warnings, "0.001", fixed = TRUE)
  expect_match(out$warnings, "block_sigma", fixed = TRUE)
  expect_length(out$unsupported, 1L)

  # A diagonal $SIGMA must stay quiet, or every ordinary model gains a warning.
  diag_only <- .extract_sigmas(empty, matrix(c(0.04, 0, 0, 0.01), 2, 2,
    dimnames = list(c("eps1", "eps2"), c("eps1", "eps2"))))
  expect_length(diag_only$warnings, 0L)
  expect_length(diag_only$unsupported, 0L)
})
