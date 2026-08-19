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

test_that("an ODE intermediate is inlined with the binding it had when written", {
  # Kills the mutation that reverts pass 2b to re-normalising the raw expression
  # with the FINAL name_map. De-shadowing makes that map time-varying, so `frac`
  # -- written when `cl` still meant the theta -- was inlined as the individual
  # parameter. Both forms parse; ferx cannot tell them apart.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ cl <- 1.0; v <- 10.0; ka <- 1; prop.err <- 0.1; eta.cl ~ 0.09 })
    model({ frac <- central/cl; cl <- cl*exp(eta.cl); v <- v; ka <- ka
            d/dt(depot)   = -ka*depot
            d/dt(central) =  ka*depot - cl/v*central - frac
            central ~ prop(prop.err) })
  }
  ir  <- suppressWarnings(rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2"))
  rhs <- vapply(ir$odes, function(o) o$rhs, "")[2]
  expect_match(rhs, "central/TVCL", fixed = TRUE)
  expect_false(grepl("central/CL", rhs, fixed = TRUE))
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
  rhs <- paste(vapply(ir$odes, function(o) o$rhs, ""), collapse = " ")
  expect_false(grepl("C_2", rhs, fixed = TRUE))   # never referenced undeclared
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

  expect_equal(vapply(ir$odes, function(o) o$state, ""), "central")
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
  # The RHS is the half that was wrong, and asserting only the two names above
  # let it stay wrong: the declaration kept CENT while every reference to it
  # resolved through name_map to the theta, emitting `d/dt(CENT) = -K * VC` --
  # the compartment amount replaced by a fixed theta, no warning, engine clean.
  # Inside [odes] the state wins; in $PK the same symbol still means the theta.
  expect_equal(ir$odes[[1]]$rhs, "-K * CENT")
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
  f <- withr::local_tempfile(fileext = ".ctl")
  writeLines(c("$PROBLEM x", "$MODEL",
               "  COMP=(CENT, DEFDOSE, DEFOBS)", "  COMP=(PERIPH)", "$PK"), f)
  out <- .extract_nm_defobs(f)
  expect_equal(out$index, 1L)
  expect_equal(out$name, "CENT")
  expect_equal(out$n_comp, 2L)
})

test_that(".extract_nm_defobs tolerates the legal $MODEL spellings", {
  mk <- function(...) {
    f <- withr::local_tempfile(fileext = ".ctl", .local_envir = parent.frame())
    writeLines(c("$PROBLEM x", "$MODEL", ..., "$PK"), f); f
  }
  # `COMP (...)` without `=`, space-separated attributes.
  expect_equal(.extract_nm_defobs(mk("  COMP (DEPOT DEFDOSE)",
                                     "  COMP (CENTRAL DEFOBS)"))$index, 2L)
  # DEFOBSERVATION, the unabbreviated spelling.
  expect_equal(.extract_nm_defobs(mk("  COMP=(A)",
                                     "  COMP=(B, DEFOBSERVATION)"))$index, 2L)
  # DEFDOSE must not be mistaken for DEFOBS -- they share the DEF prefix.
  expect_null(.extract_nm_defobs(mk("  COMP=(A, DEFDOSE)", "  COMP=(B)")))
  # A DEFOBS that only appears in a comment is not a declaration.
  expect_null(.extract_nm_defobs(mk("  COMP=(A) ; DEFOBS goes here one day",
                                    "  COMP=(B)")))
  # No $MODEL at all.
  f <- withr::local_tempfile(fileext = ".ctl")
  writeLines(c("$PROBLEM x", "$PK"), f)
  expect_null(.extract_nm_defobs(f))
})

test_that("DEFOBS decides obs_cmt AND the scaling compartment", {
  skip_if_not_installed("nonmem2rx")
  # The regression this exists for: with DEFOBS declared first, the old
  # tail(state_names, 1) guess picked PERIPH, and because the scaling lookup
  # derives its compartment number from the same guess, the correctly parsed
  # `S1 = V` was discarded with no diagnostic.
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

  # Without the hint -- the pre-fix behaviour -- both go wrong together.
  guessed <- suppressWarnings(rxui_to_ir(ui, source_format = "nonmem",
                                         scaling_hint = sc, obs_hint = NULL))
  expect_equal(guessed$structural$obs_cmt, "PERIPH")
  expect_null(guessed$scaling$obs_scale)
})

test_that("a DEFOBS that disagrees with d/dt order is refused, not trusted", {
  skip_if_not_installed("nonmem2rx")
  p  <- nm_path("defobs_not_last.ctl")
  ui <- nonmem2rx::nonmem2rx(p)
  # $MODEL says compartment 1 is NOPE; the first d/dt is for CENT. Believing the
  # index anyway would silently observe the wrong compartment, which is the
  # failure the whole fix exists to remove.
  ir <- suppressWarnings(rxui_to_ir(ui, source_format = "nonmem",
                                    obs_hint = list(index = 1L, name = "NOPE",
                                                    n_comp = 2L)))
  expect_match(ir$warnings, "orderings disagree", all = FALSE)
  expect_equal(ir$structural$obs_cmt, "PERIPH")   # fell back to the guess
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
  expect_true("bioavailability (f)" %in% ir$unsupported)
  expect_true("dose lag time (alag)" %in% ir$unsupported)
  expect_match(ir$warnings, "bioavailability for compartment 'depot'", all = FALSE)
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
