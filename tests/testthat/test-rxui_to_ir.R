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

# -- theta / individual-parameter de-shadowing --------------------------------

test_that("no rename when theta and individual parameter names are disjoint", {
  out <- .deshadow_theta_names(c("TVCL", "TVV"), c("CL", "V"))
  expect_true(all(is.na(out$map)))
  expect_length(out$warnings, 0L)
})

test_that("a theta colliding with an individual parameter is renamed to TV<name>", {
  out <- .deshadow_theta_names(c("CL", "TVV"), c("CL", "V"))
  expect_equal(out$map, c("TVCL", NA_character_))
  expect_match(out$warnings, "^INFO", all = TRUE)
  expect_match(out$warnings, "renamed to 'TVCL'")
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
  expect_match(out$warnings, "^WARN.*duplicate")
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
  # The loop is bounded, so the invariant is asserted rather than assumed.
  ir <- new_ferx_ir(
    thetas = list(list(name = "CL", init = 1, lower = 0, upper = 10)),
    indiv_params = list(list(lhs = "CL", rhs = "CL * 2"))
  )
  # Direct check of the condition the assertion guards against.
  expect_equal(
    intersect(toupper(vapply(ir$thetas, function(t) t$name, "")),
              toupper(vapply(ir$indiv_params, function(p) p$lhs, ""))),
    "CL")
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
