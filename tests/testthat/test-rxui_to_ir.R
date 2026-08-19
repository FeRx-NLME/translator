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
  ini <- rbind(theta_row("t.K", 1), eta_row("eta1", 0.09, 1L))
  st  <- function(l) vapply(suppressWarnings(rxui_to_ir(mock_ui(ini, l)))$odes,
                            function(o) o$state, "")
  above <- st(list(quote(k <- t.K * exp(eta1)), quote(central <- 0),
                   ddt("central", quote(-k * central))))
  below <- st(list(quote(k <- t.K * exp(eta1)),
                   ddt("central", quote(-k * central)), quote(central <- 0)))
  expect_equal(above, below)
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

test_that("a state whose source name is also a parameter key is refused loudly", {
  # Writing that rename into name_map would rewrite the PARAMETER's references
  # along with the state's, because both deparse to the same symbol. Refusing
  # is the only safe answer, and it must be loud.
  ini <- rbind(theta_row("t.K", 1), eta_row("ETA_X", 0.09, 1L))
  lst <- list(quote(k <- t.K * exp(ETA_X)), ddt("ETA_X", quote(-k * ETA_X)))
  ir <- suppressWarnings(rxui_to_ir(mock_ui(ini, lst)))

  expect_match(ir$warnings, "^ERROR \\| state 'ETA_X' has the same source name",
               all = FALSE)
  expect_match(ir$unsupported, "state/parameter source-name collision",
               all = FALSE)
})
