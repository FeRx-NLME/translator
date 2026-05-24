#' Convert a rxode2 UI object to a ferx intermediate representation
#'
#' Accepts the rxUI S3 object returned by `rxode2::rxode2()`, `nonmem2rx::nonmem2rx()`,
#' or `monolix2rx::monolix2rx()` and converts it to a [new_ferx_ir()] ready for
#' [emit_ferx()].
#'
#' @param ui A rxUI S3 object (environment with `$iniDf` and `$lstExpr`).
#' @param source_format One of `"nonmem"`, `"nlmixr2"`, `"monolix"`, or `NA`.
#' @param source_file Path to the source file, or `NA`.
#'
#' @return A `ferx_ir` object.
#'
#' @seealso [new_ferx_ir()], [emit_ferx()], [to_ferx()]
#'
#' @examples
#' \dontrun{
#' ui <- rxode2::rxode2(function() {
#'   ini({ tvcl <- 0.134; eta.cl ~ 0.07; err.prop ~ 0.01 })
#'   model({ cl <- tvcl * exp(eta.cl); linCmt() ~ prop(err.prop) })
#' })
#' ir <- rxui_to_ir(ui, source_format = "nlmixr2")
#' cat(emit_ferx(ir))
#' }
#' @export
rxui_to_ir <- function(ui, source_format = NA_character_, source_file = NA_character_) {
  ini  <- ui$iniDf
  lst  <- ui$lstExpr
  warn <- character()
  unsp <- character()

  theta_out <- .extract_thetas(ini)
  warn      <- c(warn, theta_out$warnings)

  omega_out <- .extract_omegas(ini)
  kappa_out <- .extract_kappas(ini)
  warn      <- c(warn, kappa_out$warnings)

  sigma_out <- .extract_sigmas(ini)

  name_map  <- .norm_map_from_ini(ini)
  expr_out  <- .parse_model_exprs(lst, name_map)
  warn      <- c(warn, expr_out$warnings)
  unsp      <- c(unsp, expr_out$unsupported)

  structural <- expr_out$structural
  if (identical(structural$type, "ode")) {
    state_names <- vapply(expr_out$odes, function(o) o$state, "")
    obs_cmt     <- tryCatch(ui$central, error = function(e) NULL)
    if (is.null(obs_cmt) || !is.character(obs_cmt)) {
      obs_cmt <- tail(state_names, 1)
      warn <- c(warn, paste0("WARN  | obs_cmt could not be inferred -- guessed '",
                             obs_cmt, "', verify in [structural_model]"))
    }
    structural$states  <- state_names
    structural$obs_cmt <- obs_cmt
  }
  if (identical(structural$type, "lincmt")) {
    pk_out <- .infer_pk_macro(expr_out$indiv_params)
    warn   <- c(warn, pk_out$warnings)
    unsp   <- c(unsp, pk_out$unsupported)
    if (is.na(pk_out$pk_call)) {
      structural <- list()
    } else {
      structural <- list(type    = "pk_macro",
                         pk_call = pk_out$pk_call,
                         pk_args = pk_out$pk_args)
    }
  }

  fit_opts <- list(method = "focei", maxiter = 500L, covariance = TRUE)
  if (length(kappa_out$kappas) > 0)
    fit_opts$iov_column <- kappa_out$iov_column

  new_ferx_ir(
    source_format = source_format,
    source_file   = source_file,
    thetas        = theta_out$thetas,
    omegas        = omega_out$omegas,
    kappas        = kappa_out$kappas,
    sigmas        = sigma_out$sigmas,
    indiv_params  = expr_out$indiv_params,
    structural    = structural,
    odes          = expr_out$odes,
    error_model   = expr_out$error_model,
    fit_options   = fit_opts,
    warnings      = warn,
    unsupported   = unsp
  )
}

# -- name normalisation -------------------------------------------------------

.norm <- function(nm) toupper(gsub(".", "_", nm, fixed = TRUE))

.norm_map_from_ini <- function(ini) {
  nms <- unique(ini$name[!is.na(ini$name)])
  setNames(vapply(nms, .norm, ""), nms)
}

# Recursively substitute known parameter names in an expression.
# Does NOT touch the function-name position of a call (call[[1]]).
# Names absent from `map` pass through unchanged (preserves state names,
# covariates, and functions like exp/log).
.normalise_expr <- function(expr, map) {
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    return(as.name(if (nm %in% names(map)) map[[nm]] else nm))
  }
  if (!is.call(expr)) return(expr)
  as.call(c(list(expr[[1]]),
            lapply(as.list(expr[-1]), .normalise_expr, map = map)))
}

# -- iniDf extractors ---------------------------------------------------------

.extract_thetas <- function(ini) {
  rows <- ini[!is.na(ini$ntheta) & is.na(ini$err), , drop = FALSE]
  warn <- character()
  thetas <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, ]
    nm  <- .norm(row$name)
    if (isTRUE(row$fix))
      warn <<- c(warn, paste0("INFO  | THETA ", nm,
                              " was FIXED in source -- treated as free in ferx"))
    list(name = nm, init = row$est, lower = row$lower, upper = row$upper)
  })
  list(thetas = thetas, warnings = warn)
}

.extract_omegas <- function(ini) {
  iiv <- ini[!is.na(ini$neta1) & ini$condition == "id", , drop = FALSE]
  if (nrow(iiv) == 0) return(list(omegas = list()))

  off   <- iiv[iiv$neta1 != iiv$neta2, , drop = FALSE]
  diag  <- iiv[iiv$neta1 == iiv$neta2, , drop = FALSE]

  if (nrow(off) == 0) {
    omegas <- lapply(seq_len(nrow(diag)), function(i) {
      list(type = "diagonal", names = .norm(diag$name[i]), values = diag$est[i])
    })
    return(list(omegas = omegas))
  }

  blocks        <- .detect_blocks(off)
  block_eta_set <- unlist(blocks)
  omegas        <- list()

  for (bg in blocks) {
    lt    <- iiv[iiv$neta1 %in% bg & iiv$neta2 %in% bg, , drop = FALSE]
    lt    <- lt[order(lt$neta1, lt$neta2), ]
    nms   <- vapply(bg, function(e) {
      row <- lt[lt$neta1 == e & lt$neta2 == e, , drop = FALSE]
      .norm(row$name[1])
    }, "")
    omegas <- c(omegas, list(list(type = "block", names = nms, values = lt$est)))
  }

  for (i in seq_len(nrow(diag))) {
    if (!diag$neta1[i] %in% block_eta_set)
      omegas <- c(omegas, list(
        list(type = "diagonal", names = .norm(diag$name[i]), values = diag$est[i])
      ))
  }

  list(omegas = omegas)
}

# Union-find over off-diagonal eta pairs; returns list of sorted integer vectors.
.detect_blocks <- function(off) {
  eta_set <- sort(unique(c(off$neta1, off$neta2)))
  parent  <- setNames(as.list(eta_set), as.character(eta_set))

  .find <- function(x) {
    while (!identical(parent[[as.character(x)]], x))
      x <- parent[[as.character(x)]]
    x
  }
  .union <- function(a, b) {
    ra <- .find(a); rb <- .find(b)
    if (!identical(ra, rb)) parent[[as.character(ra)]] <<- rb
  }

  for (i in seq_len(nrow(off))) .union(off$neta1[i], off$neta2[i])

  roots  <- vapply(eta_set, .find, 1L)
  groups <- split(eta_set, roots)
  lapply(groups, sort)
}

.extract_kappas <- function(ini) {
  iov  <- ini[!is.na(ini$neta1) & ini$condition != "id", , drop = FALSE]
  warn <- character()
  if (nrow(iov) == 0) return(list(kappas = list(), iov_column = NULL, warnings = warn))

  diag <- iov[iov$neta1 == iov$neta2, , drop = FALSE]
  off  <- iov[iov$neta1 != iov$neta2, , drop = FALSE]

  if (nrow(off) > 0)
    warn <- c(warn, "WARN  | IOV block omega detected -- only diagonal kappas emitted")

  iov_col <- unique(iov$condition)
  if (length(iov_col) > 1)
    warn <- c(warn, paste0("WARN  | Multiple IOV condition columns: ",
                           paste(iov_col, collapse = ", "), " -- using first"))
  iov_col <- iov_col[1]

  kappas <- lapply(seq_len(nrow(diag)), function(i)
    list(name = .norm(diag$name[i]), value = diag$est[i])
  )
  list(kappas = kappas, iov_column = iov_col, warnings = warn)
}

.extract_sigmas <- function(ini) {
  rows <- ini[!is.na(ini$err), , drop = FALSE]
  # rxode2 iniDf stores sigma estimates as SD, not variance.
  sigmas <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, ]
    list(name = .norm(row$name), value = row$est, scale = "sd")
  })
  list(sigmas = sigmas)
}

# -- expression classifiers ---------------------------------------------------

.is_ddt <- function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("d/dt"))
}

.is_tilde <- function(expr) {
  is.call(expr) && identical(expr[[1]], quote(`~`))
}

.is_lincmt_tilde <- function(expr) {
  .is_tilde(expr) && is.call(expr[[2]]) &&
    identical(expr[[2]][[1]], as.name("linCmt"))
}

.is_assignment <- function(expr) {
  is.call(expr) && as.character(expr[[1]]) %in% c("<-", "=", "->", "<<-")
}

# -- expression parser --------------------------------------------------------

.parse_model_exprs <- function(lst, name_map) {
  indiv_params <- list()
  odes         <- list()
  error_model  <- list()
  structural   <- list()
  warnings     <- character()
  unsupported  <- character()
  lincmt_seen  <- FALSE

  for (expr in lst) {
    if (.is_ddt(expr)) {
      state <- as.character(expr[[2]])
      rhs   <- deparse(.normalise_expr(expr[[3]], name_map), width.cutoff = 500L)
      odes  <- c(odes, list(list(state = state, rhs = rhs)))
      if (!identical(structural$type, "ode"))
        structural <- list(type = "ode")

    } else if (.is_lincmt_tilde(expr)) {
      lincmt_seen <- TRUE
      err_out     <- .parse_error_rhs(expr[[3]], name_map)
      error_model <- c(error_model,
                       list(list(dv = "DV", type = err_out$type,
                                 params = err_out$params)))
      warnings    <- c(warnings, err_out$warnings)
      if (!identical(structural$type, "ode"))
        structural <- list(type = "lincmt")

    } else if (.is_tilde(expr)) {
      err_out     <- .parse_error_rhs(expr[[3]], name_map)
      error_model <- c(error_model,
                       list(list(dv = "DV", type = err_out$type,
                                 params = err_out$params)))
      warnings    <- c(warnings, err_out$warnings)

    } else if (.is_assignment(expr)) {
      lhs_raw  <- as.character(expr[[2]])
      lhs_norm <- .norm(lhs_raw)
      name_map[lhs_raw] <- lhs_norm
      rhs_norm <- deparse(.normalise_expr(expr[[3]], name_map), width.cutoff = 500L)
      indiv_params <- c(indiv_params,
                        list(list(lhs = lhs_norm, rhs = rhs_norm)))
    }
  }

  list(
    indiv_params = indiv_params,
    odes         = odes,
    error_model  = error_model,
    structural   = structural,
    warnings     = warnings,
    unsupported  = unsupported
  )
}

.parse_error_rhs <- function(rhs, name_map) {
  warn <- character()
  if (!is.call(rhs))
    return(list(type = "proportional", params = character(), warnings = warn))

  fn <- as.character(rhs[[1]])

  if (fn == "prop") {
    params <- .norm(as.character(rhs[[2]]))
    return(list(type = "proportional", params = params, warnings = warn))
  }
  if (fn == "add") {
    params <- .norm(as.character(rhs[[2]]))
    return(list(type = "additive", params = params, warnings = warn))
  }
  if (fn == "+") {
    lhs_fn <- tryCatch(as.character(rhs[[2]][[1]]), error = function(e) "")
    rhs_fn <- tryCatch(as.character(rhs[[3]][[1]]), error = function(e) "")
    if ((lhs_fn == "add" && rhs_fn == "prop") ||
        (lhs_fn == "prop" && rhs_fn == "add")) {
      add_node  <- if (lhs_fn == "add")  rhs[[2]] else rhs[[3]]
      prop_node <- if (lhs_fn == "prop") rhs[[2]] else rhs[[3]]
      params    <- c(.norm(as.character(add_node[[2]])),
                     .norm(as.character(prop_node[[2]])))
      return(list(type = "combined", params = params, warnings = warn))
    }
  }

  warn <- c(warn, "WARN  | complex $ERROR -- classified as proportional, verify")
  params <- tryCatch(.norm(as.character(rhs[[2]])), error = function(e) character())
  list(type = "proportional", params = params, warnings = warn)
}

# -- linCmt -> pk macro -------------------------------------------------------

.infer_pk_macro <- function(indiv_params) {
  lhs_lc   <- tolower(vapply(indiv_params, function(p) p$lhs, ""))
  warn     <- character()
  unsp     <- character()

  pk_call <- if ("ka" %in% lhs_lc && "q2" %in% lhs_lc) {
    unsp <- c(unsp, "three_cpt_oral (not supported in ferx)")
    warn <- c(warn, "ERROR | three_cpt_oral detected -- not supported in ferx, structural model omitted")
    NA_character_
  } else if ("ka" %in% lhs_lc && ("q" %in% lhs_lc || "q2" %in% lhs_lc)) {
    "two_cpt_oral"
  } else if ("ka" %in% lhs_lc) {
    "one_cpt_oral"
  } else if ("q2" %in% lhs_lc) {
    "three_cpt_iv_bolus"
  } else if ("q" %in% lhs_lc) {
    "two_cpt_iv_bolus"
  } else {
    "one_cpt_iv_bolus"
  }

  if (is.na(pk_call))
    return(list(pk_call = NA_character_, pk_args = list(),
                warnings = warn, unsupported = unsp))

  arg_keys <- switch(pk_call,
    one_cpt_oral       = c("cl", "v",  "ka"),
    one_cpt_iv_bolus   = c("cl", "v"),
    two_cpt_oral       = c("cl", "v1", "q", "v2", "ka"),
    two_cpt_iv_bolus   = c("cl", "v1", "q", "v2"),
    three_cpt_iv_bolus = c("cl", "v1", "q", "v2", "q2", "v3"),
    character()
  )

  # Collect optional pk args
  opt_keys <- c()
  if ("f" %in% lhs_lc)                             opt_keys <- c(opt_keys, "f")
  if (any(c("alag", "lagtime") %in% lhs_lc))       opt_keys <- c(opt_keys, "alag")
  all_keys <- c(arg_keys, opt_keys)

  lhs_uc <- vapply(indiv_params, function(p) p$lhs, "")
  pk_args <- list()
  for (key in all_keys) {
    # Find by exact lowercase match or common alias (v <-> v1)
    idx <- which(lhs_lc == key)
    if (length(idx) == 0 && key == "v1") idx <- which(lhs_lc == "v")
    if (length(idx) == 0 && key == "v")  idx <- which(lhs_lc == "v1")
    if (length(idx) == 0 && key == "alag") idx <- which(lhs_lc == "lagtime")
    if (length(idx) > 0)
      pk_args[[key]] <- lhs_uc[idx[1]]
  }

  list(pk_call = pk_call, pk_args = pk_args, warnings = warn, unsupported = unsp)
}
