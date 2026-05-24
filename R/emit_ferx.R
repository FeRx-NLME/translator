#' Emit a .ferx file from a ferx intermediate representation
#'
#' Takes a validated `ferx_ir` and returns a single character string ready to
#' write as a `.ferx` file. Does not perform file I/O; use [write_ferx()] for that.
#'
#' @param ir A `ferx_ir` object.
#'
#' @return A single character string (the complete `.ferx` file content).
#'
#' @seealso [new_ferx_ir()], [validate_ferx_ir()], [write_ferx()]
#'
#' @examples
#' ir <- new_ferx_ir(
#'   source_format = "nonmem",
#'   thetas      = list(list(name = "TVCL", init = 0.134, lower = 0.001, upper = 10)),
#'   omegas      = list(list(type = "diagonal", names = "ETA_CL", values = 0.07)),
#'   sigmas      = list(list(name = "PROP_ERR", value = 0.01, scale = "sd")),
#'   indiv_params = list(list(lhs = "CL", rhs = "TVCL * exp(ETA_CL)")),
#'   structural  = list(type = "pk_macro", pk_call = "one_cpt_oral",
#'                      pk_args = list(cl = "CL", v = "V", ka = "KA")),
#'   error_model = list(list(dv = "DV", type = "proportional", params = "PROP_ERR")),
#'   fit_options = list(method = "foce", maxiter = 300L, covariance = TRUE)
#' )
#' cat(emit_ferx(ir))
#' @export
emit_ferx <- function(ir) {
  validate_ferx_ir(ir)

  parts <- list(
    .emit_header(ir),
    .emit_parameters_section(ir),
    if (length(ir$indiv_params)  > 0) .emit_indiv_params_section(ir),
    if (length(ir$structural)    > 0) .emit_structural_section(ir),
    if (length(ir$odes)          > 0) .emit_odes_section(ir),
    if (length(ir$diffusion)     > 0) .emit_diffusion_section(ir),
    if (length(ir$error_model)   > 0) .emit_error_model_section(ir),
    if (!is.null(ir$scaling$obs_scale)) .emit_scaling_section(ir),
    if (length(ir$fit_options)   > 0) .emit_fit_options_section(ir)
  )

  paste(Filter(Negate(is.null), parts), collapse = "\n\n")
}

# -- helpers ------------------------------------------------------------------

.fmt_num <- function(x) {
  if (!is.finite(x)) return(if (x > 0) "1e15" else "-1e15")
  s <- format(x, scientific = FALSE, trim = TRUE, digits = 15)
  if (!grepl("\\.", s)) paste0(s, ".0") else s
}

.fmt_opt <- function(v) {
  if (is.logical(v)) tolower(as.character(v)) else as.character(v)
}

.emit_header <- function(ir) {
  src  <- if (!is.na(ir$source_format)) ir$source_format else "unknown"
  file <- if (!is.na(ir$source_file))   ir$source_file   else "unknown"
  out  <- paste0("# Translated from ", src, ": ", file)
  nw   <- length(ir$warnings)
  if (nw > 0)
    out <- paste0(out, "\n# Warnings: ", nw,
                  " -- run result$warnings for details")
  nu <- length(ir$unsupported)
  if (nu > 0) {
    warn_lines <- paste0("# WARNING: ", ir$unsupported, collapse = "\n")
    out <- paste0(out, "\n", warn_lines)
  }
  out
}

.emit_parameters_section <- function(ir) {
  groups <- list(
    vapply(ir$thetas,  .emit_theta,  ""),
    unlist(lapply(ir$omegas, .emit_omega)),
    vapply(ir$kappas,  .emit_kappa,  ""),
    vapply(ir$sigmas,  .emit_sigma,  "")
  )
  groups <- Filter(function(g) length(g) > 0, groups)
  body   <- paste(unlist(lapply(groups, function(g) c(g, ""))), collapse = "\n")
  body   <- sub("\n$", "", body)
  paste0("[parameters]\n", body)
}

.emit_theta <- function(t) {
  sprintf("  theta %s(%s, %s, %s)",
          t$name, .fmt_num(t$init), .fmt_num(t$lower), .fmt_num(t$upper))
}

.emit_omega <- function(o) {
  if (identical(o$type, "block")) {
    names_str  <- paste(o$names, collapse = ", ")
    values_str <- paste(vapply(o$values, .fmt_num, ""), collapse = ", ")
    sprintf("  block_omega (%s) = [%s]", names_str, values_str)
  } else {
    sprintf("  omega %s ~ %s", o$names, .fmt_num(o$values))
  }
}

.emit_kappa <- function(k) {
  sprintf("  kappa %s ~ %s", k$name, .fmt_num(k$value))
}

.emit_sigma <- function(s) {
  suffix <- if (identical(s$scale, "sd")) " (sd)" else ""
  sprintf("  sigma %s ~ %s%s", s$name, .fmt_num(s$value), suffix)
}

.emit_indiv_params_section <- function(ir) {
  lines <- vapply(ir$indiv_params,
                  function(p) paste0("  ", p$lhs, " = ", p$rhs), "")
  paste0("[individual_parameters]\n", paste(lines, collapse = "\n"))
}

.emit_structural_section <- function(ir) {
  s    <- ir$structural
  body <- if (identical(s$type, "pk_macro")) {
    args_str <- paste(names(s$pk_args), unlist(s$pk_args), sep = "=",
                      collapse = ", ")
    paste0("  pk ", s$pk_call, "(", args_str, ")")
  } else {
    states_str <- paste(s$states, collapse = ", ")
    paste0("  ode(obs_cmt=", s$obs_cmt, ", states=[", states_str, "])")
  }
  paste0("[structural_model]\n", body)
}

.emit_odes_section <- function(ir) {
  lines <- vapply(ir$odes,
                  function(o) paste0("  d/dt(", o$state, ") = ", o$rhs), "")
  paste0("[odes]\n", paste(lines, collapse = "\n"))
}

.emit_diffusion_section <- function(ir) {
  lines <- vapply(ir$diffusion,
                  function(d) paste0("  ", d$state, " ~ ", .fmt_num(d$value)), "")
  paste0("[diffusion]\n", paste(lines, collapse = "\n"))
}

.emit_error_model_section <- function(ir) {
  lines <- vapply(ir$error_model, function(e) {
    paste0("  ", e$dv, " ~ ", e$type, "(", paste(e$params, collapse = ", "), ")")
  }, "")
  paste0("[error_model]\n", paste(lines, collapse = "\n"))
}

.emit_scaling_section <- function(ir) {
  paste0("[scaling]\n  obs_scale = ", .fmt_num(ir$scaling$obs_scale))
}

.emit_fit_options_section <- function(ir) {
  opts         <- ir$fit_options
  known_order  <- c("method", "maxiter", "covariance", "gradient_method",
                    "iov_column")
  keys <- c(intersect(known_order, names(opts)), setdiff(names(opts), known_order))
  lines <- vapply(keys,
                  function(k) paste0("  ", k, " = ", .fmt_opt(opts[[k]])), "")
  paste0("[fit_options]\n", paste(lines, collapse = "\n"))
}
