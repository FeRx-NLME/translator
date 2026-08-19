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

  # State renames are provenance, not diagnostics, so they get their own line
  # and are deliberately NOT counted as warnings -- nothing is wrong with the
  # file. They belong in it all the same: the .ferx is what gets shared, and a
  # reader holding only it cannot map `c_RTOT` back to the $MODEL compartment or
  # the A(n) index it came from.
  #
  # States only. A renamed theta is self-describing (TVCL obviously derives from
  # CL), so listing those would add a line to every de-shadowed model to say
  # something a reader can already infer.
  ren <- ir$state_renames
  if (length(ren) > 0)
    out <- paste0(out, "\n",
                  paste0("# renamed: state ", names(ren), " -> ", unname(ren),
                         collapse = "\n"))

  # INFO entries are notes about the translation, not defects in the artefact.
  # Counting them makes a clean model advertise "# Warnings: 3" for renames a
  # reader cannot act on, so the header counts only WARN and ERROR.
  actionable <- ir$warnings[!startsWith(ir$warnings, "INFO")]
  if (length(actionable) > 0)
    out <- paste0(out, "\n# Warnings: ", length(actionable),
                  " -- run result$warnings for details")

  # CLAUDE.md requires every WARN/ERROR to appear as a `# WARNING:` comment in
  # the output, not only in result$warnings. Unsupported features come first
  # because they are the ferx-core feature-gap signal.
  notes <- c(ir$unsupported, sub("^(WARN|ERROR)\\s*\\|\\s*", "", actionable))
  if (length(notes) > 0)
    out <- paste0(out, "\n", paste0("# WARNING: ", notes, collapse = "\n"))
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
  fix_str <- if (isTRUE(t$fixed)) ", FIX" else ""
  sprintf("  theta %s(%s, %s, %s%s)",
          t$name, .fmt_num(t$init), .fmt_num(t$lower), .fmt_num(t$upper), fix_str)
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
  # init() directives lead, so a reader sees each compartment's starting value
  # before its rate of change -- and because ferx parses them as [odes] lines,
  # they belong in this block rather than a section of their own.
  init <- unlist(lapply(ir$initial_conditions, function(x) {
    line <- paste0("  init(", x$state, ") = ", x$rhs)
    # A per-line note goes ABOVE the line it explains. The header warning block
    # is the index; this is the annotation a reader needs where they are.
    if (!is.null(x$note) && nzchar(x$note)) c(paste0("  # ", x$note), line)
    else                                    line
  }))
  if (is.null(init)) init <- character()
  lines <- vapply(ir$odes,
                  function(o) paste0("  d/dt(", o$state, ") = ", o$rhs), "")
  paste0("[odes]\n", paste(c(init, lines), collapse = "\n"))
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
  val <- ir$scaling$obs_scale
  out <- if (is.character(val)) val else .fmt_num(val)
  paste0("[scaling]\n  obs_scale = ", out)
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
