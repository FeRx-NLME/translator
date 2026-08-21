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
    if (length(ir$error_model)   > 0) .emit_error_model_section(ir)
    else if (length(ir$error_suggestion) > 0) paste(ir$error_suggestion, collapse = "\n"),
    if (length(ir$scaling) > 0) .emit_scaling_section(ir),
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

# -- statement rendering ------------------------------------------------------

# [individual_parameters] and [odes] each hold an ORDERED statement list, and
# both may contain an `if`. One renderer serves both, because the only thing
# that differs between them is which statement kind a bare entry defaults to.
#
# Entries carry an explicit `kind`; an entry without one is read as the kind
# that block used before statement lists existed -- `assign` for
# [individual_parameters], `ddt` for [odes], `init` for initial_conditions. That
# is what keeps the hand-built IRs in test-emit.R and test-ir.R working, and it
# is a compatibility rule rather than a guess: those shapes have no `kind` field
# at all, so there is nothing to misread.
.stmt_kind <- function(s, default) {
  if (!is.null(s$kind) && length(s$kind) == 1L && nzchar(s$kind)) s$kind else default
}

# The required fields per kind. Checked rather than assumed for the same reason
# an unknown `kind` aborts: without it, a statement of a KNOWN kind that is
# missing a field emits a broken line instead -- `list(kind = "assign", rhs =
# "1")` produced the line `   = 1`, which the identifier census cannot see
# (`unlist()` drops the NULL) and only the engine objects to, when it is run.
# Refusing one and mangling the other is the inconsistency; 5b builds these
# statements programmatically, which is when a missing field becomes reachable.
.STMT_REQUIRED <- list(
  assign = c("lhs", "rhs"),
  ddt    = c("state", "rhs"),
  init   = c("state", "rhs"),
  `if`   = c("cond", "then"))

.emit_stmt <- function(s, default_kind, indent = "  ") {
  kind <- .stmt_kind(s, default_kind)
  need <- .STMT_REQUIRED[[kind]]
  if (!is.null(need)) {
    absent <- need[vapply(need, function(f) {
      v <- s[[f]]
      length(v) == 0L || (is.character(v) && !nzchar(v[1L]))
    }, logical(1))]
    if (length(absent) > 0)
      cli::cli_abort(c(
        "{.val {kind}} statement is missing {.field {absent}}.",
        i = "A statement of a known kind must carry every field that kind emits."))
  }
  switch(kind,
    assign = paste0(indent, s$lhs, " = ", s$rhs),
    ddt    = paste0(indent, "d/dt(", s$state, ") = ", s$rhs),
    init   = {
      line <- paste0(indent, "init(", s$state, ") = ", s$rhs)
      # A per-line note goes ABOVE the line it explains. The header warning block
      # is the index; this is the annotation a reader needs where they are.
      if (!is.null(s$note) && nzchar(s$note)) c(paste0(indent, "# ", s$note), line)
      else                                    line
    },
    `if`   = .emit_if_stmt(s, default_kind, indent),
    cli::cli_abort(c(
      "Unknown statement kind {.val {kind}}.",
      i = "Expected one of {.val assign}, {.val ddt}, {.val init}, {.val if}."))
  )
}

# Braces are MANDATORY, on both arms, and this is measured rather than stylistic:
# ferx answers `E_PARSE: Expected `{` after if-condition` to the unbraced form.
# NONMEM writes its commonest conditional unbraced (`IF (DSC.LT.0.0) DSC = 0.0`,
# no THEN/ENDIF), so the one-statement case is the one that matters and it still
# has to be wrapped.
#
# A single simple statement per arm renders inline, which is what makes an
# emitted [odes] block diffable against the source $DES it came from. Anything
# longer, or a nested `if`, goes multi-line -- ferx accepts both layouts, as it
# does `else if` and `else` on its own line.
# Each arm is rendered EXACTLY ONCE, at the nested indent, and the inline form is
# derived from that text by trimming. Rendering once to choose the layout and
# again to emit it doubles the work at every nesting level -- 2^depth for the
# innermost statement -- which is avoidable rather than merely cheap at the
# depths real models use.
.emit_if_stmt <- function(s, default_kind, indent) {
  inner  <- paste0(indent, "  ")
  render <- function(b)
    unlist(lapply(b, .emit_stmt, default_kind = default_kind, indent = inner))

  thn <- render(s$then)
  els <- if (length(s$else_) > 0) render(s$else_) else NULL

  # One rendered line per arm, and neither arm is itself an `if`. Testing the
  # rendered text rather than the statement kind also catches the multi-line
  # shapes a single statement can take, such as an init() carrying a note.
  inline <- length(thn) == 1L && (is.null(els) || length(els) == 1L) &&
            !any(grepl("^\\s*if \\(", c(thn, els)))
  if (inline) {
    out <- paste0(indent, "if (", s$cond, ") { ", trimws(thn), " }")
    if (!is.null(els)) out <- paste0(out, " else { ", trimws(els), " }")
    return(out)
  }

  out <- c(paste0(indent, "if (", s$cond, ") {"), thn)
  if (!is.null(els)) out <- c(out, paste0(indent, "} else {"), els)
  c(out, paste0(indent, "}"))
}

.emit_indiv_params_section <- function(ir) {
  lines <- unlist(lapply(ir$indiv_params, .emit_stmt, default_kind = "assign"))
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
    # `obs_cmt` is omitted when a `y` readout carries the observation instead.
    # ferx ignores it silently in that case, so leaving it in would state a
    # compartment choice the file no longer makes.
    obs <- if (is.null(s$obs_cmt)) "" else paste0("obs_cmt=", s$obs_cmt, ", ")
    paste0("  ode(", obs, "states=[", states_str, "])")
  }
  # A note goes ABOVE the line it explains, as an init()'s does. `states=[...]`
  # is what fixes ferx's compartment numbering, so it is where a reader asks
  # which compartment a dose reaches -- and the header warning block is thirty
  # lines away and does not travel with the eye.
  if (!is.null(s$note) && nzchar(s$note))
    body <- paste0("  # ", s$note, "\n", body)
  paste0("[structural_model]\n", body)
}

.emit_odes_section <- function(ir) {
  # init() directives lead, so a reader sees each compartment's starting value
  # before its rate of change -- and because ferx parses them as [odes] lines,
  # they belong in this block rather than a section of their own.
  #
  # The rest of the block is emitted in list order and NOTHING here reorders it.
  # [odes] has no use-before-def check: an intermediate placed below the d/dt
  # line that reads it stays valid and reads a stale slot, and PRED collapses to
  # a constant with no diagnostic. Source order is the correctness property.
  init <- unlist(lapply(ir$initial_conditions, .emit_stmt, default_kind = "init"))
  if (is.null(init)) init <- character()
  lines <- unlist(lapply(ir$odes, .emit_stmt, default_kind = "ddt"))
  if (is.null(lines)) lines <- character()
  paste0("[odes]\n", paste(c(init, lines), collapse = "\n"))
}

.emit_diffusion_section <- function(ir) {
  lines <- vapply(ir$diffusion,
                  function(d) paste0("  ", d$state, " ~ ", .fmt_num(d$value)), "")
  paste0("[diffusion]\n", paste(lines, collapse = "\n"))
}

# Three shapes, chosen by what the entries carry, not by a mode flag:
# a plain single-endpoint list, a per-compartment `CMT=N:` list, and a
# covariate-selected if/else chain. The last entry of a selected chain is the
# bare `else` ferx requires "so every observation maps to an error model", and
# it is marked by a NULL `cond`.
.emit_error_model_section <- function(ir) {
  dv <- function(e) paste0(e$dv, " ~ ", e$type, "(",
                           paste(e$params, collapse = ", "), ")")
  has <- function(f) any(vapply(ir$error_model,
                                function(e) !is.null(e[[f]]), logical(1)))
  lines <- if (has("cmt")) {
    vapply(ir$error_model, function(e) paste0("  CMT=", e$cmt, ": ", dv(e)), "")
  } else if (has("cond")) {
    vapply(seq_along(ir$error_model), function(i) {
      e <- ir$error_model[[i]]
      if (is.null(e$cond)) paste0("  else { ", dv(e), " }")
      else paste0(if (i == 1L) "  if (" else "  else if (", e$cond, ") { ",
                  dv(e), " }")
    }, "")
  } else {
    vapply(ir$error_model, function(e) paste0("  ", dv(e)), "")
  }
  paste0("[error_model]\n", paste(lines, collapse = "\n"))
}

.emit_scaling_section <- function(ir) {
  sc <- ir$scaling
  lines <- character()
  # `[[` not `$`: on a list `$` partial-matches, so `sc$y` would return the
  # per-CMT entries for a model that has only those.
  if (!is.null(sc[["y"]]))
    lines <- c(lines, paste0("  y = ", sc[["y"]]))
  if (length(sc[["per_cmt"]]) > 0)
    lines <- c(lines, vapply(sc[["per_cmt"]], function(e)
      paste0("  y[CMT=", e$cmt, "] = ", e$expr), ""))
  if (!is.null(sc[["obs_scale"]])) {
    val <- sc[["obs_scale"]]
    lines <- c(lines, paste0("  obs_scale = ",
                             if (is.character(val)) val else .fmt_num(val)))
  }
  paste0("[scaling]\n", paste(lines, collapse = "\n"))
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
