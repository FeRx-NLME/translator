#' Create a ferx intermediate representation
#'
#' Constructs a `ferx_ir` object that captures every concept expressible in a
#' `.ferx` file. All parsers produce this; the emitter consumes it. Names are
#' stored ferx-native (uppercase, underscored) from the moment of ingestion.
#'
#' @param source_format One of `"nonmem"`, `"nlmixr2"`, `"monolix"`, `"mrgsolve"`,
#'   or `NA`.
#' @param source_file Path to the source file, or `NA`.
#' @param thetas List of theta entries. Each element is a list with fields
#'   `name` (character), `init` (numeric), `lower` (numeric), `upper` (numeric).
#' @param omegas List of omega entries. Each element is a list with `type`
#'   (`"diagonal"` or `"block"`), `names` (character vector), and `values`
#'   (numeric). For `"diagonal"`: one name, one value. For `"block"`:
#'   multiple names, lower-triangular values.
#' @param kappas List of IOV kappa entries. Each element is a list with
#'   `name` (character) and `value` (numeric).
#' @param sigmas List of sigma entries. Each element is a list with `name`
#'   (character), `value` (numeric), and `scale` (`"sd"` or `"var"`).
#' @param indiv_params List of individual parameter assignments. Each element
#'   is a list with `lhs` (character) and `rhs` (character).
#' @param structural List describing the structural model. Must have `type`:
#'   `"pk_macro"` (add `pk_call` and `pk_args`) or `"ode"` (add `obs_cmt`
#'   and `states`). May be empty during incremental construction.
#' @param odes List of ODE entries. Each element is a list with `state`
#'   (character) and `rhs` (character). Used only when
#'   `structural$type == "ode"`.
#' @param diffusion List of diffusion entries. Each element is a list with
#'   `state` (character) and `value` (numeric).
#' @param error_model List of error model entries. Each element is a list
#'   with `dv` (character), `type` (`"proportional"`, `"additive"`, or
#'   `"combined"`), and `params` (character vector of parameter names).
#' @param scaling List with `obs_scale` (numeric or `NULL`).
#' @param fit_options List with named elements such as `method`, `maxiter`,
#'   `covariance`, and (when IOV is present) `iov_column`.
#' @param warnings Character vector of diagnostic messages, each prefixed
#'   with `INFO`, `WARN`, or `ERROR`.
#' @param unsupported Character vector of features detected in the source
#'   that could not be translated.
#'
#' @return A `ferx_ir` list.
#'
#' @seealso [validate_ferx_ir()], [emit_ferx()], [rxui_to_ir()]
#'
#' @examples
#' ir <- new_ferx_ir(
#'   source_format = "nonmem",
#'   thetas = list(list(name = "TVCL", init = 0.134, lower = 0.001, upper = 10)),
#'   omegas = list(list(type = "diagonal", names = "ETA_CL", values = 0.07)),
#'   structural = list(type = "pk_macro", pk_call = "one_cpt_oral",
#'                     pk_args = list(cl = "CL", v = "V", ka = "KA")),
#'   error_model = list(list(dv = "DV", type = "proportional", params = "PROP_ERR")),
#'   fit_options = list(method = "foce", maxiter = 300L, covariance = TRUE)
#' )
#' print(ir)
#' @export
new_ferx_ir <- function(
  source_format = NA_character_,
  source_file   = NA_character_,
  thetas        = list(),
  omegas        = list(),
  kappas        = list(),
  sigmas        = list(),
  indiv_params  = list(),
  structural    = list(),
  odes          = list(),
  diffusion     = list(),
  error_model   = list(),
  scaling       = list(),
  fit_options   = list(),
  warnings      = character(),
  unsupported   = character()
) {
  structure(
    list(
      source_format = source_format,
      source_file   = source_file,
      thetas        = thetas,
      omegas        = omegas,
      kappas        = kappas,
      sigmas        = sigmas,
      indiv_params  = indiv_params,
      structural    = structural,
      odes          = odes,
      diffusion     = diffusion,
      error_model   = error_model,
      scaling       = scaling,
      fit_options   = fit_options,
      warnings      = warnings,
      unsupported   = unsupported
    ),
    class = "ferx_ir"
  )
}

#' Validate a ferx_ir object
#'
#' Stops with an informative message if the IR is structurally invalid.
#' Call this after fully populating the IR, before passing to `emit_ferx()`.
#'
#' @param ir A `ferx_ir` object.
#'
#' @return `ir` invisibly, if valid.
#'
#' @seealso [new_ferx_ir()], [emit_ferx()]
#'
#' @examples
#' ir <- new_ferx_ir(
#'   structural = list(type = "pk_macro", pk_call = "one_cpt_oral",
#'                     pk_args = list(cl = "CL", v = "V", ka = "KA"))
#' )
#' validate_ferx_ir(ir)
#' @export
validate_ferx_ir <- function(ir) {
  if (!inherits(ir, "ferx_ir"))
    cli::cli_abort("{.arg ir} must be a {.cls ferx_ir} object.")

  if (length(ir$structural) > 0 && is.null(ir$structural$type))
    cli::cli_abort(
      "structural$type is missing.",
      i = 'Must be {.val pk_macro} or {.val ode}.'
    )

  if (!is.null(ir$structural$type)) {
    valid_types <- c("pk_macro", "ode")
    if (!ir$structural$type %in% valid_types)
      cli::cli_abort(
        "structural$type must be {.or {.val {valid_types}}}, not {.val {ir$structural$type}}."
      )
  }

  if (length(ir$odes) > 0 && !identical(ir$structural$type, "ode"))
    cli::cli_abort(
      "odes is non-empty but structural$type is not {.val ode}."
    )

  if (identical(ir$structural$type, "ode")) {
    if (is.null(ir$structural$states) || length(ir$structural$states) == 0)
      cli::cli_abort(
        "structural$states must be a non-empty character vector when structural$type is {.val ode}."
      )
    if (is.null(ir$structural$obs_cmt) || !is.character(ir$structural$obs_cmt))
      cli::cli_abort(
        "structural$obs_cmt must be a character scalar when structural$type is {.val ode}."
      )
  }

  if (identical(ir$structural$type, "pk_macro")) {
    if (is.null(ir$structural$pk_call) || !nzchar(ir$structural$pk_call))
      cli::cli_abort(
        "structural$pk_call must be a non-empty string when structural$type is {.val pk_macro}."
      )
    if (!is.list(ir$structural$pk_args))
      cli::cli_abort(
        "structural$pk_args must be a named list when structural$type is {.val pk_macro}."
      )
  }

  invisible(ir)
}

#' @export
print.ferx_ir <- function(x, ...) {
  src <- if (!is.na(x$source_format)) x$source_format else "unknown"
  file_part <- if (!is.na(x$source_file)) paste0(" (", x$source_file, ")") else ""
  cli::cli_h1("ferx_ir [{src}{file_part}]")

  counts <- c(
    thetas       = length(x$thetas),
    omegas       = length(x$omegas),
    kappas       = length(x$kappas),
    sigmas       = length(x$sigmas),
    indiv_params = length(x$indiv_params),
    odes         = length(x$odes)
  )
  for (nm in names(counts[counts > 0]))
    cli::cli_bullets(c("*" = paste0(nm, ": ", counts[[nm]])))

  if (!is.null(x$structural$type))
    cli::cli_bullets(c("*" = paste0("structural: ", x$structural$type)))

  nw <- length(x$warnings)
  nu <- length(x$unsupported)
  if (nw > 0) cli::cli_alert_warning("{nw} translation warning{?s}")
  if (nu > 0) cli::cli_alert_danger("{nu} unsupported feature{?s}")

  invisible(x)
}
