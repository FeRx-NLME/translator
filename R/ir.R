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
#' @param state_renames Named character vector of ODE state renames, source
#'   name -> emitted name. Provenance, not a diagnostic: the emitted `.ferx` is
#'   the artefact that gets shared, and a reader holding only that file cannot
#'   otherwise map a sanitised state back to the source compartment it came
#'   from. Rendered as `# renamed:` comments by [emit_ferx()].
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
  unsupported   = character(),
  state_renames = character()
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
      unsupported   = unsupported,
      state_renames = state_renames
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
    # ferx rejects an observable compartment that is not a declared state
    # (E_PARSE: Observable compartment 'CENT' not in states). Emitted, it
    # produced only INFO-level warnings and an empty $unsupported -- the
    # translation looked clean and the engine refused the file.
    if (!ir$structural$obs_cmt %in% ir$structural$states)
      cli::cli_abort(c(
        "structural$obs_cmt {.val {ir$structural$obs_cmt}} is not one of structural$states.",
        i = "States are {.val {ir$structural$states}}.",
        x = "ferx rejects this with {.code E_PARSE: Observable compartment not in states}."
      ))
  }

  .validate_ir_names(ir)

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

# Every name the IR DECLARES must be a legal ferx identifier.
#
# This is the structural half of the guarantee. Legality was previously enforced
# only where each name was minted -- in .norm(), in .sanitise_state_names(), in
# .deshadow_theta_names() -- plus one corpus test that tokenised the emitted text
# of whichever models happen to be bundled. Both miss the same thing: a name that
# reaches the file through a channel nobody thought to check. An IR carrying
# `theta 1BAD`, `omega c.RTOT`, an individual parameter `A B` and a state `9CENT`
# passed validation untouched and was emitted verbatim, and the engine answered
# `E_PARSE: Expected an assignment, an 'if' block, or 'd/dt(...)', got Ident("A")`.
# Asserting it HERE means the next channel added fails the first time it emits an
# illegal name, not the first time someone bundles a model that uses one.
#
# DECLARATIONS ONLY -- deliberately not expression text. `indiv_params$rhs`,
# `odes$rhs` and `scaling$obs_scale` legitimately carry covariate references, and
# a covariate is the one name that must NOT be sanitised: ferx resolves it
# against a data column case-sensitively, so rewriting an illegal one turns a
# working reference into E_MISSING_COVARIATE at fit time. Those are reported as
# untranslatable by the translator instead. Nothing in the declared set has that
# exemption, so the check needs no allowlist -- which is why it can be an abort.
.ir_declared_names <- function(ir) {
  nm <- function(xs, f) if (length(xs) == 0) character() else unlist(lapply(xs, f))
  c(theta        = nm(ir$thetas,       function(t) t$name),
    omega        = nm(ir$omegas,       function(o) o$names),
    kappa        = nm(ir$kappas,       function(k) k$name),
    sigma        = nm(ir$sigmas,       function(s) s$name),
    indiv_param  = nm(ir$indiv_params, function(p) p$lhs),
    ode_state    = nm(ir$odes,         function(o) o$state),
    diffusion    = nm(ir$diffusion,    function(d) d$state),
    error_dv     = nm(ir$error_model,  function(e) e$dv),
    error_param  = nm(ir$error_model,  function(e) e$params),
    state        = if (is.null(ir$structural$states))  character() else ir$structural$states,
    obs_cmt      = if (is.null(ir$structural$obs_cmt)) character() else ir$structural$obs_cmt,
    pk_arg       = if (is.null(ir$structural$pk_args)) character()
                   else unlist(ir$structural$pk_args),
    state_rename = unname(ir$state_renames))
}

.validate_ir_names <- function(ir) {
  # state_renames is source name -> emitted name, and BOTH halves matter: the
  # values are emitted into `ode(states=[...])` so they must be legal, and the
  # names are rendered as the source half of a `# renamed:` comment. An unnamed
  # vector -- easy to pass to the exported constructor -- emitted
  # `# renamed: state  -> c_RTOT` with a blank source name, which is provenance
  # that provenances nothing.
  ren <- ir$state_renames
  if (length(ren) > 0) {
    if (!is.character(ren))
      cli::cli_abort("state_renames must be a character vector.")
    if (is.null(names(ren)) || !all(nzchar(names(ren))))
      cli::cli_abort(c(
        "state_renames must be NAMED, source name -> emitted name.",
        i = "It is rendered as {.code # renamed: <source> -> <emitted>}; without names the source half is blank."
      ))
  }

  declared <- .ir_declared_names(ir)
  declared <- declared[nzchar(declared)]
  bad      <- declared[!.is_ferx_ident(declared)]
  if (length(bad) > 0) {
    # Report the channel, not just the name: "c.RTOT is illegal" leaves the
    # reader hunting for which of eleven fields produced it.
    where <- sub("[0-9]*$", "", names(bad))
    cli::cli_abort(c(
      "The IR declares {length(bad)} name{?s} that {?is/are} not {?a/} legal ferx identifier{?s}.",
      setNames(paste0("{.field ", where, "}: {.val ", unname(bad), "}"),
               rep("x", length(bad))),
      i = "ferx names are letters, digits and underscore, and must not start with a digit.",
      i = "This is a ferxtranslate bug, not a source-model limitation -- the emitted file would not parse."
    ))
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
