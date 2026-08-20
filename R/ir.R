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
#' @param indiv_params Ordered list of statements. Each element is a list with
#'   a `kind`: `"assign"` (`lhs`, `rhs`, both character) or `"if"` (`cond`,
#'   character; `then`, and optionally `else_`, each a nested statement list).
#'   An element with no `kind` is read as an `"assign"`, which is the shape this
#'   field held before statement lists. Order is preserved exactly.
#' @param structural List describing the structural model. Must have `type`:
#'   `"pk_macro"` (add `pk_call` and `pk_args`) or `"ode"` (add `obs_cmt`
#'   and `states`). May be empty during incremental construction.
#' @param odes Ordered list of statements forming the `[odes]` block. Each
#'   element is a list with a `kind`: `"ddt"` (`state`, `rhs`), `"assign"`
#'   (`lhs`, `rhs`) for an ODE-block intermediate, or `"if"` (`cond`, `then`,
#'   and optionally `else_`). An element with no `kind` is read as a `"ddt"`,
#'   which is the shape this field held before statement lists. Used only when
#'   `structural$type == "ode"`.
#'
#'   Order is preserved exactly and nothing reorders it. ferx has no
#'   use-before-def check in `[odes]`: an intermediate placed below the `d/dt`
#'   line that reads it stays valid, reads a stale slot, and collapses the
#'   prediction to a constant with no diagnostic.
#' @param initial_conditions List of ODE initial-condition entries. Each element
#'   is a list with `state` (character, an emitted state name) and `rhs`
#'   (character, the expression). Rendered as `init(<state>) = <rhs>` inside the
#'   `[odes]` block. ferx resolves an init expression against individual
#'   parameters, other states and literals only -- not thetas -- so
#'   [rxui_to_ir()] drops and reports one that would reference anything else
#'   rather than emitting a file the engine rejects.
#' @param diffusion List of diffusion entries. Each element is a list with
#'   `state` (character) and `value` (numeric).
#' @param error_model List of error model entries. Each element is a list
#'   with `dv` (character), `type` (`"proportional"`, `"additive"`, or
#'   `"combined"`), and `params` (character vector of parameter names).
#' @param error_suggestion Character vector of comment lines rendered where the
#'   `[error_model]` block would have gone, when the source error expression
#'   could not be translated. The block itself is omitted so the engine rejects
#'   the file rather than accepting a guess; these lines carry the source
#'   expression and a plausible reading the user can uncomment after checking it.
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
  initial_conditions = list(),
  diffusion     = list(),
  error_model   = list(),
  error_suggestion = character(),
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
      initial_conditions = initial_conditions,
      diffusion     = diffusion,
      error_model   = error_model,
      error_suggestion = error_suggestion,
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
    # length must be checked, not just type. is.character() is TRUE for
    # character(0) and for a length-2 vector, and both then reach the `%in%`
    # below, whose `if` aborts with "argument is of length zero" / "the
    # condition has length > 1" -- a base-R error naming neither field. This is
    # the same trap rxui_to_ir() records for `ui$central`; it must not be
    # reintroduced in the exported validator.
    if (is.null(ir$structural$obs_cmt) || !is.character(ir$structural$obs_cmt) ||
        length(ir$structural$obs_cmt) != 1L)
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

  # Every error_suggestion line must be a comment. The suggestion is a plausible
  # reading of an error expression we could NOT translate, printed where the
  # block would have gone; one uncommented line makes the engine parse a guess as
  # though it were the translation, which is the exact failure this field exists
  # to avoid. Checked here rather than left to the emitter because the emitter
  # joins the lines verbatim and cannot tell the difference.
  if (length(ir$error_suggestion) > 0) {
    bad <- ir$error_suggestion[!grepl("^\\s*#", ir$error_suggestion)]
    if (length(bad) > 0)
      cli::cli_abort(c(
        "error_suggestion must contain only comment lines.",
        x = "Not a comment: {.val {bad[1]}}"
      ))
  }

  invisible(ir)
}

# Names a statement list DECLARES, walking into `if` bodies.
#
# Walking in is the whole point. A name assigned inside a branch is declared as
# much as one at the top level, and the census below is an abort rather than a
# warning -- so a channel that stops at the top level does not report a weaker
# result, it reports nothing and lets an illegal name reach the engine. That is
# the failure this census exists to prevent, reintroduced one nesting level down.
#
# `want` selects which half of the list to return, because [odes] declares two
# kinds of name that need different labels: `state` covers d/dt targets and
# init() targets, `assign` covers ODE-block intermediates. Reporting an illegal
# intermediate as an "ode state" would send the reader to $MODEL.
.stmt_declared <- function(stmts, default_kind, want) {
  out <- character()
  for (s in stmts) {
    kind <- .stmt_kind(s, default_kind)
    if (identical(kind, "if")) {
      out <- c(out, .stmt_declared(s$then,  default_kind, want),
                    .stmt_declared(s$else_, default_kind, want))
    } else if (identical(kind, "assign")) {
      if (identical(want, "assign")) out <- c(out, s$lhs)
    } else {
      if (identical(want, "state"))  out <- c(out, s$state)
    }
  }
  out
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
# Returns a data.frame of (channel, name) rather than a named vector. The names
# have to travel with a label saying which field produced them -- "c.RTOT is
# illegal" leaves the reader hunting through eleven fields -- and deriving that
# label by stripping digits off `c()`-generated suffixes does not survive real
# input: `pk_args` keys legitimately end in digits, so `pk_arg.v1` and
# `pk_arg.v2` both stripped to `pk_arg.v` and two different macro slots were
# reported under one label. A parallel vector cannot collapse that way.
.ir_declared_names <- function(ir) {
  nm  <- function(xs, f) unlist(lapply(xs, f))
  add <- function(channel, values) {
    values <- as.character(values)
    if (length(values) == 0) return(NULL)
    data.frame(channel = rep(channel, length(values)), name = values,
               stringsAsFactors = FALSE)
  }
  out <- list(
    add("theta",        nm(ir$thetas,       function(t) t$name)),
    add("omega",        nm(ir$omegas,       function(o) o$names)),
    add("kappa",        nm(ir$kappas,       function(k) k$name)),
    add("sigma",        nm(ir$sigmas,       function(s) s$name)),
    add("indiv_param",  .stmt_declared(ir$indiv_params, "assign", "assign")),
    add("ode state",    .stmt_declared(ir$odes, "ddt", "state")),
    add("ode intermediate", .stmt_declared(ir$odes, "ddt", "assign")),
    add("initial condition", nm(ir$initial_conditions, function(x) x$state)),
    add("diffusion",    nm(ir$diffusion,    function(d) d$state)),
    add("error dv",     nm(ir$error_model,  function(e) e$dv)),
    add("error param",  nm(ir$error_model,  function(e) e$params)),
    add("state",        ir$structural$states),
    add("obs_cmt",      ir$structural$obs_cmt),
    add("state rename", unname(ir$state_renames)))
  # pk_args is a NAMED list and both halves are emitted (`cl=CL`), so the slot
  # name is part of the label rather than a separate channel.
  pk <- ir$structural$pk_args
  if (length(pk) > 0) {
    keys <- if (is.null(names(pk))) rep("?", length(pk)) else names(pk)
    for (i in seq_along(pk))
      out <- c(out, list(add(paste0("pk_arg ", keys[i]), pk[[i]])))
  }
  do.call(rbind, Filter(Negate(is.null), out))
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
  if (!is.null(declared)) declared <- declared[nzchar(declared$name), , drop = FALSE]
  bad <- if (is.null(declared)) declared[0, ]
         else declared[!.is_ferx_ident(declared$name), , drop = FALSE]
  if (!is.null(bad) && nrow(bad) > 0) {
    # The offending name is passed as a cli VALUE, never pasted into the format
    # string. cli evaluates `{...}` in the format string as R code, and the
    # names reaching this line are exactly those .is_ferx_ident() rejected --
    # which rejects `{` -- so brace-carrying names are the whole population
    # here. Interpolated, a theta named `X{Sys.getenv("USER")}` had the
    # expression EVALUATED into the abort text, and `CL{foo}` replaced the
    # diagnostic with "Could not evaluate cli {} expression". The repo already
    # guards this elsewhere (test-translate.R, "engine text with braces is never
    # evaluated as cli syntax"); the same rule applies to any attacker- or
    # source-controlled string, and every name here is source-controlled.
    detail <- setNames(
      vapply(seq_len(nrow(bad)),
             function(i) paste0("{.field ", bad$channel[i], "}: {.val {bad$name[", i, "]}}"),
             character(1)),
      rep("x", nrow(bad)))
    cli::cli_abort(c(
      "The IR declares {nrow(bad)} name{?s} that {?is/are} not {?a/} legal ferx identifier{?s}.",
      detail,
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
