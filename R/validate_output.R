# Validate emitted .ferx text with the ferx engine before it reaches the user.
#
# Without this the translator returns successfully on output that
# ferx_model_validate() rejects immediately, which is how a dotted state name,
# an undefined $DES intermediate and an orphaned ALAG parameter all reached a
# user instead of a test (issue #6, defect 7).
#
# Two limits are deliberate and are reported rather than hidden:
#
#   * ferx is a Suggests dependency. When it is absent, validation is skipped
#     with an INFO warning -- an optional dependency must never turn a working
#     translation into a failure.
#   * Model-only validation cannot see covariate problems. Every block that
#     accepts a covariate treats an unknown identifier as one, so
#     `y = ... TYPO_NAME` validates ok = TRUE with no dataset and only becomes
#     E_MISSING_COVARIATE once data is supplied. When the dataset cannot be
#     resolved, say so instead of implying a clean bill of health.

# Run ferx_model_validate() over `text`. Returns NULL when validation did not
# run, otherwise a list with `ok`, `diagnostics`, `data_file`, `warnings` and
# `unsupported`.
# Indirection so the ferx-absent branch is testable. That branch is how the fast
# PR job runs -- only the `engine` job installs ferx -- so it must never turn a
# working translation into a failure, and that has to be asserted, not assumed.
.has_ferx <- function() requireNamespace("ferx", quietly = TRUE)

.validate_ferx_text <- function(text, data_file = NA_character_) {
  if (!.has_ferx())
    return(list(
      ok = NA, diagnostics = .empty_diagnostics(), data_file = NA_character_,
      warnings = paste0(
        "INFO  | ferx is not installed -- emitted .ferx was NOT validated. ",
        "Install it to have translation output checked by the engine."),
      unsupported = character()
    ))

  tmp <- tempfile(fileext = ".ferx")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(text, tmp)

  # ferx_model_validate() prints its own report unconditionally; capture it so
  # the translator speaks with one voice.
  res <- tryCatch(
    utils::capture.output(
      out <- if (is.na(data_file)) ferx::ferx_model_validate(tmp)
             else                  ferx::ferx_model_validate(tmp, data = data_file)
    ),
    error = function(e) e
  )
  # The engine could not run. That is a tooling problem, not a defect in the
  # model, so it must not abort the translation or be reported as one -- ferx is
  # a Suggests dependency with no version floor, and blaming the user's .ferx
  # for a signature change would be actively misleading.
  if (inherits(res, "error"))
    return(list(
      ok = NA, diagnostics = .empty_diagnostics(), data_file = data_file,
      warnings = paste0("WARN  | ferx_model_validate() could not run, so the ",
                        "emitted .ferx was NOT validated: ",
                        .one_line(conditionMessage(res))),
      unsupported = character()
    ))

  diags <- out$diagnostics
  if (!is.data.frame(diags)) diags <- .empty_diagnostics()

  # E_DATA means the engine could not read the DATASET -- a file the translator
  # neither produces nor controls, and whose NONMEM $DATA options (IGNORE=, and
  # a headerless file described by $INPUT) it cannot convey. Failing the
  # translation on it would reject models that are perfectly correct, and would
  # make the outcome depend on whether the CSV happens to sit next to the .ctl.
  # Drop back to model-only validation and say so.
  if (any(grepl("^E_DATA", as.character(diags$code)))) {
    msg <- .one_line(diags$message[grepl("^E_DATA", as.character(diags$code))][1])
    res2 <- tryCatch(utils::capture.output(out <- ferx::ferx_model_validate(tmp)),
                     error = function(e) e)
    if (inherits(res2, "error"))
      return(list(ok = NA, diagnostics = .empty_diagnostics(),
                  data_file = NA_character_,
                  warnings = paste0("WARN  | emitted .ferx was NOT validated: ",
                                    .one_line(conditionMessage(res2))),
                  unsupported = character()))
    diags     <- if (is.data.frame(out$diagnostics)) out$diagnostics
                 else .empty_diagnostics()
    data_file <- NA_character_
    data_note <- paste0("WARN  | ferx could not read the dataset, so covariate ",
                        "references and endpoint coverage were NOT checked: ", msg)
  } else data_note <- character()

  warn <- data_note
  if (is.na(data_file) && length(data_note) == 0L)
    warn <- c(warn, paste0(
      "INFO  | validated without data -- covariate references and endpoint ",
      "coverage were NOT checked (an unknown name is read as a covariate)"))
  else
    warn <- c(warn, paste0("INFO  | validated against data: ", data_file))

  if (nrow(diags) > 0) {
    sev  <- as.character(diags$severity)
    warn <- c(warn, vapply(seq_len(nrow(diags)), function(i) paste0(
      if (identical(sev[i], "error")) "ERROR | " else "WARN  | ",
      "ferx ", diags$code[i], ": ", .one_line(diags$message[i])), character(1)))
  }

  is_err <- if (nrow(diags) > 0) as.character(diags$severity) == "error" else logical()
  unsupported <- if (any(is_err))
    paste0("ferx rejected the emitted .ferx (", diags$code[is_err], "): ",
           .one_line(diags$message[is_err]))
  else character()

  list(ok = isTRUE(out$ok), diagnostics = diags, data_file = data_file,
       warnings = warn, unsupported = unsupported)
}

# Emit validation findings to the console at translation time, and stop when
# strict and the engine rejected the output.
.report_validation <- function(val, strict = TRUE, result = NULL) {
  if (is.null(val)) return(invisible(NULL))

  # Engine text is DATA, never a cli format string. cli glue-interpolates its
  # first argument, so a diagnostic or a $DATA path containing braces would be
  # evaluated as R code -- replacing the real message with a cli parse error, or
  # silently rewriting the path. The rest of the package already uses "{w}".
  for (w in val$warnings) {
    if (startsWith(w, "ERROR")) next          # reported together below
    if (startsWith(w, "WARN"))  cli::cli_warn("{w}") else cli::cli_inform("{w}")
  }

  # Gate on the engine's own verdict as well as on error-severity rows: ferx has
  # two independent notions of validity, and a future severity level or an
  # unexpected return shape would otherwise let an invalid model through.
  n_err <- length(val$unsupported)
  if (n_err == 0L && !isFALSE(val$ok)) return(invisible(NULL))

  detail <- if (n_err > 0L) val$unsupported
            else "ferx reported the model as invalid without an error-severity diagnostic"
  msg <- c("The emitted {.field .ferx} is not valid: {length(detail)} problem{?s} from {.pkg ferx}.",
           stats::setNames(.cli_escape(detail), rep("x", length(detail))))
  if (isTRUE(strict))
    # The result rides on the condition so a caller can still reach $ferx_text
    # and $unsupported after catching it, rather than getting an exception and
    # nothing to act on.
    cli::cli_abort(c(msg, i = paste0("Pass {.code strict = FALSE} to return the ",
                                     "result anyway, or {.code validate = FALSE} ",
                                     "to skip validation.")),
                   class = "ferxtranslate_invalid_output", result = result,
                   call = NULL)
  cli::cli_warn(msg)
  invisible(NULL)
}

.empty_diagnostics <- function() {
  data.frame(severity = character(), code = character(), message = character(),
             stringsAsFactors = FALSE)
}

# Neutralise glue braces so engine text can be used where cli expects a format
# string. Doubling is glue's own escape: `{{` renders as a literal `{`.
.cli_escape <- function(x) {
  gsub("}", "}}", gsub("{", "{{", as.character(x), fixed = TRUE), fixed = TRUE)
}

# Collapse a multi-line engine message so one diagnostic stays one warning.
.one_line <- function(x) gsub("\\s+", " ", trimws(as.character(x)))
