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
  if (inherits(res, "error"))
    return(list(
      ok = FALSE, diagnostics = .empty_diagnostics(), data_file = data_file,
      warnings = paste0("ERROR | ferx_model_validate() failed on the emitted ",
                        ".ferx: ", conditionMessage(res)),
      unsupported = paste0("emitted .ferx could not be validated: ",
                           conditionMessage(res))
    ))

  diags <- out$diagnostics
  if (!is.data.frame(diags)) diags <- .empty_diagnostics()

  warn <- character()
  if (is.na(data_file))
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
.report_validation <- function(val, strict = TRUE) {
  if (is.null(val)) return(invisible(NULL))

  n_err <- length(val$unsupported)
  for (w in val$warnings) {
    if (startsWith(w, "ERROR")) next          # reported together below
    if (startsWith(w, "WARN"))  cli::cli_warn(w) else cli::cli_inform(w)
  }
  if (n_err == 0L) return(invisible(NULL))

  msg <- c(
    "The emitted {.field .ferx} is not valid: {n_err} error{?s} from {.pkg ferx}.",
    stats::setNames(val$unsupported, rep("x", n_err))
  )
  if (isTRUE(strict))
    cli::cli_abort(c(msg,
      i = "Pass {.code strict = FALSE} to return the result anyway, or {.code validate = FALSE} to skip validation."),
      call = NULL)
  cli::cli_warn(msg)
  invisible(NULL)
}

.empty_diagnostics <- function() {
  data.frame(severity = character(), code = character(), message = character(),
             stringsAsFactors = FALSE)
}

# Collapse a multi-line engine message so one diagnostic stays one warning.
.one_line <- function(x) gsub("\\s+", " ", trimws(as.character(x)))
