#' Translate a pharmacometric model to ferx format
#'
#' Parses `source` via the appropriate intermediary (`nonmem2rx`, `rxode2`, or
#' `monolix2rx`), converts to a [new_ferx_ir()] object, emits a `.ferx` string,
#' and optionally writes it to disk.
#'
#' @param source For `"nonmem"`: path to a `.ctl` file. For `"nlmixr2"`: an
#'   nlmixr2/rxode2 model function. For `"monolix"`: path to a `.mlxtran` file.
#' @param format Source format. One of `"nonmem"`, `"nlmixr2"`, `"monolix"`.
#' @param output Optional path for the output `.ferx` file. If `NULL` (default)
#'   the result is returned but not written.
#' @param overwrite If `output` is set and the file already exists, pass `TRUE`
#'   to overwrite it. Default `FALSE`.
#' @param validate Check the emitted `.ferx` with `ferx::ferx_model_validate()`
#'   before returning. Default `TRUE`. Requires the `ferx` package; when it is
#'   not installed the check is skipped and an `INFO` warning says so. For a
#'   NONMEM source the `$DATA` file is resolved relative to the control stream
#'   and passed to the validator when it exists - without a dataset the engine
#'   cannot check covariate references, because every block that accepts a
#'   covariate reads an unknown name as one.
#' @param strict If validation reports any error-severity diagnostic, abort.
#'   Default `TRUE`. Pass `FALSE` to downgrade to a warning and return the
#'   result anyway. Ignored when `validate = FALSE`.
#'
#' @return A `ferx_translate_result` invisibly, with fields `$ferx_text`,
#'   `$warnings`, `$unsupported`, and `$validation` (`NULL` when
#'   `validate = FALSE`, otherwise a list with `ok`, `diagnostics`, the
#'   `data_file` used, and the full `warnings` including the `INFO` notes on how
#'   validation ran). Only engine findings about the model reach `$warnings` and
#'   the emitted file; notes about validation coverage stay in `$validation`.
#'
#' @seealso [write_ferx()], [rxui_to_ir()], [emit_ferx()]
#'
#' @examples
#' \dontrun{
#' # From NONMEM
#' result <- nm_to_ferx("run001.ctl", output = "run001.ferx")
#' result$warnings
#'
#' # From nlmixr2
#' f <- function() {
#'   ini({ tvcl <- 0.134; eta.cl ~ 0.07; err.prop ~ 0.01 })
#'   model({ cl <- tvcl * exp(eta.cl); linCmt() ~ prop(err.prop) })
#' }
#' result <- nlmixr2_to_ferx(f)
#' cat(result$ferx_text)
#'
#' # From Monolix
#' result <- mlx_to_ferx("project.mlxtran")
#' }
#' @export
to_ferx <- function(source,
                    format    = c("nonmem", "nlmixr2", "monolix"),
                    output    = NULL,
                    overwrite = FALSE,
                    validate  = TRUE,
                    strict    = TRUE) {
  format <- match.arg(format)

  if (format == "monolix" && !requireNamespace("monolix2rx", quietly = TRUE))
    cli::cli_abort(c(
      "monolix2rx is required for Monolix translation.",
      i = "Install it with: {.code install.packages('monolix2rx')}"
    ))

  src_label <- if (is.character(source)) source else "<model function>"
  rxui <- tryCatch(
    switch(format,
      nonmem  = nonmem2rx::nonmem2rx(source),
      nlmixr2 = rxode2::rxode2(source),
      monolix = monolix2rx::monolix2rx(source)
    ),
    error = function(e)
      cli::cli_abort(
        c(paste0("Failed to parse ", src_label, " as ", format, "."),
          "i" = paste0("Original error: ", conditionMessage(e))),
        call = NULL
      )
  )

  src_file     <- if (is.character(source)) source else NA_character_
  scaling_hint <- if (format == "nonmem" && is.character(source) && file.exists(source))
    .extract_nm_scaling(source) else list()
  ir <- rxui_to_ir(rxui, source_format = format, source_file = src_file,
                   scaling_hint = scaling_hint)
  text <- emit_ferx(ir)

  val <- if (isTRUE(validate)) {
    data_file <- if (format == "nonmem" && is.character(source))
      .extract_nm_data_path(source) else NA_character_
    .validate_ferx_text(text, data_file = data_file)
  } else NULL

  if (!is.null(val)) {
    # Only findings ABOUT THE MODEL go into the file. The INFO notes about how
    # validation ran (skipped, with or without data) belong to the translation
    # session, not the artefact -- folding them in would make every clean model
    # emit "# Warnings: 1" for a note that says nothing is wrong.
    findings <- val$warnings[!startsWith(val$warnings, "INFO")]
    ir$warnings <- c(ir$warnings, findings)
    ir$unsupported <- c(ir$unsupported, val$unsupported)
    # Re-emit only when something was added, so a clean translation stays
    # byte-identical. The emitter is pure, so this cannot change the model.
    if (length(findings) > 0 || length(val$unsupported) > 0)
      text <- emit_ferx(ir)
  }

  result <- new_ferx_translate_result(text, ir, validation = val)

  # Surface diagnostics before writing, so a strict failure cannot leave an
  # invalid file on disk that the user believes was checked.
  .report_validation(val, strict = strict)

  if (!is.null(output)) write_ferx(result, output, overwrite)
  invisible(result)
}

#' @rdname to_ferx
#' @param ctl_file Path to a NONMEM control stream file.
#' @param ... Additional arguments passed to [to_ferx()] (e.g. `overwrite`).
#' @export
nm_to_ferx <- function(ctl_file, output = NULL, ...) {
  to_ferx(ctl_file, "nonmem", output, ...)
}

#' @rdname to_ferx
#' @param model_fn An nlmixr2 / rxode2 model function.
#' @export
nlmixr2_to_ferx <- function(model_fn, output = NULL, ...) {
  to_ferx(model_fn, "nlmixr2", output, ...)
}

#' @rdname to_ferx
#' @param mlxtran Path to a Monolix `.mlxtran` project file.
#' @export
mlx_to_ferx <- function(mlxtran, output = NULL, ...) {
  to_ferx(mlxtran, "monolix", output, ...)
}

# -- ferx_translate_result ----------------------------------------------------

#' Create a ferx_translate_result object
#'
#' @param text Character string: the complete `.ferx` file content.
#' @param ir A `ferx_ir` object from which `warnings` and `unsupported` are
#'   pulled.
#' @param validation Result of `.validate_ferx_text()`, or `NULL` when
#'   validation was not run.
#'
#' @return A `ferx_translate_result` list with fields `ferx_text`, `warnings`,
#'   `unsupported`, `validation`, `source_format`, and `source_file`.
#'
#' @seealso [to_ferx()], [write_ferx()]
#' @noRd
new_ferx_translate_result <- function(text, ir, validation = NULL) {
  structure(
    list(
      ferx_text     = text,
      warnings      = ir$warnings,
      unsupported   = ir$unsupported,
      validation    = validation,
      source_format = ir$source_format,
      source_file   = ir$source_file
    ),
    class = "ferx_translate_result"
  )
}

#' @export
print.ferx_translate_result <- function(x, ...) {
  src  <- if (!is.na(x$source_format)) x$source_format else "unknown"
  file <- if (!is.na(x$source_file))   x$source_file   else "unknown"
  cli::cli_h1("Translated from {src}: {file}")

  cat(x$ferx_text, "\n", sep = "")

  if (length(x$warnings) > 0) {
    cli::cli_h2("Translation warnings")
    for (w in x$warnings) cli::cli_alert_warning("{w}")
  }

  if (length(x$unsupported) > 0) {
    cli::cli_h2("Unsupported features (manual fix required)")
    for (u in x$unsupported) cli::cli_alert_danger("{u}")
  }

  cli::cli_h2("Engine validation")
  if (is.null(x$validation)) {
    cli::cli_alert_info("not run ({.code validate = FALSE})")
  } else if (is.na(x$validation$ok)) {
    cli::cli_alert_info("skipped -- {.pkg ferx} is not installed")
  } else {
    scope <- if (is.na(x$validation$data_file)) "model only, no data"
             else paste0("with data: ", x$validation$data_file)
    if (isTRUE(x$validation$ok)) cli::cli_alert_success("valid ({scope})")
    else                         cli::cli_alert_danger("INVALID ({scope})")
  }

  invisible(x)
}

# -- write_ferx ---------------------------------------------------------------

#' Write a ferx_translate_result to a .ferx file
#'
#' @param result A `ferx_translate_result`.
#' @param path Output file path.
#' @param overwrite If the file exists and this is `FALSE`, an error is raised.
#'
#' @return `result` invisibly, for pipe use.
#'
#' @seealso [to_ferx()]
#'
#' @examples
#' \dontrun{
#' result <- nm_to_ferx("run001.ctl")
#' write_ferx(result, "run001.ferx")
#' # or via pipe:
#' nm_to_ferx("run001.ctl") |> write_ferx("run001.ferx")
#' }
#' @export
write_ferx <- function(result, path, overwrite = FALSE) {
  if (!overwrite && file.exists(path))
    cli::cli_abort(c(
      "{.path {path}} already exists.",
      i = "Pass {.code overwrite = TRUE} to replace it."
    ))
  writeLines(result$ferx_text, path)
  invisible(result)
}
