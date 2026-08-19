# -- helpers ------------------------------------------------------------------

minimal_ir <- function(...) {
  new_ferx_ir(
    source_format = "nonmem",
    source_file   = "run001.ctl",
    thetas        = list(list(name = "TVCL", init = 0.134, lower = 0.001, upper = 10)),
    omegas        = list(list(type = "diagonal", names = "ETA_CL", values = 0.07)),
    sigmas        = list(list(name = "PROP_ERR", value = 0.01, scale = "sd")),
    structural    = list(type = "pk_macro", pk_call = "one_cpt_oral",
                         pk_args = list(cl = "CL", v = "V", ka = "KA")),
    error_model   = list(list(dv = "DV", type = "proportional", params = "PROP_ERR")),
    fit_options   = list(method = "foce", maxiter = 300L, covariance = TRUE),
    ...
  )
}

# -- new_ferx_translate_result ------------------------------------------------

test_that("result stores ferx_text, warnings, unsupported from ir", {
  ir  <- minimal_ir(warnings = c("WARN | check"), unsupported = c("MIXTURE"))
  res <- new_ferx_translate_result(emit_ferx(ir), ir)

  expect_s3_class(res, "ferx_translate_result")
  expect_type(res$ferx_text, "character")
  expect_match(res$ferx_text, "[parameters]", fixed = TRUE)
  expect_equal(res$warnings,      c("WARN | check"))
  expect_equal(res$unsupported,   c("MIXTURE"))
  expect_equal(res$source_format, "nonmem")
  expect_equal(res$source_file,   "run001.ctl")
})

test_that("result with no warnings or unsupported has empty vectors", {
  ir  <- minimal_ir()
  res <- new_ferx_translate_result(emit_ferx(ir), ir)
  expect_equal(res$warnings,    character())
  expect_equal(res$unsupported, character())
})

# -- print.ferx_translate_result ----------------------------------------------

test_that("print runs without error for clean result", {
  ir  <- minimal_ir()
  res <- new_ferx_translate_result(emit_ferx(ir), ir)
  expect_no_error(print(res))
})

test_that("print includes ferx_text in output", {
  ir  <- minimal_ir()
  res <- new_ferx_translate_result(emit_ferx(ir), ir)
  out <- capture.output(print(res))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "[parameters]",    fixed = TRUE)
  expect_match(combined, "one_cpt_oral",    fixed = TRUE)
})

test_that("print with warnings runs without error and object holds warnings", {
  ir  <- minimal_ir(warnings = c("WARN | check theta bounds"))
  res <- new_ferx_translate_result(emit_ferx(ir), ir)
  expect_no_error(print(res))
  expect_equal(res$warnings, c("WARN | check theta bounds"))
})

test_that("print with unsupported runs without error and object holds unsupported", {
  ir  <- minimal_ir(unsupported = c("MIXTURE model"))
  res <- new_ferx_translate_result(emit_ferx(ir), ir)
  expect_no_error(print(res))
  expect_equal(res$unsupported, c("MIXTURE model"))
})

test_that("result with no warnings/unsupported has empty vectors", {
  ir  <- minimal_ir()
  res <- new_ferx_translate_result(emit_ferx(ir), ir)
  expect_equal(res$warnings,    character())
  expect_equal(res$unsupported, character())
})

test_that("print returns result invisibly", {
  ir  <- minimal_ir()
  res <- new_ferx_translate_result(emit_ferx(ir), ir)
  expect_invisible(print(res))
})

test_that("print handles unknown source gracefully", {
  ir  <- new_ferx_ir()
  res <- new_ferx_translate_result("", ir)
  expect_no_error(print(res))
})

# -- write_ferx ---------------------------------------------------------------

test_that("write_ferx creates file with ferx_text content", {
  ir   <- minimal_ir()
  res  <- new_ferx_translate_result(emit_ferx(ir), ir)
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  write_ferx(res, path)

  expect_true(file.exists(path))
  written <- paste(readLines(path), collapse = "\n")
  expect_match(written, "[parameters]",   fixed = TRUE)
  expect_match(written, "one_cpt_oral",   fixed = TRUE)
})

test_that("write_ferx returns result invisibly", {
  ir   <- minimal_ir()
  res  <- new_ferx_translate_result(emit_ferx(ir), ir)
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))
  expect_invisible(write_ferx(res, path))
})

test_that("write_ferx errors when file exists and overwrite = FALSE", {
  ir   <- minimal_ir()
  res  <- new_ferx_translate_result(emit_ferx(ir), ir)
  path <- tempfile(fileext = ".ferx")
  writeLines("existing content", path)
  on.exit(unlink(path))

  expect_error(write_ferx(res, path, overwrite = FALSE), "already exists")
})

test_that("write_ferx overwrites when overwrite = TRUE", {
  ir   <- minimal_ir()
  res  <- new_ferx_translate_result(emit_ferx(ir), ir)
  path <- tempfile(fileext = ".ferx")
  writeLines("old content", path)
  on.exit(unlink(path))

  write_ferx(res, path, overwrite = TRUE)
  written <- readLines(path)
  expect_false(any(grepl("old content", written, fixed = TRUE)))
  expect_true(any(grepl("[parameters]", written, fixed = TRUE)))
})

test_that("write_ferx pipe: result unchanged after write", {
  ir   <- minimal_ir()
  res  <- new_ferx_translate_result(emit_ferx(ir), ir)
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  returned <- write_ferx(res, path)
  expect_identical(returned, res)
})

# -- to_ferx (requires nonmem2rx / rxode2) ------------------------------------

test_that("to_ferx NONMEM errors without nonmem2rx", {
  skip_if_not(
    !requireNamespace("nonmem2rx", quietly = TRUE),
    "nonmem2rx is installed -- skipping absent-package check"
  )
  expect_error(to_ferx("run001.ctl", "nonmem"))
})

test_that("to_ferx nlmixr2 round-trip", {
  skip_if_not_installed("rxode2")
  f_1cpt <- function() {
    ini({
      tvcl <- 0.134; tvv <- 8.1; tvka <- 1.0
      eta.cl ~ 0.07; eta.v ~ 0.02
      prop.err <- 0.01
    })
    model({
      cl <- tvcl * exp(eta.cl)
      v  <- tvv  * exp(eta.v)
      ka <- tvka
      linCmt() ~ prop(prop.err)
    })
  }

  res <- to_ferx(f_1cpt, "nlmixr2")

  expect_s3_class(res, "ferx_translate_result")
  expect_match(res$ferx_text, "[parameters]",   fixed = TRUE)
  expect_match(res$ferx_text, "one_cpt_oral",   fixed = TRUE)
  expect_match(res$ferx_text, "[error_model]",  fixed = TRUE)
})

test_that("to_ferx writes file when output is given", {
  skip_if_not_installed("rxode2")
  f_1cpt <- function() {
    ini({ tvcl <- 0.134; tvv <- 8.1; eta.cl ~ 0.07; prop.err <- 0.01 })
    model({ cl <- tvcl * exp(eta.cl); v <- tvv; linCmt() ~ prop(prop.err) })
  }
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  to_ferx(f_1cpt, "nlmixr2", output = path)
  expect_true(file.exists(path))
})

test_that("mlx_to_ferx errors cleanly without monolix2rx", {
  skip_if(requireNamespace("monolix2rx", quietly = TRUE),
          "monolix2rx installed")
  expect_error(mlx_to_ferx("project.mlxtran"), "monolix2rx")
})

# -- $DATA path resolution ----------------------------------------------------

test_that("$DATA path resolves relative to the control stream", {
  dir  <- tmp_ctl_dir()
  data <- file.path(dir, "mydata.csv")
  writeLines("ID,TIME,DV", data)
  ctl <- file.path(dir, "run1.ctl")
  writeLines(c("$PROBLEM x", "$DATA mydata.csv IGNORE=@", "$PK"), ctl)
  expect_equal(normalizePath(.extract_nm_data_path(ctl)), normalizePath(data))
})

test_that("$DATA accepts a quoted path containing spaces", {
  dir  <- tmp_ctl_dir()
  data <- file.path(dir, "my data.csv")
  writeLines("ID,TIME,DV", data)
  ctl <- file.path(dir, "run1.ctl")
  writeLines(c("$PROBLEM x", '$DATA "my data.csv" IGNORE=@'), ctl)
  expect_equal(normalizePath(.extract_nm_data_path(ctl)), normalizePath(data))
})

test_that("$DATA returns NA when the dataset is absent or undeclared", {
  dir <- tmp_ctl_dir()
  no_file <- file.path(dir, "a.ctl")
  writeLines(c("$PROBLEM x", "$DATA notthere.csv IGNORE=@"), no_file)
  expect_true(is.na(.extract_nm_data_path(no_file)))

  no_data <- file.path(dir, "b.ctl")
  writeLines(c("$PROBLEM x", "$PK"), no_data)
  expect_true(is.na(.extract_nm_data_path(no_data)))

  expect_true(is.na(.extract_nm_data_path("/nonexistent/path.ctl")))
})

test_that("$DATA ignores a trailing comment", {
  dir  <- tmp_ctl_dir()
  data <- file.path(dir, "d.csv"); writeLines("ID", data)
  ctl  <- file.path(dir, "c.ctl")
  writeLines(c("$PROBLEM x", "$DATA d.csv ; the dataset"), ctl)
  expect_equal(normalizePath(.extract_nm_data_path(ctl)), normalizePath(data))
})

# -- output validation --------------------------------------------------------

test_that("validate = FALSE skips validation and reports nothing", {
  skip_if_not_installed("nonmem2rx")
  result <- nm_to_ferx(nm_path("1cpt_oral.ctl"), validate = FALSE)
  expect_null(result$validation)
  expect_false(any(grepl("validated", result$warnings)))
})

test_that("a valid model validates clean and records the outcome", {
  skip_if_not_installed("nonmem2rx")
  skip_if_not_installed("ferx")
  result <- suppressMessages(nm_to_ferx(nm_path("1cpt_oral.ctl")))
  expect_true(result$validation$ok)
  expect_equal(nrow(result$validation$diagnostics), 0L)
  expect_length(result$unsupported, 0L)
})

test_that("validating without a dataset says so rather than implying a clean bill", {
  skip_if_not_installed("nonmem2rx")
  skip_if_not_installed("ferx")
  # The bundled .ctl files name datasets that are not shipped, so this is the
  # no-data path. It matters because an unknown name is read as a covariate in
  # every block that accepts one, and only data turns that into an error.
  result <- suppressMessages(nm_to_ferx(nm_path("1cpt_oral.ctl")))
  expect_true(is.na(result$validation$data_file))
  expect_true(any(grepl("^INFO.*without data", result$validation$warnings)))
  # The note is about the session, not the model, so it must not leak into the
  # emitted file and make a clean translation advertise a warning.
  expect_false(any(grepl("without data", result$warnings)))
  expect_false(grepl("# Warnings:", result$ferx_text))
})

test_that("invalid output aborts under strict and warns otherwise", {
  skip_if_not_installed("ferx")
  # An INCOMPLETE model, not a misspelt one. This used to use a dotted state
  # name (`c.BAD`) and stopped working the moment validate_ferx_ir() began
  # rejecting illegal declared identifiers -- the same rot the test below this
  # one already records for `c.RTOT`. The subject here is the strict/non-strict
  # branch of .report_validation(), so the vehicle must be something no phase of
  # the translator will ever make valid: an IR that declares an ODE and simply
  # has no [individual_parameters] and no [error_model] is not a defect being
  # exercised, it is a model that is genuinely missing required blocks, and
  # ferx answers E_MISSING_BLOCK. Every name in it is a legal ferx identifier
  # on purpose, so the fixture does not depend on a name rule holding still.
  bad <- new_ferx_ir(
    thetas     = list(list(name = "TVCL", init = 1, lower = 0, upper = 10)),
    structural = list(type = "ode", obs_cmt = "CENT", states = "CENT"),
    odes       = list(list(state = "CENT", rhs = "-TVCL * CENT"))
  )
  val <- .validate_ferx_text(emit_ferx(bad))
  expect_false(val$ok)
  expect_gt(length(val$unsupported), 0L)
  expect_error(.report_validation(val, strict = TRUE), "not valid")
  expect_warning(.report_validation(val, strict = FALSE), "not valid")
})

test_that("a missing ferx skips validation without failing the translation", {
  skip_if_not_installed("nonmem2rx")
  # This is how the fast PR job runs: ferx is a Suggests dependency and only the
  # engine job installs it. An optional dependency must never turn a working
  # translation into an error, so assert the degradation rather than assume it.
  local_mocked_bindings(.has_ferx = function() FALSE)
  result <- expect_no_error(
    suppressMessages(nm_to_ferx(nm_path("1cpt_oral.ctl")))
  )
  expect_true(is.na(result$validation$ok))
  expect_match(result$validation$warnings, "ferx is not installed", all = FALSE)
  expect_length(result$unsupported, 0L)
  # Nothing about the skipped check may leak into the emitted file.
  expect_false(grepl("# Warnings:", result$ferx_text))
})

test_that("a strict abort leaves no output file behind", {
  skip_if_not_installed("ferx")
  skip_if_not_installed("nonmem2rx")
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "bad.ctl")
  # $DES references a name that is defined nowhere -- not a theta, not a state,
  # not a $PK variable, not a data column. `[odes]` rejects undefined names
  # outright (it does not even fall back to reading them as covariates), so this
  # is invalid for a reason no phase of the translator can remove: the SOURCE is
  # broken, there is nothing to translate correctly.
  #
  # It used to use a dotted state name (`c.RTOT`, issue #6 defect 1), which
  # stopped failing the moment that defect was fixed. The subject under test is
  # the validate-before-write ordering, so it must not rest on a live defect.
  writeLines(c(
    "$PROBLEM undefined name in DES", "$INPUT ID TIME AMT DV MDV",
    "$DATA d.csv IGNORE=@", "$SUBROUTINES ADVAN13 TOL=9",
    "$MODEL", "  COMP=(CENT, DEFDOSE, DEFOBS)",
    "$PK", "  KEL = THETA(1)*EXP(ETA(1))", "  VC = THETA(2)",
    "$DES", "  DADT(1) = -KEL*A(1)*NOSUCHNAME",
    "$ERROR", "  Y = A(1)/VC*(1+EPS(1))",
    "$THETA (0,0.05) (0,3)", "$OMEGA 0.09", "$SIGMA 0.04", "$EST METHOD=1"), ctl)
  writeLines("ID,TIME,AMT,DV,MDV\n1,0,100,0,1", file.path(dir, "d.csv"))
  out <- file.path(dir, "out.ferx")

  expect_error(suppressMessages(nm_to_ferx(ctl, output = out)), "not valid")
  expect_false(file.exists(out))

  # strict = FALSE still writes, but says what is wrong in the file itself.
  res <- suppressWarnings(suppressMessages(
    nm_to_ferx(ctl, output = out, strict = FALSE)))
  expect_true(file.exists(out))
  expect_false(res$validation$ok)
  # Engine rejections are carried as ERROR-prefixed warnings, not as
  # $unsupported (which is the ferx-core feature-gap signal), and the emitter
  # renders every WARN/ERROR into the file's # WARNING: block.
  expect_match(res$ferx_text, "# WARNING: ferx E_PARSE")
  expect_length(res$unsupported, 0L)
})

# -- behaviours a mutation campaign found unasserted ---------------------------

test_that("the resolved $DATA path actually reaches the validator", {
  skip_if_not_installed("ferx")
  skip_if_not_installed("nonmem2rx")
  # Stubbing data_file to NA in to_ferx() previously left the whole suite green:
  # the $DATA tests exercised .extract_nm_data_path() in isolation and nothing
  # asserted the wiring, so the documented reason for resolving it at all
  # (covariate references cannot be checked without data) was untested.
  dir <- tmp_ctl_dir()
  file.copy(nm_path("1cpt_oral.ctl"), file.path(dir, "m.ctl"))
  writeLines(c("ID,TIME,AMT,DV,MDV,EVID,CMT", "1,0,100,.,1,1,1"),
             file.path(dir, "1cpt_oral.csv"))
  res <- suppressWarnings(suppressMessages(
    nm_to_ferx(file.path(dir, "m.ctl"), strict = FALSE)))
  expect_false(is.na(res$validation$data_file))
  expect_match(res$validation$data_file, "1cpt_oral\\.csv$")

  # $validation$data_file alone proves nothing: .validate_ferx_text() echoes
  # back whatever path it was handed, so it reads the same whether or not the
  # engine ever saw it. Assert on a diagnostic only a data-aware run can
  # produce -- a covariate reference is invisible without a dataset.
  ir <- suppressWarnings(rxui_to_ir(nonmem2rx::nonmem2rx(nm_path("1cpt_oral.ctl")),
                                    source_format = "nonmem"))
  ir$indiv_params <- c(ir$indiv_params,
                       list(list(lhs = "ZZ", rhs = "NO_SUCH_COLUMN * 2")))
  txt <- emit_ferx(ir)
  with_data <- .validate_ferx_text(txt, data_file = file.path(dir, "1cpt_oral.csv"))
  no_data   <- .validate_ferx_text(txt, data_file = NA_character_)
  expect_true(any(grepl("NO_SUCH_COLUMN", with_data$warnings)))
  expect_false(any(grepl("NO_SUCH_COLUMN", no_data$warnings)))
})

test_that("engine severities route to the right console channel", {
  # Swapping cli_warn and cli_inform in .report_validation() previously changed
  # nothing any test could see.
  val <- list(ok = FALSE, diagnostics = .empty_diagnostics(),
              data_file = NA_character_,
              warnings = c("INFO  | a note", "WARN  | an engine warning"),
              unsupported = character())
  expect_message(suppressWarnings(.report_validation(val, strict = FALSE)), "a note")
  expect_warning(suppressMessages(.report_validation(val, strict = FALSE)),
                 "an engine warning")
})

test_that("engine text with braces is never evaluated as cli syntax", {
  # cli glue-interpolates its first argument, so raw engine text used as a
  # format string would execute R code or die with a cli parse error.
  val <- list(ok = FALSE, diagnostics = .empty_diagnostics(),
              data_file = NA_character_,
              warnings = "WARN  | value must be in {0, 1}",
              unsupported = 'rejected: {Sys.getenv("USER")}')
  w <- tryCatch(suppressMessages(.report_validation(val, strict = FALSE)),
                warning = function(w) conditionMessage(w))
  expect_match(w, "value must be in {0, 1}", fixed = TRUE)
  err <- tryCatch(.report_validation(val, strict = TRUE), error = function(e) e)
  expect_match(conditionMessage(err), 'Sys.getenv("USER")', fixed = TRUE)
})

test_that("a strict abort carries the result on the condition", {
  val <- list(ok = FALSE, diagnostics = .empty_diagnostics(),
              data_file = NA_character_, warnings = character(),
              unsupported = "rejected")
  err <- tryCatch(.report_validation(val, strict = TRUE, result = "THE-RESULT"),
                  error = function(e) e)
  expect_s3_class(err, "ferxtranslate_invalid_output")
  expect_equal(err$result, "THE-RESULT")
})

test_that("an invalid verdict aborts even with no error-severity diagnostic", {
  # ferx has two notions of validity; gating only on severity == "error" let an
  # ok = FALSE model be written to disk as a success.
  val <- list(ok = FALSE, diagnostics = .empty_diagnostics(),
              data_file = NA_character_, warnings = character(),
              unsupported = character())
  expect_error(.report_validation(val, strict = TRUE), "not valid")
})

test_that("a data-read failure degrades to model-only instead of failing", {
  skip_if_not_installed("ferx")
  skip_if_not_installed("nonmem2rx")
  # E_DATA is about the dataset, which the translator neither produces nor
  # controls, and whose $DATA IGNORE= options it cannot convey to the engine.
  # A model that is itself fine must not be rejected because of it.
  dir <- tmp_ctl_dir()
  file.copy(nm_path("1cpt_oral.ctl"), file.path(dir, "m.ctl"))
  writeLines(c("# a comment line NONMEM would IGNORE", "nonsense,columns"),
             file.path(dir, "1cpt_oral.csv"))
  res <- expect_no_error(suppressWarnings(suppressMessages(
    nm_to_ferx(file.path(dir, "m.ctl")))))
  expect_true(res$validation$ok)
  expect_true(any(grepl("could not read the dataset", res$validation$warnings)))
})

test_that("a directory is never offered to the engine as a dataset", {
  dir <- tmp_ctl_dir()
  dir.create(file.path(dir, "dat"))
  ctl <- file.path(dir, "m.ctl")
  writeLines(c("$PROBLEM x", "$DATA dat IGNORE=#"), ctl)
  expect_true(is.na(.extract_nm_data_path(ctl)))
})

test_that("an absolute $DATA path is used as given", {
  dir  <- tmp_ctl_dir()
  data <- file.path(dir, "abs.csv"); writeLines("ID", data)
  ctl  <- file.path(dir, "m.ctl")
  writeLines(c("$PROBLEM x", paste("$DATA", normalizePath(data))), ctl)
  expect_equal(normalizePath(.extract_nm_data_path(ctl)), normalizePath(data))
  expect_true(.is_abs_path("/a/b.csv"))
  expect_false(.is_abs_path("a/b.csv"))
})

test_that("$DATA strips a comment with no preceding space", {
  dir  <- tmp_ctl_dir()
  data <- file.path(dir, "d.csv"); writeLines("ID", data)
  ctl  <- file.path(dir, "c.ctl")
  # The earlier test used "d.csv ; comment", which the first-token split already
  # handles -- so it could not detect the comment strip being removed.
  writeLines(c("$PROBLEM x", "$DATA d.csv;the dataset"), ctl)
  expect_equal(normalizePath(.extract_nm_data_path(ctl)), normalizePath(data))
})

test_that("a multi-line engine message collapses to one line", {
  expect_equal(.one_line("a\n  b\tc "), "a b c")
})

test_that("print() reports each engine-validation outcome distinctly", {
  ir <- new_ferx_ir(structural = list(type = "pk_macro", pk_call = "one_cpt_iv",
                                      pk_args = list(cl = "CL", v = "V")))
  # cli alerts go to stderr, so capture the message stream, not stdout.
  shown <- function(res) paste(utils::capture.output(print(res), type = "message"),
                              collapse = "\n")
  expect_match(shown(new_ferx_translate_result("x", ir)), "not run")
  # ok = NA covers two outcomes; with ferx installed it means the engine failed,
  # not that it is absent.
  expect_match(shown(new_ferx_translate_result("x", ir,
    validation = list(ok = NA, data_file = NA_character_))),
    if (.has_ferx()) "could not run" else "not installed")
  expect_match(shown(new_ferx_translate_result("x", ir,
    validation = list(ok = TRUE, data_file = NA_character_))), "valid")
  expect_match(shown(new_ferx_translate_result("x", ir,
    validation = list(ok = FALSE, data_file = NA_character_))), "INVALID")
})

test_that("a data-read failure leaves the emitted .ferx machine-independent", {
  skip_if_not_installed("ferx")
  skip_if_not_installed("nonmem2rx")
  # The note is INFO, not WARN: it describes how validation ran, not a defect in
  # the model. A WARN would be written into the file, so the same control stream
  # would emit different text depending on whether a CSV happened to sit beside
  # it -- complete with an embedded absolute path.
  dir <- tmp_ctl_dir()
  file.copy(nm_path("1cpt_oral.ctl"), file.path(dir, "m.ctl"))
  bare <- suppressWarnings(suppressMessages(nm_to_ferx(file.path(dir, "m.ctl"))))
  writeLines(c("# comment NONMEM would IGNORE", "nonsense,columns"),
             file.path(dir, "1cpt_oral.csv"))
  with_csv <- suppressWarnings(suppressMessages(nm_to_ferx(file.path(dir, "m.ctl"))))

  expect_equal(with_csv$ferx_text, bare$ferx_text)
  expect_false(grepl("could not read the dataset", with_csv$ferx_text))
  # ...and the note must not contradict itself by also claiming data coverage.
  expect_true(any(grepl("could not read the dataset", with_csv$validation$warnings)))
  expect_false(any(grepl("validated against data", with_csv$validation$warnings)))
})

test_that("a data error caused by the model is not absorbed as a file problem", {
  skip_if_not_installed("ferx")
  skip_if_not_installed("nonmem2rx")
  # E_DATA is ferx-core's catch-all for any read failure, including ones the
  # translator causes: an emitted iov_column naming a column the dataset lacks
  # aborts the read before the engine's own E_IOV_MISSING_OCC can fire. Absorbing
  # that would return a clean bill of health for a model that cannot be fit.
  ir <- suppressWarnings(rxui_to_ir(nonmem2rx::nonmem2rx(nm_path("1cpt_oral.ctl")),
                                    source_format = "nonmem"))
  ir$fit_options$iov_column <- "OCC"     # no OCC column in the dataset
  data <- system.file("testdata", "1cpt_oral_concordance.csv",
                      package = "ferxtranslate", mustWork = TRUE)
  val <- .validate_ferx_text(emit_ferx(ir), data_file = data)

  expect_false(isTRUE(val$ok))
  expect_gt(length(val$unsupported), 0L)
  expect_match(val$unsupported, "iov_column", all = FALSE)
  # ...and it must not simultaneously claim the run was data-backed.
  expect_false(any(grepl("validated against data", val$warnings)))
})
