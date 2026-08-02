# -- helpers ------------------------------------------------------------------

nm_path_t <- function(file) {
  system.file("testmodels", "nonmem", file, package = "ferxtranslate",
              mustWork = TRUE)
}

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

# Base-R scratch directory; withr is not a declared dependency of this package.
# R clears tempdir() at session end, so no explicit cleanup is needed.
tmp_ctl_dir <- function() {
  d <- tempfile("ferxtr-")
  dir.create(d)
  d
}

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
  result <- nm_to_ferx(nm_path_t("1cpt_oral.ctl"), validate = FALSE)
  expect_null(result$validation)
  expect_false(any(grepl("validated", result$warnings)))
})

test_that("a valid model validates clean and records the outcome", {
  skip_if_not_installed("nonmem2rx")
  skip_if_not_installed("ferx")
  result <- suppressMessages(nm_to_ferx(nm_path_t("1cpt_oral.ctl")))
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
  result <- suppressMessages(nm_to_ferx(nm_path_t("1cpt_oral.ctl")))
  expect_true(is.na(result$validation$data_file))
  expect_true(any(grepl("^INFO.*without data", result$validation$warnings)))
  # The note is about the session, not the model, so it must not leak into the
  # emitted file and make a clean translation advertise a warning.
  expect_false(any(grepl("without data", result$warnings)))
  expect_false(grepl("# Warnings:", result$ferx_text))
})

test_that("invalid output aborts under strict and warns otherwise", {
  skip_if_not_installed("ferx")
  bad <- new_ferx_ir(
    thetas     = list(list(name = "TVCL", init = 1, lower = 0, upper = 10)),
    structural = list(type = "ode", obs_cmt = "c.BAD", states = "c.BAD"),
    odes       = list(list(state = "c.BAD", rhs = "-TVCL * c.BAD"))
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
    suppressMessages(nm_to_ferx(nm_path_t("1cpt_oral.ctl")))
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
  # $MODEL names a compartment that collides with an $ERROR variable, so
  # nonmem2rx prefixes it to `c.RTOT` and the emitted state identifier carries a
  # dot that ferx cannot parse (issue #6, defect 1).
  writeLines(c(
    "$PROBLEM dotted state", "$INPUT ID TIME AMT DV MDV",
    "$DATA d.csv IGNORE=@", "$SUBROUTINES ADVAN13 TOL=9",
    "$MODEL", "  COMP=(CENT, DEFDOSE, DEFOBS)", "  COMP=(RTOT)",
    "$PK", "  KEL = THETA(1)*EXP(ETA(1))", "  VC = THETA(2)",
    "$DES", "  DADT(1) = -KEL*A(1)", "  DADT(2) = -KEL*A(2)",
    "$ERROR", "  RTOT = A(2)", "  Y = A(1)/VC*(1+EPS(1))",
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
  expect_match(res$ferx_text, "# WARNING: ferx rejected")
})
