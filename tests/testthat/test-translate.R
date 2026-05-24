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
