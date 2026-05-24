test_that("new_ferx_ir constructs with defaults", {
  ir <- new_ferx_ir()
  expect_s3_class(ir, "ferx_ir")
  expect_true(is.na(ir$source_format))
  expect_true(is.na(ir$source_file))
  expect_equal(ir$thetas, list())
  expect_equal(ir$omegas, list())
  expect_equal(ir$warnings, character())
})

test_that("new_ferx_ir stores supplied values", {
  ir <- new_ferx_ir(
    source_format = "nonmem",
    source_file   = "run001.ctl",
    thetas        = list(list(name = "TVCL", init = 0.134, lower = 0.001, upper = 10)),
    omegas        = list(list(type = "diagonal", names = "ETA_CL", values = 0.07)),
    warnings      = c("WARN | something")
  )
  expect_equal(ir$source_format, "nonmem")
  expect_equal(ir$thetas[[1]]$name, "TVCL")
  expect_equal(ir$omegas[[1]]$type, "diagonal")
  expect_length(ir$warnings, 1L)
})

test_that("validate_ferx_ir accepts empty IR", {
  ir <- new_ferx_ir()
  expect_invisible(validate_ferx_ir(ir))
})

test_that("validate_ferx_ir accepts IR with valid structural type", {
  ir <- new_ferx_ir(structural = list(type = "pk_macro", pk_call = "one_cpt_oral", pk_args = list()))
  expect_invisible(validate_ferx_ir(ir))
})

test_that("validate_ferx_ir rejects structural with missing type", {
  ir <- new_ferx_ir(structural = list(pk_call = "one_cpt_oral"))
  expect_error(validate_ferx_ir(ir), "structural\\$type is missing")
})

test_that("validate_ferx_ir rejects unknown structural type", {
  ir <- new_ferx_ir(structural = list(type = "compartmental"))
  expect_error(validate_ferx_ir(ir), "pk_macro")
})

test_that("validate_ferx_ir rejects odes without structural type ode", {
  ir <- new_ferx_ir(
    structural = list(type = "pk_macro", pk_call = "one_cpt_oral", pk_args = list()),
    odes       = list(list(state = "depot", rhs = "-KA * depot"))
  )
  expect_error(validate_ferx_ir(ir), "structural\\$type is not")
})

test_that("validate_ferx_ir rejects non-ferx_ir input", {
  expect_error(validate_ferx_ir(list()), "ferx_ir")
})

test_that("print.ferx_ir runs without error", {
  ir <- new_ferx_ir(
    source_format = "nonmem",
    source_file   = "run001.ctl",
    thetas        = list(list(name = "TVCL", init = 0.134, lower = 0.001, upper = 10)),
    omegas        = list(list(type = "diagonal", names = "ETA_CL", values = 0.07)),
    structural    = list(type = "pk_macro", pk_call = "one_cpt_oral", pk_args = list()),
    warnings      = c("WARN | check theta bounds"),
    unsupported   = c("MIXTURE model")
  )
  expect_no_error(print(ir))
})

test_that("print.ferx_ir handles empty IR without error", {
  expect_no_error(print(new_ferx_ir()))
})
