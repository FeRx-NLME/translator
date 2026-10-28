# $SIGMA BLOCK -> block_sigma (issue #39)

.sig_ir <- function(...) {
  new_ferx_ir(
    sigmas = list(
      list(name = "EPS1", value = 0.2, scale = "sd", fixed = FALSE, block = 1L, block_idx = 1L),
      list(name = "EPS2", value = 0.1, scale = "sd", fixed = FALSE, block = 1L, block_idx = 2L)),
    sigma_blocks = list(list(id = 1L,
                             cov = matrix(c(0.04, 0.001, 0.001, 0.01), 2, 2),
                             fixed = FALSE)),
    ...)
}

test_that("a sigma block emits one block_sigma on the variance scale", {
  txt <- emit_ferx(.sig_ir())
  expect_match(txt, "  block_sigma (EPS1, EPS2) = [0.04, 0.001, 0.01]", fixed = TRUE)
  expect_no_match(txt, "  sigma EPS", fixed = TRUE)
})

test_that("FIX goes on the whole block_sigma", {
  ir <- .sig_ir()
  ir$sigma_blocks[[1]]$fixed <- TRUE
  expect_match(emit_ferx(ir), "= [0.04, 0.001, 0.01] FIX", fixed = TRUE)
})

test_that("the lower triangle follows the members' declaration order", {
  # ferx reads the triangle positionally, so reordering the names without
  # permuting the triangle would swap the two variances.
  ir <- .sig_ir()
  ir$sigmas <- rev(ir$sigmas)
  expect_match(emit_ferx(ir),
               "block_sigma (EPS2, EPS1) = [0.01, 0.001, 0.04]", fixed = TRUE)
})

test_that("a plain sigma beside a block keeps its own line", {
  ir <- .sig_ir()
  ir$sigmas <- c(ir$sigmas, list(list(name = "EPS3", value = 0.3, scale = "sd")))
  txt <- emit_ferx(ir)
  expect_match(txt, "block_sigma (EPS1, EPS2)", fixed = TRUE)
  expect_match(txt, "  sigma EPS3 ~ 0.3 (sd)", fixed = TRUE)
})

test_that("validate_ferx_ir rejects a block that does not match its members", {
  ir <- .sig_ir()
  ir$sigma_blocks[[1]]$cov <- matrix(0.04, 1, 1)
  expect_error(validate_ferx_ir(ir), "sigma_blocks")
})

test_that("a $SIGMA BLOCK translates to block_sigma with no gap reported", {
  res <- suppressWarnings(nm_to_ferx(nm_path("sigma_block.ctl")))
  expect_match(res$ferx_text, "block_sigma (EPS1, EPS2) = [0.04, 0.001, 0.01]",
               fixed = TRUE)
  expect_length(res$unsupported, 0L)
  expect_false(any(grepl("residual covariances", res$warnings)))
})

test_that("the engine accepts the emitted block_sigma and it changes the fit", {
  skip_if_not_installed("ferx")
  skip_on_cran()
  skip_if(Sys.getenv("FERXTRANSLATE_ENGINE_TESTS") != "true" &&
            !identical(Sys.getenv("NOT_CRAN"), "true"))
  res  <- suppressWarnings(nm_to_ferx(nm_path("sigma_block.ctl")))
  data <- system.file("testdata", "ode_1cpt_oral_concordance.csv",
                      package = "ferxtranslate")
  d <- tempfile("sigblk-"); dir.create(d)
  ofv_of <- function(txt) {
    txt <- sub("maxiter = 500", "maxiter = 0", txt, fixed = TRUE)
    txt <- sub("covariance = true", "covariance = false", txt, fixed = TRUE)
    f <- tempfile(tmpdir = d, fileext = ".ferx")
    writeLines(txt, f)
    invisible(utils::capture.output(fit <- suppressMessages(
      ferx::ferx_fit(f, data = data))))
    fit$ofv
  }
  block <- ofv_of(res$ferx_text)
  flat  <- ofv_of(sub("block_sigma (EPS1, EPS2) = [0.04, 0.001, 0.01]",
                      "sigma EPS1 ~ 0.2 (sd)\n  sigma EPS2 ~ 0.1 (sd)",
                      res$ferx_text, fixed = TRUE))
  expect_true(is.finite(block))
  # Discriminating: dropping the covariance is what the translator used to do,
  # and must give a different objective on this fixture.
  expect_false(isTRUE(all.equal(block, flat, tolerance = 1e-6)))
})
