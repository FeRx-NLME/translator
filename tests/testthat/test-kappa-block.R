# Correlated IOV -> block_kappa (issue #40)

.kap_ir <- function(...) {
  new_ferx_ir(
    kappas = list(
      list(name = "IOV_CL", value = 0.05, fixed = FALSE, block = 1L, block_idx = 1L),
      list(name = "IOV_V",  value = 0.03, fixed = FALSE, block = 1L, block_idx = 2L)),
    kappa_blocks = list(list(id = 1L,
                             cov = matrix(c(0.05, 0.01, 0.01, 0.03), 2, 2),
                             fixed = FALSE)),
    ...)
}

test_that("a kappa block emits one block_kappa on the variance scale", {
  txt <- emit_ferx(.kap_ir())
  expect_match(txt, "  block_kappa (IOV_CL, IOV_V) = [0.05, 0.01, 0.03]", fixed = TRUE)
  expect_no_match(txt, "  kappa IOV", fixed = TRUE)
})

test_that("FIX goes on the whole block_kappa", {
  ir <- .kap_ir()
  ir$kappa_blocks[[1]]$fixed <- TRUE
  expect_match(emit_ferx(ir), "= [0.05, 0.01, 0.03] FIX", fixed = TRUE)
})

test_that("the kappa triangle follows the members' order", {
  ir <- .kap_ir()
  ir$kappas <- rev(ir$kappas)
  expect_match(emit_ferx(ir),
               "block_kappa (IOV_V, IOV_CL) = [0.03, 0.01, 0.05]", fixed = TRUE)
})

test_that("a plain kappa beside a block keeps its own line", {
  ir <- .kap_ir()
  ir$kappas <- c(ir$kappas, list(list(name = "IOV_KA", value = 0.2)))
  txt <- emit_ferx(ir)
  expect_match(txt, "block_kappa (IOV_CL, IOV_V)", fixed = TRUE)
  expect_match(txt, "  kappa IOV_KA ~ 0.2", fixed = TRUE)
})

test_that("validate_ferx_ir rejects a kappa block that does not match its members", {
  ir <- .kap_ir()
  ir$kappa_blocks[[1]]$cov <- matrix(0.05, 1, 1)
  expect_error(validate_ferx_ir(ir), "kappa_blocks")
})

test_that("an nlmixr2 IOV block translates to block_kappa with no gap reported", {
  skip_if_not_installed("rxode2")
  source(r2_path("iov_block_nlmixr2.R"), local = TRUE)
  res <- suppressWarnings(nlmixr2_to_ferx(f_iov_block))
  expect_match(res$ferx_text,
               "block_kappa (IOV_CL, IOV_V) = [0.05, 0.01, 0.03]", fixed = TRUE)
  expect_length(res$unsupported, 0L)
  expect_false(any(grepl("IOV block omega", res$warnings)))
})

test_that("the engine accepts block_kappa and it changes the fit", {
  skip_if_not_installed("ferx")
  skip_if_not_installed("rxode2")
  skip_on_cran()
  skip_if(Sys.getenv("FERXTRANSLATE_ENGINE_TESTS") != "true" &&
            !identical(Sys.getenv("NOT_CRAN"), "true"))
  source(r2_path("iov_block_nlmixr2.R"), local = TRUE)
  res <- suppressWarnings(nlmixr2_to_ferx(f_iov_block))

  # Synthetic data: two occasions per subject, a dose at the start of each.
  set.seed(1)
  rows <- do.call(rbind, lapply(1:12, function(id) {
    do.call(rbind, lapply(1:2, function(occ) {
      t0 <- (occ - 1) * 24
      cl <- 0.134 * exp(rnorm(1, 0, 0.2)); v <- 8.1 * exp(rnorm(1, 0, 0.2))
      tt <- t0 + c(0.5, 2, 6, 12)
      conc <- 100 / v * exp(-cl / v * (tt - t0)) * exp(rnorm(4, 0, 0.1))
      rbind(data.frame(ID = id, TIME = t0, DV = 0, EVID = 1L, AMT = 100, CMT = 1L, MDV = 1L, OCC = occ),
            data.frame(ID = id, TIME = tt, DV = conc, EVID = 0L, AMT = 0, CMT = 1L, MDV = 0L, OCC = occ))
    }))
  }))
  d <- tempfile("kapblk-"); dir.create(d)
  data_file <- file.path(d, "d.csv")
  write.csv(rows, data_file, row.names = FALSE, quote = FALSE)

  ofv_of <- function(txt) {
    txt <- sub("maxiter = 500", "maxiter = 0", txt, fixed = TRUE)
    txt <- sub("covariance = true", "covariance = false", txt, fixed = TRUE)
    txt <- sub("iov_column = occ", "iov_column = OCC", txt, fixed = TRUE)
    f <- tempfile(tmpdir = d, fileext = ".ferx")
    writeLines(txt, f)
    invisible(utils::capture.output(fit <- suppressMessages(
      ferx::ferx_fit(f, data = data_file))))
    fit$ofv
  }
  block <- ofv_of(res$ferx_text)
  flat  <- ofv_of(sub("block_kappa (IOV_CL, IOV_V) = [0.05, 0.01, 0.03]",
                      "kappa IOV_CL ~ 0.05\n  kappa IOV_V ~ 0.03",
                      res$ferx_text, fixed = TRUE))
  expect_true(is.finite(block))
  # Discriminating: dropping the covariance is what the translator used to do.
  expect_false(isTRUE(all.equal(block, flat, tolerance = 1e-6)))
})
