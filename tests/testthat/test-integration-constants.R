# rxode2 built-in constants (#76): the Tier 2 half.
#
# Full-pipeline translations, and the checks of .RX_CONSTANTS against rxode2's
# parser and R's Rmath.h, which read files. The inline unit tests for
# .normalise_expr() and the covariate classifiers are in test-rxui_to_ir.R; the
# weighted [covariates] case is in test-compartment-free.R.

test_that(".RX_CONSTANTS holds exactly the constants rxode2's parser knows", {
  # The other direction (review round 1, finding 5). A constant a later rxode2
  # adds is emitted verbatim -- a data column to ferx -- until it is in the
  # table, and a reserved name that is NOT a constant (`tlast`, `rxFlag`) must
  # never be in it: `.rxIsReservedName()` says yes to both, so it cannot tell.
  # The parser's own list is the reference: the `M_*` names
  # isReservedVariable() skips in parseVars.h, plus `pi`, which rxode2 handles
  # apart because a model may assign it.
  skip_if_not_installed("rxode2")
  hdr <- system.file("include", "parseVars.h", package = "rxode2")
  skip_if(!nzchar(hdr), "this rxode2 build ships no include/parseVars.h")
  src <- readLines(hdr, warn = FALSE)
  hits <- regmatches(src, regexpr('strcmp\\("M_[A-Za-z0-9_]+"', src))
  parser_m <- unique(sub('^strcmp\\("', "", sub('"$', "", hits)))
  # A header that exists but yields nothing means the pattern no longer fits
  # rxode2's source -- fail, so the guard cannot go quiet exactly when a new
  # constant could slip through.
  expect_gt(length(parser_m), 0L)
  expect_true(any(grepl('strcmp("pi", s)', src, fixed = TRUE)))
  table_m <- grep("^M_", names(.RX_CONSTANTS), value = TRUE)
  expect_setequal(table_m, parser_m)
  expect_setequal(setdiff(names(.RX_CONSTANTS), table_m), "pi")
})

test_that("every .RX_CONSTANTS value matches the Rmath.h macro of its name", {
  # rxode2 compiles these names to the C macros, so Rmath.h is the reference --
  # and an independent one: a wrong formula in the table (log2(10) for
  # M_LOG10_2, say) is caught here and nowhere else.
  hdr <- file.path(R.home("include"), "Rmath.h")
  skip_if_not(file.exists(hdr), "no Rmath.h in R.home('include')")
  defs <- regmatches(readLines(hdr),
                     regexec("^#define[[:space:]]+(M_[A-Za-z0-9_]+)[[:space:]]+([0-9.]+)",
                             readLines(hdr)))
  defs <- defs[lengths(defs) == 3L]
  rmath <- setNames(as.numeric(vapply(defs, `[`, "", 3L)),
                    vapply(defs, `[`, "", 2L))
  m <- names(.RX_CONSTANTS)[startsWith(names(.RX_CONSTANTS), "M_")]
  expect_true(all(m %in% names(rmath)))
  expect_equal(.RX_CONSTANTS[m], rmath[m], tolerance = 1e-15)
  expect_identical(.RX_CONSTANTS[["pi"]], pi)
})

test_that("pi in an nlmixr2 model is emitted as a literal, not a column", {
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ e0 <- 1; eta.e ~ 0.1; sd <- 0.5 })
    model({ E0 <- e0 + eta.e + sin(pi * TIME) + M_E
            pred <- E0; pred ~ add(sd) })
  }
  txt <- suppressWarnings(suppressMessages(
    nlmixr2_to_ferx(f, validate = FALSE)))$ferx_text
  expect_match(txt, "E0 = TVE0 + ETA_E + sin(3.14159265358979 * TIME) + 2.71828182845905",
               fixed = TRUE)
  expect_no_match(txt, "\\bpi\\b|M_E")
})

test_that("a pi the nlmixr2 model assigns keeps meaning the model's variable", {
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ e0 <- 1; eta.e ~ 0.1; sd <- 0.5 })
    model({ pi <- 2; E0 <- e0 + eta.e + sin(pi * TIME)
            pred <- E0; pred ~ add(sd) })
  }
  txt <- suppressWarnings(suppressMessages(
    nlmixr2_to_ferx(f, validate = FALSE)))$ferx_text
  expect_match(txt, "PI = 2", fixed = TRUE)
  expect_match(txt, "sin(PI * TIME)", fixed = TRUE)
  expect_no_match(txt, "3.14159")
})

test_that("pi in an ODE right-hand side no longer reads as an undeclared name", {
  # On main this raised `ERROR | [odes] references PI` and aborted under strict.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ tcl <- log(2); tv <- log(10); eta.cl ~ 0.1; sd <- 0.1 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(centr) <- -cl / v * centr * (1 + 0.1 * sin(2 * pi * time / 24))
            cp <- centr / v; cp ~ add(sd) })
  }
  ir <- suppressWarnings(suppressMessages(
    rxui_to_ir(rxode2::rxode2(f), source_format = "nlmixr2")))
  expect_no_match(ir$warnings, "\\[odes\\] references")
  expect_length(ir$unsupported, 0L)
  txt <- emit_ferx(ir)
  expect_match(txt, "sin(2 * 3.14159265358979 * time/24)", fixed = TRUE)
})

test_that("a STATE named pi survives the passes that re-parse emitted text", {
  # Review round 1, finding 1. A state keeps its lowercase name, so the partial
  # maps of the re-parse passes do not hold it, and substituting there turned
  # first-order elimination into zero-order. Two fixtures, one per pass: the
  # theta read from [odes] goes through .scope_odes_to_params(), the dose-shaped
  # `F1` through .rewrite_syms(). main got both right.
  skip_if_not_installed("rxode2")
  carrier <- function() {
    ini({ tcl <- log(2); v <- 10; eta.cl ~ 0.1; sd <- 0.1 })
    model({ k <- exp(tcl + eta.cl)
            d/dt(pi) <- -exp(tcl) / v * pi * k
            cp <- pi / v; cp ~ add(sd) })
  }
  dose_attr <- function() {
    ini({ cl <- 1; v <- 10; eta.cl ~ 0.1; sd <- 0.1 })
    model({ cl <- cl * exp(eta.cl); F1 <- 2
            d/dt(pi) <- -cl / v * pi * F1
            cp <- pi / v; cp ~ add(sd) })
  }
  odes <- function(f) {
    txt <- suppressWarnings(suppressMessages(
      nlmixr2_to_ferx(f, validate = FALSE)))$ferx_text
    grep("d/dt(pi)", strsplit(txt, "\n")[[1]], fixed = TRUE, value = TRUE)
  }
  expect_identical(trimws(odes(carrier)), "d/dt(pi) = -exp(TCL)/V * pi * K")
  expect_identical(trimws(odes(dose_attr)), "d/dt(pi) = -CL/V * pi * F1_PAR")
})

test_that("pi in an initial condition is emitted as a literal", {
  # Review round 1, finding 2. main dropped this statement with an ERROR, so the
  # compartment started at 0 instead of 100*pi.
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({ cl <- 1; v <- 10; eta.cl ~ 0.1; sd <- 0.1 })
    model({ cl <- cl * exp(eta.cl); centr(0) <- 100 * pi
            d/dt(centr) <- -cl / v * centr
            cp <- centr / v; cp ~ add(sd) })
  }
  r <- suppressWarnings(suppressMessages(nlmixr2_to_ferx(f, validate = FALSE)))
  expect_match(r$ferx_text, "init(centr) = 100 * 3.14159265358979", fixed = TRUE)
  expect_no_match(r$ferx_text, "\\bpi\\b|\\bPI\\b")
  expect_length(r$unsupported, 0L)
})
