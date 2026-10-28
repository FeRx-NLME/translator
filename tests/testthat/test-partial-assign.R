# A name assigned on only some paths through an `if` and read afterwards is
# refused by ferx-core (E_PARSE). .default_partial_assignments() gives it the
# default NONMEM would have left it with.

A <- function(lhs, rhs) list(lhs = lhs, rhs = rhs)
IF <- function(cond, then, else_ = list())
  list(kind = "if", cond = cond, then = then, else_ = else_)

.lhs <- function(stmts) unlist(lapply(stmts, function(s) s$lhs))

# -- Tier 1: .assigned_every_path ------------------------------------------------

test_that(".assigned_every_path needs an else, and every arm to assign", {
  expect_equal(.assigned_every_path(list(A("x", "1"))), "X")
  expect_equal(.assigned_every_path(list(IF("a > 0", list(A("x", "1"))))),
               character())
  expect_equal(
    .assigned_every_path(list(IF("a > 0", list(A("x", "1"), A("y", "2")),
                                 list(A("x", "3"))))),
    "X")
  # an else-less nested if does not count as assigning on every path
  expect_equal(
    .assigned_every_path(list(IF("a > 0",
                                 list(IF("b > 0", list(A("x", "1")))),
                                 list(A("x", "3"))))),
    character())
})

# -- Tier 1: .default_partial_assignments ----------------------------------------

test_that("an else-less if whose name is read later gets a default before it", {
  st  <- list(IF("sex == 0", list(A("tvcl", "t1"))), A("cl", "tvcl * 2"))
  out <- .default_partial_assignments(st)
  expect_equal(.lhs(out$stmts[1]), "tvcl")
  expect_equal(out$stmts[[1]]$rhs, "0")
  expect_identical(out$stmts[[2]], st[[1]])
  expect_length(out$stmts, 3L)
  expect_match(out$warnings, "^WARN ", all = FALSE)
  expect_length(out$warnings, 1L)
})

test_that("two else-less ifs on one name get ONE default, before the first", {
  st <- list(IF("s == 0", list(A("tvcl", "a"))),
             IF("s == 1", list(A("tvcl", "b"))),
             A("cl", "tvcl"))
  out <- .default_partial_assignments(st)
  expect_equal(vapply(out$stmts, function(s) s$kind %||% "assign", ""),
               c("assign", "if", "if", "assign"))
  expect_equal(out$stmts[[1]]$lhs, "tvcl")
  expect_length(out$warnings, 1L)
})

test_that("an if/else assigning on every path needs no default", {
  st <- list(IF("s == 0", list(A("tvcl", "a")), list(A("tvcl", "b"))),
             A("cl", "tvcl"))
  out <- .default_partial_assignments(st)
  expect_identical(out$stmts, st)
  expect_length(out$warnings, 0L)
})

test_that("a name already assigned before the if is left alone", {
  st <- list(A("tvcl", "a"), IF("s == 1", list(A("tvcl", "b"))), A("cl", "tvcl"))
  out <- .default_partial_assignments(st)
  expect_identical(out$stmts, st)
  expect_length(out$warnings, 0L)
})

test_that("a name nothing reads is left alone, one read elsewhere is not", {
  st <- list(IF("s == 1", list(A("tmp", "1"))), A("cl", "2"))
  expect_identical(.default_partial_assignments(st)$stmts, st)
  # read by another block (an ODE right-hand side, say)
  out <- .default_partial_assignments(st, external_reads = "TMP")
  expect_equal(out$stmts[[1]]$lhs, "tmp")
})

test_that("a name read only inside a later if condition counts as read", {
  st <- list(IF("s == 1", list(A("flag", "1"))),
             IF("flag > 0", list(A("cl", "3")), list(A("cl", "4"))))
  out <- .default_partial_assignments(st)
  expect_equal(out$stmts[[1]]$lhs, "flag")
})

test_that("the default for a name in a nested else-less if is placed before the outer if", {
  st <- list(IF("a > 0", list(IF("b > 0", list(A("x", "1"))), A("x2", "2")),
                list(A("x", "3"))),
             A("cl", "x + x2"))
  out <- .default_partial_assignments(st)
  expect_setequal(.lhs(out$stmts[1:2]), c("x", "x2"))
  expect_equal(out$stmts[[3]]$kind, "if")
})

test_that("dose attributes default to the engine defaults: F -> 1, lag -> 0", {
  st <- list(IF("s == 1", list(A("F1", "0.8"), A("ALAG1", "0.5"),
                               A("D1", "2"))))
  out <- .default_partial_assignments(st)
  vals <- setNames(vapply(out$stmts[1:3], function(s) s$rhs, ""),
                   vapply(out$stmts[1:3], function(s) s$lhs, ""))
  expect_equal(vals[["F1"]], "1")
  expect_equal(vals[["ALAG1"]], "0")
  expect_equal(vals[["D1"]], "0")
  # dose attributes are read by the engine even when no statement reads them,
  # and an INFO, not a WARN, says so
  expect_true(all(grepl("^INFO ", out$warnings)))
  expect_length(out$warnings, 3L)
})

test_that("F is never defaulted to 0", {
  st <- list(IF("s == 1", list(A("F", "0.5"), A("F2", "0.5"))))
  out <- .default_partial_assignments(st)
  expect_false(any(vapply(out$stmts[1:2], function(s) s$rhs == "0", NA)))
})

# -- Tier 2: the bundled model ---------------------------------------------------

test_that("the two-IF idiom translates with a default, once, and says so", {
  res <- suppressWarnings(nm_to_ferx(nm_path("if_no_else.ctl")))
  txt <- strsplit(res$ferx_text, "\n")[[1]]
  trim <- trimws(txt)
  i_def <- which(trim == "TVCL = 0")
  i_if  <- grep("^if \\(SEX == 0\\)", trim)
  expect_length(i_def, 1L)
  expect_length(i_if, 1L)
  expect_lt(i_def, i_if)
  expect_equal(sum(grepl("^WARN.*TVCL is assigned in only some branches",
                         res$warnings)), 1L)
  expect_length(res$unsupported, 0L)
})

test_that("the translated two-IF model is accepted by the engine", {
  skip_if_not_installed("ferx")
  res <- suppressWarnings(nm_to_ferx(nm_path("if_no_else.ctl")))
  tmp <- tempfile(fileext = ".ferx")
  writeLines(res$ferx_text, tmp)
  on.exit(unlink(tmp))
  out <- NULL
  utils::capture.output(out <- ferx::ferx_model_validate(tmp))
  expect_true(out$ok, info = paste(out$diagnostics$message, collapse = "\n"))
})

test_that("without the default the engine (once it checks) rejects the same model", {
  skip_if_not_installed("ferx")
  res <- suppressWarnings(nm_to_ferx(nm_path("if_no_else.ctl")))
  bare <- sub("(?m)^\\s*TVCL = 0\\n", "", res$ferx_text, perl = TRUE)
  expect_false(identical(bare, res$ferx_text))
  tmp <- tempfile(fileext = ".ferx")
  writeLines(bare, tmp)
  on.exit(unlink(tmp))
  out <- NULL
  utils::capture.output(out <- ferx::ferx_model_validate(tmp))
  # The engine pin used in CI predates this check; the assertion only bites on
  # an engine that has it, and that is when it proves the default is needed.
  if (isTRUE(out$ok))
    skip("installed ferx predates the every-path assignment check")
  expect_match(paste(out$diagnostics$message, collapse = " "),
               "every path", fixed = TRUE)
})
