# Issue #26: a $MODEL compartment declared after every DADT-bearing one makes
# nonmem2rx drop it and name the rest `rxddta<N>` (N = DADT/A index). All
# fixtures here are self-written.

.trail_ctl <- function(comps, des, y = "A(2)", s = "S2 = V",
                       dir = tmp_ctl_dir(), name = "trail.ctl") {
  ctl <- file.path(dir, name)
  writeLines(c(
    "$PROBLEM trailing $MODEL compartment with no DADT",
    "$INPUT ID TIME DV AMT EVID MDV CMT",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN13 TOL=6",
    "$MODEL", paste0("  ", comps),
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    paste0("  ", s),
    "$DES", paste0("  ", des),
    "$ERROR",
    paste0("  Y = ", y, "*(1 + EPS(1))"),
    "$THETA (0,5) (0,50) (0,1)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)
  ctl
}

.trail_translate <- function(...)
  suppressWarnings(suppressMessages(nm_to_ferx(.trail_ctl(...), validate = FALSE)))

.states_line <- function(r)
  grep("ode(", strsplit(r$ferx_text, "\n")[[1]], fixed = TRUE, value = TRUE)

test_that("a trailing $MODEL compartment keeps every state name", {
  skip_if_not_installed("nonmem2rx")
  r <- .trail_translate(
    c("COMP=(DEPOT, DEFDOSE)", "COMP=(CENTRAL, DEFOBS)", "COMP=(DUMMY)"),
    c("DADT(1) = -KA*A(1)", "DADT(2) =  KA*A(1) - (CL/V)*A(2)"))
  expect_false(grepl("rxddta", r$ferx_text, ignore.case = TRUE))
  expect_match(.states_line(r), "obs_cmt=CENTRAL", fixed = TRUE)
  expect_match(.states_line(r), "states=[DEPOT, CENTRAL, DUMMY]", fixed = TRUE)
  expect_match(r$ferx_text, "d/dt(DEPOT) = -KA * DEPOT", fixed = TRUE)
  expect_match(r$ferx_text, "d/dt(CENTRAL) = KA * DEPOT - (CL/V) * CENTRAL",
               fixed = TRUE)
  # The trailing compartment has no equation; it is kept so a CMT naming it is
  # still a compartment.
  expect_match(r$ferx_text, "d/dt(DUMMY) = 0", fixed = TRUE)
})

test_that("the false DEFOBS contradiction and the numbering ERROR are gone", {
  skip_if_not_installed("nonmem2rx")
  r <- .trail_translate(
    c("COMP=(DEPOT, DEFDOSE)", "COMP=(CENTRAL, DEFOBS)", "COMP=(DUMMY)"),
    c("DADT(1) = -KA*A(1)", "DADT(2) =  KA*A(1) - (CL/V)*A(2)"))
  expect_length(grep("as DEFOBS", r$warnings), 0L)
  expect_length(grep("^ERROR", r$warnings), 0L)
  expect_length(r$unsupported, 0L)
  note <- grep("lost the compartment names", r$warnings, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "^INFO")
  expect_match(note, "'DUMMY'", fixed = TRUE)
})

test_that("the numbering is $MODEL order when $DES is written in another order", {
  skip_if_not_installed("nonmem2rx")
  # DADT(2) precedes DADT(1): d/dt statement order is [CENTRAL, DEPOT], COMP
  # order is [DEPOT, CENTRAL]. Reading the placeholders by statement position
  # would swap the two names.
  r <- .trail_translate(
    c("COMP=(DEPOT, DEFDOSE)", "COMP=(CENTRAL, DEFOBS)", "COMP=(DUMMY)"),
    c("DADT(2) =  KA*A(1) - (CL/V)*A(2)", "DADT(1) = -KA*A(1)"))
  expect_match(.states_line(r), "states=[DEPOT, CENTRAL, DUMMY]", fixed = TRUE)
  expect_match(.states_line(r), "obs_cmt=CENTRAL", fixed = TRUE)
  expect_match(r$ferx_text, "d/dt(DEPOT) = -KA * DEPOT", fixed = TRUE)
  expect_match(r$ferx_text, "d/dt(CENTRAL) = KA * DEPOT - (CL/V) * CENTRAL",
               fixed = TRUE)
})

test_that("a COMP order that differs from the model's natural order is respected", {
  skip_if_not_installed("nonmem2rx")
  # CENTRAL is compartment 1 here, so DADT(1) is the central equation. The
  # equations are the same model as above with the numbers swapped.
  r <- .trail_translate(
    c("COMP=(CENTRAL, DEFOBS)", "COMP=(DEPOT, DEFDOSE)", "COMP=(DUMMY)"),
    c("DADT(2) = -KA*A(2)", "DADT(1) =  KA*A(2) - (CL/V)*A(1)"),
    y = "A(1)", s = "S1 = V")
  expect_match(.states_line(r), "states=[CENTRAL, DEPOT, DUMMY]", fixed = TRUE)
  expect_match(.states_line(r), "obs_cmt=CENTRAL", fixed = TRUE)
  expect_match(r$ferx_text, "d/dt(DEPOT) = -KA * DEPOT", fixed = TRUE)
  expect_match(r$ferx_text, "d/dt(CENTRAL) = KA * DEPOT - (CL/V) * CENTRAL",
               fixed = TRUE)
})

test_that("several trailing compartments are all restored", {
  skip_if_not_installed("nonmem2rx")
  r <- .trail_translate(
    c("COMP=(DEPOT, DEFDOSE)", "COMP=(CENTRAL, DEFOBS)", "COMP=(SPARE1)",
      "COMP=(SPARE2)"),
    c("DADT(1) = -KA*A(1)", "DADT(2) =  KA*A(1) - (CL/V)*A(2)"))
  expect_match(.states_line(r), "states=[DEPOT, CENTRAL, SPARE1, SPARE2]",
               fixed = TRUE)
  expect_match(r$ferx_text, "d/dt(SPARE1) = 0", fixed = TRUE)
  expect_match(r$ferx_text, "d/dt(SPARE2) = 0", fixed = TRUE)
  expect_match(grep("lost the compartment names", r$warnings, value = TRUE),
               "'SPARE1', 'SPARE2'", fixed = TRUE)
})

test_that("the restored model validates against the engine", {
  skip_if_not_installed("nonmem2rx")
  skip_if_not_installed("ferx")
  skip_if_not(nzchar(Sys.getenv("FERXTRANSLATE_ENGINE_TESTS")) ||
              !nzchar(Sys.getenv("CI")))
  ctl <- .trail_ctl(
    c("COMP=(DEPOT, DEFDOSE)", "COMP=(CENTRAL, DEFOBS)", "COMP=(DUMMY)"),
    c("DADT(2) =  KA*A(1) - (CL/V)*A(2)", "DADT(1) = -KA*A(1)"))
  r <- suppressWarnings(suppressMessages(nm_to_ferx(ctl, validate = TRUE)))
  expect_true(r$validation$ok)
})

test_that(".restore_placeholder_states declines unless the mapping is total", {
  # `d/dt(x)` parses as `/`(d, dt(x)); build it explicitly so the fixture is the
  # shape nonmem2rx produces.
  mk <- function(n, rhs) call("<-", call("/", as.name("d"),
                                         call("dt", as.name(paste0("rxddta", n)))),
                              rhs)
  lst <- list(mk(1, quote(-ka * rxddta1)),
              mk(2, quote(ka * rxddta1 - k * rxddta2)))
  ok <- .restore_placeholder_states(lst, c("DEPOT", "CENTRAL"))
  expect_identical(.state_raw_names(ok$lst), c("DEPOT", "CENTRAL"))
  expect_length(ok$warnings, 1L)

  same <- function(r) { expect_identical(r$lst, lst); expect_length(r$warnings, 0L) }
  # No COMP list (a non-NONMEM source).
  same(.restore_placeholder_states(lst, NULL))
  # A placeholder index past the COMP list: the mapping is not total.
  same(.restore_placeholder_states(lst, "DEPOT"))
  # Two COMP names that collapse to one identifier.
  same(.restore_placeholder_states(lst, c("A-B", "A.B")))
  same(.restore_placeholder_states(lst, c("DEPOT", "depot")))
  # A state that is not a placeholder: names were not lost.
  named <- list(call("<-", call("/", as.name("d"), call("dt", as.name("depot"))), 0))
  expect_identical(.restore_placeholder_states(named, c("DEPOT", "CENTRAL"))$lst,
                   named)
})
