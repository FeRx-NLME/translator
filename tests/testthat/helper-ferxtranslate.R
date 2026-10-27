# Shared test helpers. testthat sources helper-*.R before every test file, so
# these live here rather than being duplicated per file (nm_path / nm_path_t).

nm_path <- function(file) {
  system.file("testmodels", "nonmem", file, package = "ferxtranslate",
              mustWork = TRUE)
}

r2_path <- function(file) {
  system.file("testmodels", "nlmixr2", file, package = "ferxtranslate",
              mustWork = TRUE)
}

# Every bundled NONMEM model. Filtered by extension so a README, dataset or
# .lst dropped into the directory does not get fed to nonmem2rx().
.bundled_nm_models <- function() {
  list.files(system.file("testmodels", "nonmem", package = "ferxtranslate"),
             pattern = "\\.(ctl|mod)$", full.names = TRUE)
}

# Base-R scratch directory; withr is not a declared dependency of this package,
# and adding one for a two-line helper would put a WARNING in
# `R CMD check --as-cran` ("'::' or ':::' import not declared"). R clears
# tempdir() at session end, so no explicit cleanup is needed.
#
# Lives here rather than in one test file because more than one file needs it:
# defined inside test-translate.R it was invisible to test-rxui_to_ir.R, which
# is the kind of breakage that shows up as an opaque "could not find function".
tmp_ctl_dir <- function() {
  d <- tempfile("ferxtr-")
  dir.create(d)
  d
}

# Minimal iniDf for the .extract_* unit tests. Only the columns those extractors
# read; building the real 12-column frame would assert nothing extra and would
# drift with upstream.
#
# Here rather than in one test file because issues #39, #40 and #41 each need a
# synthetic iniDf too, and a helper defined mid-file gets copy-pasted instead of
# reused -- which is what the note at the top of this file exists to prevent.
.mk_ini <- function(...) {
  rows <- list(...)
  do.call(rbind, lapply(rows, function(r) data.frame(
    ntheta = r$ntheta %||% NA_integer_, neta1 = r$neta1 %||% NA_integer_,
    neta2  = r$neta2  %||% NA_integer_, name = r$name, est = r$est,
    fix = r$fix, err = r$err %||% NA_character_,
    condition = r$condition %||% NA_character_,
    stringsAsFactors = FALSE)))
}

# One-compartment IV bolus control stream whose $ERROR is whatever the test
# supplies. `solved = TRUE` is ADVAN1; otherwise the same model as an ODE.
.readout_ctl <- function(err, solved = FALSE, pk_extra = character(),
                         theta_extra = character(),
                         omega_extra = character(),
                         y_line = "  Y = IPRED*(1+EPS(1))") {
  ctl <- file.path(tmp_ctl_dir(), "m.ctl")
  writeLines(c(
    "$PROBLEM readout arithmetic",
    "$INPUT ID TIME DV AMT EVID MDV CMT",
    "$DATA d.csv IGNORE=@",
    if (solved) "$SUBROUTINES ADVAN1 TRANS2" else c("$SUBROUTINES ADVAN6 TOL=6",
                                                    "$MODEL COMP=(CENT)"),
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    pk_extra,
    if (!solved) c("$DES", "  DADT(1) = -CL/V*A(1)"),
    "$ERROR", err,
    y_line,
    "$THETA", "  (0,3)", "  (0,30)", theta_extra,
    "$OMEGA 0.09", omega_extra, "$SIGMA 0.01", "$EST METHOD=1 INTER"), ctl)
  ctl
}

.scaling_lines <- function(txt) {
  ln <- strsplit(txt, "\n")[[1]]
  i  <- match("[scaling]", ln)
  if (is.na(i)) return(character())
  ln[i:(i + match(TRUE, !nzchar(ln[-seq_len(i)]), nomatch = length(ln)) - 1L)][-1L]
}
