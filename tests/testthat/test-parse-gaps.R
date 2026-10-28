# Control streams nonmem2rx cannot parse come back as a gap, not a raw parse
# error. The fixtures are written inline: they cannot be bundled, because the
# corpus sweeps translate every file in inst/testmodels/nonmem/.

.gap_ctl <- function(name, body) {
  d <- tmp_ctl_dir()
  p <- file.path(d, name)
  writeLines(body, p)
  p
}

.gap_streams <- list(
  advan9 = c("$PROBLEM ADVAN9 one-compartment",
    "$INPUT ID TIME AMT DV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN9 TOL=6",
    "$MODEL COMP=(CENT) COMP=(DEPOT)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    "  S1 = V",
    "$DES",
    "  DADT(1) = KA*A(2) - CL/V*A(1)",
    "  DADT(2) = -KA*A(2)",
    "$ERROR",
    "  IPRED = A(1)/V",
    "  Y = IPRED*(1+EPS(1))",
    "$THETA (0,3) (0,30) (0,1)",
    "$OMEGA 0.04",
    "$SIGMA 0.02",
    "$ESTIMATION METHOD=1 MAXEVAL=99"),
  advan10 = c("$PROBLEM ADVAN10 michaelis-menten",
    "$INPUT ID TIME AMT DV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN10 TRANS1",
    "$PK",
    "  VM = THETA(1)*EXP(ETA(1))",
    "  KM = THETA(2)",
    "  V  = THETA(3)",
    "  S1 = V",
    "$ERROR",
    "  IPRED = F",
    "  Y = IPRED*(1+EPS(1))",
    "$THETA (0,3) (0,30) (0,10)",
    "$OMEGA 0.04",
    "$SIGMA 0.02",
    "$ESTIMATION METHOD=1 MAXEVAL=99"),
  solved_periph = c("$PROBLEM solved system reading peripheral",
    "$INPUT ID TIME AMT DV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN4 TRANS4",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V2 = THETA(2)",
    "  Q  = THETA(3)",
    "  V3 = THETA(4)",
    "  KA = THETA(5)",
    "  S2 = V2",
    "$ERROR",
    "  IPRED = F",
    "  PER = A(3)/V3",
    "  Y = IPRED*(1+EPS(1))",
    "$THETA (0,3) (0,30) (0,2) (0,40) (0,1)",
    "$OMEGA 0.04",
    "$SIGMA 0.02",
    "$ESTIMATION METHOD=1 MAXEVAL=99")
)

# -- Tier 1: classification -----------------------------------------------------

test_that("each known upstream message maps to its NONMEM feature", {
  expect_equal(.classify_nm_parse_error(
    "Differential Algebra Equations are not supported in translation (ADVAN9 or ADVAN15)")$id,
    "dae")
  expect_equal(.classify_nm_parse_error(
    "Michelis Menton model translation not supported (ADVAN10)")$id, "mm")
  # the corrected spelling, should upstream fix it
  expect_equal(.classify_nm_parse_error(
    "Michaelis Menten model translation not supported (ADVAN10)")$id, "mm")
  expect_equal(.classify_nm_parse_error(paste0(
    "can only request depot and central compartments for solved systems ",
    "in rxode2 translations"))$id, "solved_cmt")
})

test_that("an unrecognised message is not classified", {
  expect_null(.classify_nm_parse_error("rxode2 syntax error: :007: df <- cl/v"))
  expect_null(.classify_nm_parse_error("something nobody has seen"))
  expect_null(.classify_nm_parse_error(character()))
})

# -- Tier 2: through nm_to_ferx() ----------------------------------------------

test_that("each unparseable class returns a gap result under strict = FALSE", {
  skip_if_not_installed("nonmem2rx")
  expected <- c(advan9 = "ADVAN9", advan10 = "ADVAN10", solved_periph = "ADVAN4")
  for (nm in names(.gap_streams)) {
    p <- .gap_ctl(paste0(nm, ".ctl"), .gap_streams[[nm]])
    res <- suppressMessages(suppressWarnings(
      nm_to_ferx(p, validate = FALSE, strict = FALSE)))
    expect_s3_class(res, "ferx_translate_result")
    expect_length(res$unsupported, 1L)
    expect_match(res$unsupported, "unparseable by nonmem2rx", info = nm)
    expect_match(res$unsupported, expected[[nm]], info = nm)
    expect_true(any(startsWith(res$warnings, "ERROR")), info = nm)
    expect_match(res$warnings, "nonmem2rx", all = FALSE, info = nm)
    expect_false(grepl("^\\[", res$ferx_text))
  }
})

test_that("the default strict = TRUE aborts, with the result on the condition", {
  skip_if_not_installed("nonmem2rx")
  p <- .gap_ctl("advan10.ctl", .gap_streams$advan10)
  err <- tryCatch(suppressMessages(nm_to_ferx(p, validate = FALSE)),
                  ferxtranslate_unparseable = function(e) e)
  expect_s3_class(err, "ferxtranslate_unparseable")
  expect_s3_class(err$result, "ferx_translate_result")
  expect_match(err$result$unsupported, "ADVAN10")
})

test_that("a recognised gap never writes the output file", {
  skip_if_not_installed("nonmem2rx")
  p   <- .gap_ctl("advan10.ctl", .gap_streams$advan10)
  out <- file.path(dirname(p), "out.ferx")
  suppressMessages(suppressWarnings(
    nm_to_ferx(p, output = out, validate = FALSE, strict = FALSE)))
  expect_false(file.exists(out))
})

test_that("an unrecognised parse failure still raises", {
  skip_if_not_installed("nonmem2rx")
  # A stream nonmem2rx rejects for a reason this package does not classify:
  # without this half, classification would absorb the next upstream change.
  p <- .gap_ctl("garbage.ctl", c("$PROBLEM nothing", "$INPUT ID TIME DV",
                                 "$DATA d.csv", "$SUBROUTINE ADVAN99"))
  expect_error(suppressMessages(suppressWarnings(
    nm_to_ferx(p, validate = FALSE, strict = FALSE))),
    "Failed to parse")
})

test_that("a translatable model is unaffected", {
  skip_if_not_installed("nonmem2rx")
  res <- suppressWarnings(nm_to_ferx(nm_path("1cpt_oral.ctl"), validate = FALSE))
  expect_length(grep("unparseable", res$unsupported), 0L)
  expect_match(res$ferx_text, "[individual_parameters]", fixed = TRUE)
})

# -- Tier 1: names rxode2 rejects ----------------------------------------------

.reserved_stream <- c(
  "$PROBLEM reserved name in DES",
  "$INPUT ID TIME AMT DV",
  "$DATA d.csv IGNORE=@",
  "$SUBROUTINE ADVAN13 TOL=6",
  "$MODEL COMP=(DEPOT) COMP=(CENT)",
  "$PK",
  "  CL = THETA(1)*EXP(ETA(1))",
  "  V  = THETA(2)",
  "  KA = THETA(3)",
  "  S2 = V",
  "$DES",
  "  DF = CL/V",
  "  DADT(1) = -KA*A(1)",
  "  DADT(2) = KA*A(1) - DF*A(2)",
  "$ERROR",
  "  IPRED = F",
  "  Y = IPRED*(1+EPS(1))",
  "$THETA (0,3) (0,30) (0,1)",
  "$OMEGA 0.04",
  "$SIGMA 0.02",
  "$ESTIMATION METHOD=1 MAXEVAL=99"
)
.with_name <- function(old, new) gsub(old, new, .reserved_stream, fixed = TRUE)

test_that(".nm_assigned_names reads assignments, not comparisons, comments or DADT", {
  l <- c("$PK", "  IF (SEX.EQ.1) DF = 2 ; X = 5", "  A1 = 1", "  IF (A1 .GE. 1) THEN",
         "  B2 = 3", "  ENDIF", "$DES", "  DADT(1) = -K*A(1)", "  C = A1 == 1",
         "$THETA (0,1)", "  NOTCODE = 4")
  expect_setequal(.nm_assigned_names(l), c("DF", "A1", "B2", "C"))
})

test_that("only the names rxode2 itself rejects are renamed", {
  skip_if_not_installed("rxode2")
  expect_equal(.nm_reserved_renames(.gap_ctl("df.ctl", .reserved_stream)),
               c(DF = "DF_VAR"))
  # the control from the issue: same stream, a name rxode2 accepts
  expect_equal(.nm_reserved_renames(.gap_ctl("ke.ctl", .with_name("DF", "KE"))),
               character())
})

test_that("the rename touches code records only, whole words only, and no comment", {
  l <- c("$INPUT ID TIME DF", "$PK", "  DF = 1 ; DF stays in a comment", "  DFX = DF + ADF",
         "$DES", "  DADT(1) = -DF*A(1)", "$TABLE ID DF")
  out <- .nm_rename_identifiers(l, c(DF = "DF_VAR"))
  expect_identical(out[c(1, 7)], l[c(1, 7)])
  expect_equal(out[3], "  DF_VAR = 1 ; DF stays in a comment")
  expect_equal(out[4], "  DFX = DF_VAR + ADF")
  expect_equal(out[6], "  DADT(1) = -DF_VAR*A(1)")
})

test_that("a stream with nothing to rename is returned byte for byte", {
  l <- .with_name("DF", "KE")
  expect_identical(.nm_rename_identifiers(l, character()), l)
})

test_that("the new name does not collide with one already in the stream", {
  skip_if_not_installed("rxode2")
  p <- .gap_ctl("c.ctl", sub("  IPRED = F", "  DF_VAR = 1\n  IPRED = F",
                              .reserved_stream, fixed = TRUE))
  expect_equal(unname(.nm_reserved_renames(p)), "DF_VARX")
})

# -- Tier 2: through nm_to_ferx() -----------------------------------------------

test_that("a $DES intermediate named DF translates, renamed, with an INFO", {
  skip_if_not_installed("nonmem2rx")
  skip_if_not_installed("rxode2")
  res <- suppressMessages(suppressWarnings(
    nm_to_ferx(.gap_ctl("df.ctl", .reserved_stream), validate = FALSE)))
  expect_match(res$ferx_text, "DF_VAR = CL/V", fixed = TRUE)
  expect_match(res$ferx_text, "- DF_VAR * CENT", fixed = TRUE)
  expect_false(grepl("\\bDF\\b", gsub("#[^\n]*", "", res$ferx_text)))
  expect_equal(sum(grepl("^INFO.*DF was renamed to DF_VAR", res$warnings)), 1L)
  expect_length(res$unsupported, 0L)
})

test_that("the renamed model is accepted by the engine", {
  skip_if_not_installed("nonmem2rx")
  skip_if_not_installed("rxode2")
  skip_if_not_installed("ferx")
  res <- suppressMessages(suppressWarnings(
    nm_to_ferx(.gap_ctl("df.ctl", .reserved_stream), validate = FALSE)))
  tmp <- tempfile(fileext = ".ferx")
  writeLines(res$ferx_text, tmp)
  on.exit(unlink(tmp))
  out <- NULL
  utils::capture.output(out <- ferx::ferx_model_validate(tmp))
  expect_true(out$ok)
})

test_that("the control (KE) is not renamed and says nothing about it", {
  skip_if_not_installed("nonmem2rx")
  skip_if_not_installed("rxode2")
  res <- suppressMessages(suppressWarnings(
    nm_to_ferx(.gap_ctl("ke.ctl", .with_name("DF", "KE")), validate = FALSE)))
  expect_match(res$ferx_text, "KE = CL/V", fixed = TRUE)
  expect_length(grep("renamed", res$warnings), 0L)
})
