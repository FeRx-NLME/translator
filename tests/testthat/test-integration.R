# Integration tests -- require nonmem2rx and/or rxode2.
# All tests are gated with skip_if_not_installed(); they are skipped in CI
# unless the relevant packages are available.  On first run with the packages
# present, expect_snapshot() creates the _snaps/ files.  Subsequent runs
# compare against those snapshots.

# -- helpers ------------------------------------------------------------------

# Strip machine-specific installed paths from header comment so snapshots are
# portable across machines and CI.  Reduces e.g.
#   "# Translated from nonmem: /Library/.../1cpt_oral.ctl"
# to
#   "# Translated from nonmem: 1cpt_oral.ctl"
norm_snap <- function(txt) {
  sub("(# Translated from [^:]+: ).*/([^/\n]+)", "\\1\\2", txt, perl = TRUE)
}

# -- NONMEM models ------------------------------------------------------------

test_that("1-cpt oral NONMEM: snapshot + no unsupported", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("1cpt_oral.ctl")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "one_cpt_oral", fixed = TRUE)
})

test_that("2-cpt oral with covariates: snapshot + no unsupported", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("2cpt_oral_cov.ctl")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "two_cpt_oral", fixed = TRUE)
})

test_that("2-cpt IV: infers two_cpt_iv", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("2cpt_iv.ctl")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "two_cpt_iv", fixed = TRUE)
})

test_that("3-cpt IV: translates to three_cpt_iv pk macro", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("3cpt_iv.ctl")))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "three_cpt_iv", fixed = TRUE)
})

test_that("ODE warfarin: full $DES path, [odes] section present", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("ode_warfarin.ctl")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "[odes]",       fixed = TRUE)
  expect_match(result$ferx_text, "d/dt(DEPOT)",  fixed = TRUE)
  expect_length(result$unsupported, 0L)
})

test_that("block omega: block_omega line in output", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("block_omega.ctl")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "block_omega", fixed = TRUE)
})

test_that("IOV model: KAPPA_CL emitted as omega + flattening warning (nonmem2rx treats IOV as IIV)", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("iov.ctl")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "KAPPA_CL", fixed = TRUE)
  # nonmem2rx flattens the ETA-coded IOV to IIV; the translator must warn so the
  # silent loss of occasion structure is visible to the user.
  expect_true(any(grepl("inter-occasion", result$warnings, fixed = TRUE)))
})

# -- nlmixr2 models -----------------------------------------------------------

test_that("1-cpt oral nlmixr2: snapshot + one_cpt_oral", {
  skip_if_not_installed("rxode2")
  fn     <- source(r2_path("1cpt_oral_nlmixr2.R"))$value
  result <- suppressWarnings(nlmixr2_to_ferx(fn))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_length(result$unsupported, 0L)
  expect_match(result$ferx_text, "one_cpt_oral", fixed = TRUE)
})

test_that("nlmixr2 source: a KAPPA-named IIV eta does NOT trigger the NONMEM-only IOV warning", {
  skip_if_not_installed("rxode2")
  fn     <- source(r2_path("iov_kappa_nlmixr2.R"))$value
  # suppressWarnings() silences rxode2's benign "non-mu referenced" parse note
  # for exp(eta.cl + kappa.cl); it does not touch result$warnings (the
  # translator's own channel), which is what the assertions below check.
  result <- suppressWarnings(nlmixr2_to_ferx(fn))
  # The eta is present in the IIV block (so the helper would match its name)...
  expect_match(result$ferx_text, "omega KAPPA_CL", fixed = TRUE)
  # ...but the flattening warning is nonmem2rx-specific and must stay silent here.
  expect_false(any(grepl("inter-occasion", result$warnings, fixed = TRUE)))
})

test_that("ODE nlmixr2: d/dt expressions produce [odes] section", {
  skip_if_not_installed("rxode2")
  fn     <- source(r2_path("ode_nlmixr2.R"))$value
  result <- suppressWarnings(nlmixr2_to_ferx(fn))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "[odes]",      fixed = TRUE)
  expect_match(result$ferx_text, "d/dt(depot)", fixed = TRUE)
  expect_length(result$unsupported, 0L)
})

# -- output writing -----------------------------------------------------------

test_that("nm_to_ferx writes file when output path given", {
  skip_if_not_installed("nonmem2rx")
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))
  nm_to_ferx(nm_path("1cpt_oral.ctl"), output = path)
  expect_true(file.exists(path))
  expect_match(paste(readLines(path), collapse = "\n"), "[parameters]", fixed = TRUE)
})

# -- amp.sim example models ---------------------------------------------------

test_that("amp.sim 1-cpt oral ODE: [odes] section + obs_cmt inferred", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("pk_1cmt_oral.mod")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "[odes]",        fixed = TRUE)
  expect_match(result$ferx_text, "obs_cmt=",      fixed = TRUE)
  expect_match(result$ferx_text, "d/dt(",         fixed = TRUE)
  expect_match(result$ferx_text, "proportional",  fixed = TRUE)
  expect_length(result$unsupported, 0L)
})

test_that("pk_1cmt_oral.mod: S2=V scaling emits [scaling] obs_scale = V", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("pk_1cmt_oral.mod")))
  expect_match(result$ferx_text, "[scaling]",      fixed = TRUE)
  expect_match(result$ferx_text, "obs_scale = V",  fixed = TRUE)
  expect_true(any(grepl("S2 = V", result$warnings, fixed = TRUE)))
})

test_that("amp.sim PKPD indirect response: 4-state ODE + additive error", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("pkpd_ir.mod")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_match(result$ferx_text, "[odes]",   fixed = TRUE)
  expect_match(result$ferx_text, "obs_cmt=", fixed = TRUE)
  expect_match(result$ferx_text, "d/dt(",    fixed = TRUE)
  expect_match(result$ferx_text, "additive", fixed = TRUE)
  n_odes <- length(regmatches(result$ferx_text,
                              gregexpr("d/dt\\(", result$ferx_text))[[1]])
  expect_equal(n_odes, 4L)
})

test_that("pk_1cmt_oral_ampsim: fixed-effect V passthrough appears in pk macro", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("pk_1cmt_oral_ampsim.ctl")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  # Fixed-effect V (no ETA) must appear as passthrough in [individual_parameters]
  # and be passed to the pk macro, otherwise ferx predicts zero concentration.
  # The passthrough reads `V = TVV`, not `V = V`: the source names the theta V
  # too, and a theta silently shadows an identically named individual parameter,
  # so the theta is renamed.
  expect_match(result$ferx_text, "V = TVV",            fixed = TRUE)
  expect_no_match(result$ferx_text, "V = V\n",        perl  = TRUE)
  expect_match(result$ferx_text, "v=V",                fixed = TRUE)
  expect_match(result$ferx_text, "one_cpt_oral",       fixed = TRUE)
  expect_length(result$unsupported, 0L)
})

test_that("TMDD with a FLAG dispatch: Form C readout + selected error model", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("qss_tmdd.mod")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_length(result$unsupported, 0L)
  # The three halves of defect 5, asserted rather than left to the snapshot: a
  # snapshot records whatever is emitted, including a silently dropped dispatch.
  expect_match(result$ferx_text, "y = if (FLAG == 2) c_RTOT else CENT/VC",
               fixed = TRUE)
  expect_match(result$ferx_text, "if (FLAG == 2) { DV ~ proportional(EPS2) }",
               fixed = TRUE)
  expect_match(result$ferx_text, "else { DV ~ proportional(EPS1) }", fixed = TRUE)
  # Defects 12 and 15: both inferences the readout replaces are gone. obs_scale
  # would double-scale every prediction, obs_cmt is ignored.
  expect_no_match(result$ferx_text, "obs_scale", fixed = TRUE)
  expect_no_match(result$ferx_text, "obs_cmt", fixed = TRUE)
  # Defect 6: the indicators are consumed, not left dead in the parameter block.
  expect_no_match(result$ferx_text, "W1", fixed = TRUE)
})

test_that("PKPD with a CMT dispatch: per-CMT readout + per-CMT error model", {
  skip_if_not_installed("nonmem2rx")
  result <- suppressWarnings(nm_to_ferx(nm_path("pkpd_cmt.mod")))
  expect_snapshot(cat(norm_snap(result$ferx_text)))
  expect_length(result$unsupported, 0L)
  # ferx does not expose CMT as a covariate, so this source cannot use Form C --
  # measured, `Model references covariate(s) not found in data: CMT`.
  expect_match(result$ferx_text, "y[CMT=1] = CENT/VC", fixed = TRUE)
  expect_match(result$ferx_text, "y[CMT=2] = PD", fixed = TRUE)
  expect_match(result$ferx_text, "CMT=1: DV ~ proportional(EPS1)", fixed = TRUE)
  # The two endpoints carry different error TYPES, which ferx allows per branch.
  expect_match(result$ferx_text, "CMT=2: DV ~ additive(EPS2)", fixed = TRUE)
  expect_no_match(result$ferx_text, "obs_cmt", fixed = TRUE)
})

# -- corpus-wide sweep (engine-free, runs on every PR) ------------------------

# Translates every bundled model once and asserts three things: that nothing
# crashes, that no theta shadows an individual parameter, and that every emitted
# name is a legal ferx identifier.
#
# The first two the engine cannot check for us. A theta sharing a name with an
# individual parameter shadows it in every ferx block where thetas are in scope,
# so the individual definition is written and never read -- silently, with no
# diagnostic from ferx_model_validate(). $unsupported is printed as a gap table,
# not failed.
#
# The identifier check the engine WOULD catch, but only for whichever model
# happens to be validated, and only as an E_PARSE that says nothing about which
# name is at fault. Asserting it here names the model and the identifier, and it
# runs without ferx installed. Covariates are deliberately excluded: they must
# keep the case and spelling of the data column, so an illegal one is reported
# as unsupported rather than renamed.
#
# The engine half of this sweep (validate each emitted .ferx, fail on
# error-severity diagnostics) lives in test-concordance.R behind the ferx gate.
test_that("every bundled model translates, without a shadowing theta", {
  skip_if_not_installed("nonmem2rx")

  models <- .bundled_nm_models()
  expect_gt(length(models), 0L)

  crashed   <- character()
  offenders <- character()
  illegal   <- character()
  gaps      <- character()

  for (m in models) {
    # Through the SAME path a user gets, hints and all. Calling rxui_to_ir()
    # bare swept every model through the pre-fix guessing path: defobs_not_last.ctl,
    # the fixture bundled specifically to guard DEFOBS, came out with
    # obs_cmt=PERIPH and no [scaling] here while nm_to_ferx() gave obs_cmt=CENT
    # and obs_scale=V. The gate was inspecting output no user could obtain.
    ir <- tryCatch(suppressWarnings(
                     rxui_to_ir(nonmem2rx::nonmem2rx(m), source_format = "nonmem",
                                scaling_hint = .extract_nm_scaling(m),
                                obs_hint     = .extract_nm_defobs(m))),
                   error = function(e) conditionMessage(e))
    if (is.character(ir)) {
      crashed <- c(crashed, paste0(basename(m), " -- ", ir))
      next
    }
    clash <- intersect(
      toupper(vapply(ir$thetas, function(t) t$name, "")),
      toupper(vapply(ir$indiv_params, function(p) p$lhs, ""))
    )
    if (length(clash) > 0)
      offenders <- c(offenders, paste0(basename(m), ": ", paste(clash, collapse = ", ")))

    # Tokenise the EMITTED TEXT, not the IR's declaration fields. An illegal
    # name is accepted where it is declared -- `ode(states=[c.RTOT])` parses --
    # and rejected at every reference, so a declaration-only check inspects the
    # half that works. Enumerating IR fields also silently misses whichever
    # channel is added next: reading the artefact cannot.
    # Guarded, because emit_ferx() opens with validate_ferx_ir() and that
    # cli_abort()s: uncaught, ONE bad model ends the sweep at that model and the
    # crashed/offenders/illegal/gaps tables -- all message()d after the loop --
    # are never printed, which is the opposite of a sweep that reports per-model
    # gaps. Classed rather than bare, unlike the rxui_to_ir() guard above: there
    # success is a list and failure a string, so is.character() tells them
    # apart; here both are one string.
    txt <- tryCatch(emit_ferx(ir),
                    error = function(e) structure(conditionMessage(e),
                                                  class = "emit_failure"))
    if (inherits(txt, "emit_failure")) {
      crashed <- c(crashed, paste0(basename(m), " -- emit_ferx: ", unclass(txt)))
      next
    }
    lines <- strsplit(txt, "\n")[[1]]
    lines <- lines[!grepl("^\\s*#", lines)]        # comments carry source names
    # TWO checks, because they can see different things.
    #
    # (a) Declaration positions, extracted STRUCTURALLY. The token scan below
    # cannot do this job: its pattern admits `.` as the only illegal
    # continuation character, so every OTHER illegal character simply splits a
    # bad name into two legal tokens -- `ode(obs_cmt=A-B, states=[A-B])` scans
    # as `ode, obs_cmt, A, B, states`, all legal, and the assertion passes on an
    # unparseable file. `has space` and `2CPT` fail the same way (a token cannot
    # begin with a digit, so `2CPT` is only ever seen as `CPT`). Reading the
    # declaration positions by shape catches all of them.
    decl <- c(
      unlist(lapply(regmatches(lines, gregexpr("states=\\[[^]]*\\]", lines)),
                    function(m) unlist(strsplit(sub("^states=\\[(.*)\\]$", "\\1", m), ",")))),
      unlist(regmatches(lines, gregexpr("(?<=obs_cmt=)[^,)]+", lines, perl = TRUE))),
      unlist(regmatches(lines, gregexpr("(?<=d/dt\\()[^)]+", lines, perl = TRUE))))
    decl <- trimws(decl)
    decl <- decl[nzchar(decl)]
    #
    # (b) The token scan, which still earns its place for EXPRESSION text: a
    # dotted covariate reference (`WT.KG`) is the realistic illegal spelling and
    # is invisible to (a). Declaration legality itself is now guaranteed
    # structurally by validate_ferx_ir(), which emit_ferx() always runs -- (a)
    # is the belt to that braces, and it reads the artefact rather than the IR.
    toks  <- unique(unlist(regmatches(
      lines, gregexpr("[A-Za-z_][A-Za-z0-9_.]*", lines))))
    bad <- unique(c(decl[!.is_ferx_ident(decl)], toks[!.is_ferx_ident(toks)]))
    # No setdiff for "true"/"false" here: .is_ferx_ident("true") is already TRUE,
    # so filtering them was dead code that read as though a real exemption
    # existed.
    #
    # A name the translator already reported as untranslatable is not a silent
    # leak -- covariates in particular must keep the data column's exact
    # spelling, so an illegal one is declared unsupported rather than renamed.
    # Matched as a whole token, NOT as a substring of any unsupported string:
    # `grepl(b, ir$unsupported, fixed = TRUE)` exempted a name because some
    # unrelated entry happened to contain those characters, which already
    # green-lit a genuinely unparseable `ode(states=[c.RTOT])`.
    # Only names reported as an untranslatable COVARIATE are exempt, and the
    # name is taken from that entry's fixed prefix rather than by tokenising the
    # whole string. Harvesting every token from every $unsupported entry is
    # self-defeating now that entries can carry deparsed source: a model with
    # `q(A.B) <- ...` yields "unsupported assignment target: q(A.B)", whose
    # tokens include `A.B`, which then exempted `A.B` everywhere else in the
    # file -- including a genuinely unparseable `ode(states=[A.B])`.
    cov_pref   <- "covariate name is not a legal ferx identifier: "
    exempt     <- sub(cov_pref, "",
                      grep(cov_pref, ir$unsupported, fixed = TRUE, value = TRUE),
                      fixed = TRUE)
    bad <- setdiff(bad, exempt)
    if (length(bad) > 0)
      illegal <- c(illegal, paste0(basename(m), ": ", paste(bad, collapse = ", ")))

    if (length(ir$unsupported) > 0)
      gaps <- c(gaps, paste0(basename(m), " -- ", ir$unsupported))
  }

  if (length(gaps) > 0)
    message("\ntranslation gap report (", length(gaps), " gap(s) across ",
            length(models), " models):\n", paste(gaps, collapse = "\n"))
  else
    message("translation gap report: no unsupported features across ",
            length(models), " models")

  expect_equal(crashed, character())
  expect_equal(offenders, character())
  expect_equal(illegal, character())
})

test_that("a covariate conditional in $PK does not block a $ERROR dispatch", {
  skip_if_not_installed("nonmem2rx")
  # The ordinary shape of a PKPD model: a covariate effect on a PK parameter
  # AND a two-endpoint error model. Both are conditionals, and the readout
  # divides by the parameter the covariate acts on, so the $PK conditional sits
  # inside the readout's backward closure. Treating every conditional there as
  # part of the dispatch reads SEX beside FLAG and fails the whole model on
  # "more than one column"; only the FLAG conditionals select an endpoint.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "cov_disp.mod")
  writeLines(c(
    "$PROBLEM covariate effect plus endpoint dispatch",
    "$INPUT ID TIME AMT EVID MDV CMT SEX FLAG DV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINES ADVAN13 TOL=9",
    "$MODEL",
    "  COMP=(CENT, DEFDOSE, DEFOBS)",
    "  COMP=(PERI)",
    "$PK",
    "  KEL = THETA(1)*EXP(ETA(1))",
    "  VC  = THETA(2)",
    "  IF (SEX.EQ.1) VC = THETA(2)*THETA(3)",
    "  KCP = THETA(4)",
    "$DES",
    "  DADT(1) = -KEL*A(1) - KCP*A(1)",
    "  DADT(2) =  KCP*A(1)",
    "$ERROR",
    "  CONC = A(1)/VC",
    "  PERIF = A(2)",
    "  IPRED = CONC",
    "  IF (FLAG.EQ.2) IPRED = PERIF",
    "  W1 = 0",
    "  W2 = 0",
    "  IF (FLAG.EQ.1) W1 = 1",
    "  IF (FLAG.EQ.2) W2 = 1",
    "  Y = IPRED*(1 + W1*EPS(1) + W2*EPS(2))",
    "$THETA",
    "  (0.001, 0.1, 1.0)",
    "  (0.5, 5.0, 50.0)",
    "  (0.1, 1.2, 3.0)",
    "  (0.01, 0.3, 5.0)",
    "$OMEGA",
    "  0.09",
    "$SIGMA",
    "  0.0225",
    "  0.04",
    "$ESTIMATION METHOD=1 INTER MAXEVAL=9999"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  expect_match(result$ferx_text, "y = if (FLAG == 2) PERI else CENT/VC",
               fixed = TRUE)
  expect_match(result$ferx_text, "if (FLAG == 2) { DV ~ proportional(EPS2) }",
               fixed = TRUE)
  # The covariate conditional stays where it belongs, referencing the theta it
  # was written against, and `y` reads VC as an ordinary individual parameter.
  expect_match(result$ferx_text, "if (SEX == 1) { VC = THETA2 * THETA3 }",
               fixed = TRUE)
  expect_match(result$ferx_text, "  VC = THETA2\n", fixed = TRUE)
  expect_no_match(result$ferx_text, "more than one column", fixed = TRUE)
})

test_that("a non-dispatch conditional on the readout is reported, not silently dropped", {
  skip_if_not_installed("nonmem2rx")
  # A clamp on the prediction. It defines a name the readout must resolve per
  # case, and it is not an equality test, so no dispatch can be built from it.
  # Before this phase it reached the generic "nothing reads it, so it has no
  # effect and is dropped" INFO -- which is wrong twice over: it does have an
  # effect in NONMEM, and it is not dead. Nothing else about the output changes,
  # so the report is the whole of the fix.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "clamped.mod")
  writeLines(c(
    "$PROBLEM clamped readout",
    "$INPUT ID TIME AMT EVID MDV CMT FLAG DV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINES ADVAN13 TOL=9",
    "$MODEL",
    "  COMP=(CENT, DEFDOSE, DEFOBS)",
    "  COMP=(PERI)",
    "$PK",
    "  KEL = THETA(1)*EXP(ETA(1))",
    "  VC  = THETA(2)",
    "  KCP = THETA(3)",
    "$DES",
    "  DADT(1) = -KEL*A(1) - KCP*A(1)",
    "  DADT(2) =  KCP*A(1)",
    "$ERROR",
    "  CONC = A(1)/VC",
    "  IPRED = CONC",
    "  IF (CONC.LT.0.0) IPRED = 0",
    "  W1 = 0",
    "  IF (FLAG.EQ.1) W1 = 1",
    "  Y = IPRED*(1 + W1*EPS(1))",
    "$THETA",
    "  (0.001, 0.1, 1.0)",
    "  (0.5, 5.0, 50.0)",
    "  (0.01, 0.3, 5.0)",
    "$OMEGA",
    "  0.09",
    "$SIGMA",
    "  0.0225",
    "$ESTIMATION METHOD=1 INTER MAXEVAL=9999"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  blocked <- grep("^ERROR .*\\$ERROR block is conditional", result$warnings, value = TRUE)
  expect_length(blocked, 1L)
  expect_match(blocked, "CONC < 0", fixed = TRUE)
  expect_match(blocked, "not a `<column> == <number>` test", fixed = TRUE)
  # Reported, not acted on: no readout is emitted and the pre-6b path runs, so
  # the engine still gets the file it would have got before. Anchored on the
  # section header, because the ERROR text itself reaches the file as a
  # `# WARNING:` comment and names `[scaling]` inside it.
  expect_no_match(result$ferx_text, "\n[scaling]\n", fixed = TRUE)
  expect_match(result$ferx_text, "ode(obs_cmt=CENT, states=[CENT, PERI])",
               fixed = TRUE)
})

test_that("an unexpressible endpoint reports once, names the right sigma, and keeps the readout", {
  skip_if_not_installed("nonmem2rx")
  # A clean FLAG dispatch where only the FLAG == 2 endpoint carries a scaled
  # sigma. Before this the model fell back to the single-endpoint path and
  # produced TWO ERRORs: the per-case one, and a second that re-diagnosed the
  # un-substituted `Y` and reported "`EPS1` is weighted by `IPRED * W1` ... an
  # indicator-weighted ... error model". EPS1 is fine, EPS2 is the problem, and
  # "indicator-weighted" is what this phase had just resolved. The file also
  # carried a single-endpoint `DV ~ combined(EPS1, EPS2)` suggestion, which is
  # the wrong SHAPE for a model that needs one entry per endpoint.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "gap.mod")
  writeLines(c(
    "$PROBLEM one endpoint not expressible",
    "$INPUT ID TIME AMT EVID MDV CMT FLAG DV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINES ADVAN13 TOL=9",
    "$MODEL",
    "  COMP=(CENT, DEFDOSE, DEFOBS)",
    "  COMP=(PERI)",
    "$PK",
    "  KEL = THETA(1)*EXP(ETA(1))",
    "  VC  = THETA(2)",
    "  KCP = THETA(3)",
    "$DES",
    "  DADT(1) = -KEL*A(1) - KCP*A(1)",
    "  DADT(2) =  KCP*A(1)",
    "$ERROR",
    "  CONC = A(1)/VC",
    "  PERIF = A(2)",
    "  IPRED = CONC",
    "  IF (FLAG.EQ.2) IPRED = PERIF",
    "  W1 = 0",
    "  W2 = 0",
    "  IF (FLAG.EQ.1) W1 = 1",
    "  IF (FLAG.EQ.2) W2 = 1",
    "  Y = IPRED*(1 + W1*EPS(1)) + W2*EPS(2)*THETA(3)",
    "$THETA",
    "  (0.001, 0.1, 1.0)",
    "  (0.5, 5.0, 50.0)",
    "  (0.01, 0.3, 5.0)",
    "$OMEGA",
    "  0.09",
    "$SIGMA",
    "  0.0225",
    "  0.04",
    "$ESTIMATION METHOD=1 INTER MAXEVAL=9999"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  errs <- grep("^ERROR", result$warnings, value = TRUE)
  expect_length(errs, 1L)
  expect_match(errs, "the endpoint for FLAG == 2", fixed = TRUE)
  expect_match(errs, "`EPS2` is weighted by `THETA3`", fixed = TRUE)
  # The prediction named is eps-free.
  expect_match(errs, "the prediction `PERI`", fixed = TRUE)
  # The epsilon that is FINE is not blamed, and the indicator complaint this
  # phase resolved is not raised.
  expect_no_match(errs, "EPS1` is weighted", fixed = TRUE)

  # The readout is emitted for real and is complete.
  expect_match(result$ferx_text, "\n[scaling]\n  y = if (FLAG == 2) PERI else CENT/VC",
               fixed = TRUE)
  # The error model is a dispatch-shaped comment with one marked gap.
  expect_match(result$ferx_text, "#   if (FLAG == 2) { DV ~ ??? }", fixed = TRUE)
  expect_match(result$ferx_text, "#   else { DV ~ proportional(EPS1) }", fixed = TRUE)
  expect_no_match(result$ferx_text, "combined(EPS1, EPS2)", fixed = TRUE)
  # ferx still rejects the file, which is what keeps this loud.
  expect_no_match(result$ferx_text, "\n[error_model]\n", fixed = TRUE)
})

test_that("a per-CMT gap emits keyed y[CMT=N] plus a keyed suggestion, wherever it falls", {
  skip_if_not_installed("nonmem2rx")
  # The per-CMT partial path had no pipeline coverage at all, and the unit
  # fixtures for it both put the failing endpoint FIRST. Here CMT == 2 is seen
  # first, so the broken endpoint (CMT == 1, a scaled sigma) is the LAST listed
  # value -- the position that used to send the whole model back to the
  # single-endpoint path and bring back both the two-ERROR report and a
  # `DV ~ combined(...)` suggestion of the wrong shape.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "pcgap.mod")
  writeLines(c(
    "$PROBLEM per-CMT dispatch, last listed endpoint unexpressible",
    "$INPUT ID TIME AMT EVID MDV CMT DV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINES ADVAN13 TOL=9",
    "$MODEL",
    "  COMP=(CENT, DEFDOSE, DEFOBS)",
    "  COMP=(PD)",
    "$PK",
    "  KEL = THETA(1)*EXP(ETA(1))",
    "  VC  = THETA(2)",
    "  KCP = THETA(3)",
    "$DES",
    "  DADT(1) = -KEL*A(1) - KCP*A(1)",
    "  DADT(2) =  KCP*A(1)",
    "$ERROR",
    "  CONC = A(1)/VC",
    "  RESP = A(2)",
    "  IPRED = CONC",
    "  IF (CMT.EQ.2) IPRED = RESP",
    "  W1 = 0",
    "  W2 = 0",
    "  IF (CMT.EQ.1) W1 = 1",
    "  IF (CMT.EQ.2) W2 = 1",
    "  Y = IPRED*(1 + W2*EPS(2)) + W1*EPS(1)*THETA(3)",
    "$THETA",
    "  (0.001, 0.1, 1.0)",
    "  (0.5, 5.0, 50.0)",
    "  (0.01, 0.3, 5.0)",
    "$OMEGA",
    "  0.09",
    "$SIGMA",
    "  0.0225",
    "  0.04",
    "$ESTIMATION METHOD=1 INTER MAXEVAL=9999"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  errs <- grep("^ERROR", result$warnings, value = TRUE)
  expect_length(errs, 1L)
  expect_match(errs, "the endpoint for CMT == 1", fixed = TRUE)
  expect_match(errs, "`EPS1` is weighted by `THETA3`", fixed = TRUE)

  # Both readouts emitted and keyed; ferx checks per-CMT coverage against the
  # data, so an incomplete set is rejected by name rather than silently.
  expect_match(result$ferx_text, "  y[CMT=1] = CENT/VC", fixed = TRUE)
  expect_match(result$ferx_text, "  y[CMT=2] = PD", fixed = TRUE)
  # The suggestion keeps its CMT keys and marks only the endpoint that failed.
  expect_match(result$ferx_text, "#   CMT=1: DV ~ ???", fixed = TRUE)
  expect_match(result$ferx_text, "#   CMT=2: DV ~ proportional(EPS2)", fixed = TRUE)
  expect_no_match(result$ferx_text, "combined(", fixed = TRUE)
  # The indicators are consumed, not left dead in the parameter block.
  expect_no_match(result$ferx_text, "W1", fixed = TRUE)
  expect_no_match(result$ferx_text, "obs_cmt", fixed = TRUE)
})

# -- Phase 6c: which compartment is observed ----------------------------------
#
# These are tier 2 rather than tier 1 because the evidence they weigh -- $MODEL
# DEFOBS and $PK's `S<n>` -- is read off the raw control stream by nm_to_ferx()
# and reaches rxui_to_ir() only as hints. Reaching the lower tiers from an
# inline rxode2 model is not possible either: rxode2 requires the endpoint to
# name a defined variable, and any such variable resolves to a compartment at
# tier 1, so the cascade never gets past its first step.

test_that("the observed compartment is taken from S<n> when nothing else names it", {
  skip_if_not_installed("nonmem2rx")
  # s_scaling_not_last.ctl declares PERIPH last and scales CENTRAL. Before this
  # phase the cascade ran out of evidence and took `tail(states)`, which is
  # PERIPH -- and then dropped [scaling] as well, because the scaled
  # compartment was not the one it had decided was observed.
  m <- system.file("testmodels/nonmem/s_scaling_not_last.ctl",
                   package = "ferxtranslate")
  skip_if(m == "")
  result <- suppressWarnings(nm_to_ferx(m, validate = FALSE))

  expect_match(result$ferx_text,
               "ode(obs_cmt=CENTRAL, states=[DEPOT, CENTRAL, PERIPH])",
               fixed = TRUE)
  # Both halves, because the failure took both: the wrong compartment AND the
  # scaling that went missing with it.
  expect_match(result$ferx_text, "obs_scale = V", fixed = TRUE)
  expect_length(result$unsupported, 0L)
  expect_length(grep("^ERROR .*compartment could be inferred", result$warnings), 0L)
})

test_that("$MODEL DEFOBS outranks $PK scaling when the two disagree", {
  skip_if_not_installed("nonmem2rx")
  # DEFOBS states which compartment is observed; `S<n>` says which compartment's
  # amount is converted to the data's scale. The first is a statement and the
  # second an inference from purpose, so DEFOBS wins. The fixture makes them
  # name different compartments, which is the only arrangement that can show the
  # order being wrong.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "defobs_beats_scaling.ctl")
  writeLines(c(
    "$PROBLEM DEFOBS and S2 disagree",
    "$INPUT ID TIME DV AMT EVID MDV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(DEPOT)",
    "  COMP=(CENTRAL)",
    "  COMP=(PERIPH, DEFOBS)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    "  Q  = THETA(4)",
    "  V3 = THETA(5)",
    "  S2 = V",
    "$DES",
    "  DADT(1) = -KA*A(1)",
    "  DADT(2) =  KA*A(1) - (CL/V)*A(2) - (Q/V)*A(2) + (Q/V3)*A(3)",
    "  DADT(3) =  (Q/V)*A(2) - (Q/V3)*A(3)",
    "$ERROR",
    "  IPRED = F",
    "  Y     = IPRED*(1 + EPS(1))",
    "$THETA (0,5) (0,50) (0,1) (0,8) (0,60)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  expect_match(result$ferx_text, "obs_cmt=PERIPH", fixed = TRUE)
  # The S<n> tier did not fire, so it must not have announced that it did.
  expect_length(grep("taken from \\$PK's S", result$warnings), 0L)
})

test_that("scaling for more than one compartment identifies none", {
  skip_if_not_installed("nonmem2rx")
  # A source that scales several compartments has named none of them, so the
  # tier declines rather than taking the lowest number. PERIPH is declared last
  # so declining is visible: the cascade falls through to the guess and says so,
  # where picking S1 would have answered DEPOT silently.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "two_scalings.ctl")
  writeLines(c(
    "$PROBLEM S1 and S2 both present",
    "$INPUT ID TIME DV AMT EVID MDV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(DEPOT)",
    "  COMP=(CENTRAL)",
    "  COMP=(PERIPH)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    "  Q  = THETA(4)",
    "  V3 = THETA(5)",
    "  S1 = V",
    "  S2 = V",
    "$DES",
    "  DADT(1) = -KA*A(1)",
    "  DADT(2) =  KA*A(1) - (CL/V)*A(2) - (Q/V)*A(2) + (Q/V3)*A(3)",
    "  DADT(3) =  (Q/V)*A(2) - (Q/V3)*A(3)",
    "$ERROR",
    "  IPRED = F",
    "  Y     = IPRED*(1 + EPS(1))",
    "$THETA (0,5) (0,50) (0,1) (0,8) (0,60)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  expect_length(grep("taken from \\$PK's S", result$warnings), 0L)
  expect_length(grep("^ERROR .*compartment could be inferred", result$warnings), 1L)
  expect_no_match(result$ferx_text, "obs_cmt=DEPOT", fixed = TRUE)
})

test_that("a source that names no compartment reports an ERROR and an unsupported entry", {
  skip_if_not_installed("nonmem2rx")
  # No compartment in the DV expression, no DEFOBS, no scaling. The answer is
  # declaration order, which is position and not evidence -- so it is reported
  # at ERROR and listed as a gap, not announced as an inference.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "no_evidence.ctl")
  writeLines(c(
    "$PROBLEM nothing names the observed compartment",
    "$INPUT ID TIME DV AMT EVID MDV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(DEPOT)",
    "  COMP=(CENTRAL)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    "$DES",
    "  DADT(1) = -KA*A(1)",
    "  DADT(2) =  KA*A(1) - (CL/V)*A(2)",
    "$ERROR",
    "  IPRED = F",
    "  Y     = IPRED*(1 + EPS(1))",
    "$THETA (0,5) (0,50) (0,1)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  err <- grep("^ERROR .*compartment could be inferred", result$warnings, value = TRUE)
  expect_length(err, 1L)
  expect_match(err, "declared last, which is position and not evidence", fixed = TRUE)
  expect_length(grep("obs_cmt guessed", result$unsupported), 1L)
  # Reported, not fatal: the file is still emitted and still names a compartment.
  expect_match(result$ferx_text, "obs_cmt=CENTRAL", fixed = TRUE)
})

test_that("a one-compartment model is not guessing and says nothing", {
  skip_if_not_installed("nonmem2rx")
  # `tail(states)` and "the only compartment" are the same answer here, so there
  # is no ambiguity to report. Without this carve-out the commonest ODE shape in
  # pharmacometrics collects an ERROR about a choice it never had.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "single_cmt.ctl")
  writeLines(c(
    "$PROBLEM one compartment, no DEFOBS, no scaling",
    "$INPUT ID TIME DV AMT EVID MDV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(CENTRAL)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "$DES",
    "  DADT(1) = -(CL/V)*A(1)",
    "$ERROR",
    "  IPRED = F",
    "  Y     = IPRED*(1 + EPS(1))",
    "$THETA (0,5) (0,50)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  expect_match(result$ferx_text, "ode(obs_cmt=CENTRAL, states=[CENTRAL])", fixed = TRUE)
  expect_length(grep("compartment could be inferred", result$warnings), 0L)
  expect_length(result$unsupported, 0L)
})

test_that("a reordered $DES is renumbered into $MODEL order, not declined", {
  skip_if_not_installed("nonmem2rx")
  # `n` in `S<n>` is a $MODEL COMP ordinal; d/dt order is not. nonmem2rx keeps
  # $DES statement order, so a block writing DADT(2) first yields states
  # [CENTRAL, DEPOT] while S2 still means CENTRAL.
  #
  # Phase 6c DECLINED here, because it had no way to reconcile the two. Issue
  # #25 supplies one: `states=[...]` is what ferx numbers compartments by, so
  # emitting it in $MODEL COMP order makes S2 mean CENTRAL again -- and makes
  # the source's own CMT column keep selecting the compartments it meant.
  # Measured on ferx 0.3.0: before the renumbering this file's CMT=1 dose landed
  # in CENTRAL rather than DEPOT, and predictions moved ~40% at t=1.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "des_reordered.ctl")
  writeLines(c(
    "$PROBLEM DADT written in reverse index order",
    "$INPUT ID TIME DV AMT EVID MDV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(DEPOT)",
    "  COMP=(CENTRAL)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    "  S2 = V",
    "$DES",
    "  DADT(2) =  KA*A(1) - (CL/V)*A(2)",
    "  DADT(1) = -KA*A(1)",
    "$ERROR",
    "  IPRED = F",
    "  Y     = IPRED*(1 + EPS(1))",
    "$THETA (0,5) (0,50) (0,1)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  # states=[...] carries NONMEM's numbering, NOT the $DES statement order.
  expect_match(result$ferx_text, "states=[DEPOT, CENTRAL]", fixed = TRUE)
  # The [odes] block keeps source order -- deliberately. ferx has no
  # use-before-def check, so an intermediate moved below the d/dt line that
  # reads it would silently read a stale slot; see .emit_odes_section().
  expect_lt(regexpr("d/dt(CENTRAL)", result$ferx_text, fixed = TRUE),
            regexpr("d/dt(DEPOT)",   result$ferx_text, fixed = TRUE))
  # S2 resolves to CENTRAL again, so the scaling is emitted rather than dropped.
  expect_match(result$ferx_text, "obs_cmt=CENTRAL", fixed = TRUE)
  expect_match(result$ferx_text, "obs_scale = V",   fixed = TRUE)
  # Renumbering is a correctness change to the emitted file, so it is announced.
  ren <- grep("^INFO .*put in \\$MODEL compartment order", result$warnings,
              value = TRUE)
  expect_length(ren, 1L)
  expect_match(ren, "DEPOT, CENTRAL", fixed = TRUE)
  # Nothing is declined any more, so neither the 6c disagreement WARN nor the
  # positional fallback fires.
  expect_length(grep("orderings disagree", result$warnings), 0L)
  expect_length(grep("obs_cmt guessed", result$unsupported), 0L)
})

test_that("two S<n> entries bind to the right one when $DES is reordered", {
  skip_if_not_installed("nonmem2rx")
  # Issue #25's headline case, and the reason it needs TWO scaling entries.
  #
  # $ERROR names A(2) outright, so obs_cmt resolves at the top of the cascade and
  # is CENTRAL either way -- this is not an obs_cmt defect. The defect is that
  # obs_cmt_num was a d/dt POSITION handed to a lookup keyed by NONMEM
  # compartment NUMBER: CENTRAL sits at position 1 of [CENTRAL, DEPOT], so the
  # lookup read S1 and emitted the DEPOT's scale variable, announcing it at INFO
  # as a success. Every prediction was then off by V/VD, here a factor of ~7.1.
  #
  # With only S2 present the same bug drops the scaling LOUDLY instead of
  # substituting the wrong variable, so a one-scaling fixture cannot tell "wrong
  # variable" from "no variable" and proves nothing. VD must differ from V.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "two_scalings_reordered.ctl")
  writeLines(c(
    "$PROBLEM two S<n> entries with DADT written in reverse index order",
    "$INPUT ID TIME DV AMT EVID MDV CMT",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(DEPOT)",
    "  COMP=(CENTRAL)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    "  VD = THETA(4)",
    "  S1 = VD",
    "  S2 = V",
    "$DES",
    "  DADT(2) =  KA*A(1) - (CL/V)*A(2)",
    "  DADT(1) = -KA*A(1)",
    "$ERROR",
    # `A(2)/S2`, not a bare `A(2)`. Both name compartment 2 outright, so either
    # resolves obs_cmt through the explicit tier, which is what this test is
    # about -- but only the dividing spelling is SCALED. NONMEM applies S<n> to
    # `F` and not to a bare `A(n)` (anchored in tests/nonmem-anchor/), so with
    # the bare spelling the correct output has no [scaling] block at all and
    # there is no binding left to get right. See issue #32.
    "  Y = A(2)/S2*(1 + EPS(1))",
    "$THETA (0,5) (0,50) (0,1) (0,7)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  expect_match(result$ferx_text, "obs_cmt=CENTRAL", fixed = TRUE)
  # The observed compartment is CENTRAL and NONMEM scales it by S2 = V. VD is
  # the depot's, and is what the defect emitted.
  expect_match(result$ferx_text, "obs_scale = V\n", perl = TRUE)
  expect_no_match(result$ferx_text, "obs_scale = VD", fixed = TRUE)
  # Announced against the right compartment number, too.
  expect_length(grep("INFO .*S2 = V detected", result$warnings), 1L)
  expect_length(grep("S1 = VD detected", result$warnings), 0L)
})

test_that("S<n> is applied to F but not to a bare A(n)", {
  skip_if_not_installed("nonmem2rx")
  # Anchored against NONMEM 7.6.0 (tests/nonmem-anchor/): two control streams
  # differing only in $ERROR, with S2 = V and V = 10, give a ratio of exactly 10
  # at every timepoint. `F` is A(n)/S<n>; a bare `A(n)` is the raw amount and
  # S<n> does not touch it.
  #
  # Three shapes, and the middle one is the defect: emitting obs_scale there
  # divides by a scale NONMEM never applied, putting every prediction low by a
  # factor of S<n>, silently, in a file that validates and fits.
  shapes <- list(
    list(err = "Y = F*(1 + EPS(1))",       scaled = TRUE,  what = "F"),
    list(err = "Y = A(2)*(1 + EPS(1))",    scaled = FALSE, what = "bare A(n)"),
    list(err = "Y = A(2)/S2*(1 + EPS(1))", scaled = TRUE,  what = "A(n)/S<n>"))

  dir <- tmp_ctl_dir()
  for (sh in shapes) {
    ctl <- file.path(dir, paste0("shape_", make.names(sh$what), ".ctl"))
    writeLines(c(
      "$PROBLEM S<n> application by DV shape",
      "$INPUT ID TIME DV AMT EVID MDV CMT",
      "$DATA d.csv IGNORE=@",
      "$SUBROUTINE ADVAN6 TOL=6",
      "$MODEL",
      "  COMP=(DEPOT, DEFDOSE)",
      "  COMP=(CENTRAL, DEFOBS)",
      "$PK",
      "  KA = THETA(1)",
      "  CL = THETA(2)",
      "  V  = THETA(3)",
      "  S2 = V",
      "$DES",
      "  DADT(1) = -KA*A(1)",
      "  DADT(2) =  KA*A(1) - (CL/V)*A(2)",
      "$ERROR",
      paste0("  ", sh$err),
      "$THETA (0,1) (0,3) (0,10)",
      "$OMEGA 0.09",
      "$SIGMA 0.04",
      "$EST METHOD=1"), ctl)

    result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
    body   <- grep("^#", strsplit(result$ferx_text, "\n")[[1]],
                   value = TRUE, invert = TRUE)
    has_scaling <- any(grepl("obs_scale", body, fixed = TRUE))
    expect_equal(has_scaling, sh$scaled,
                 label = paste0("obs_scale emitted for ", sh$what))
    # Silence is not acceptable for the unscaled case: the source declares S<n>,
    # so a reader diffing the two files will ask where the block went.
    if (!sh$scaled)
      expect_length(grep("reads the compartment amount", result$warnings), 1L)
  }
})

.scale_shape_ctl <- function(dir, name, pk_extra, err) {
  ctl <- file.path(dir, name)
  writeLines(c(
    "$PROBLEM scaling by DV shape",
    "$INPUT ID TIME DV AMT EVID MDV CMT",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(DEPOT, DEFDOSE)",
    "  COMP=(CENTRAL, DEFOBS)",
    "$PK",
    "  KA = THETA(1)",
    "  CL = THETA(2)",
    "  V  = THETA(3)",
    pk_extra,
    "$DES",
    "  DADT(1) = -KA*A(1)",
    "  DADT(2) =  KA*A(1) - (CL/V)*A(2)",
    "$ERROR",
    paste0("  ", err),
    "$THETA (0,1) (0,3) (0,10)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)
  ctl
}

test_that("an unscaled readout does not warn that its scaling block is missing", {
  skip_if_not_installed("nonmem2rx")
  # The no-entry WARN exists for a scaling that SHOULD have been emitted and was
  # not. When the DV is an unscaled amount the absence is correct, and warning
  # about it sends the user looking for a bug that is not there.
  #
  # The fixture needs `S<n>` for a compartment OTHER than the observed one, or
  # the no-entry branch is never reached and the test proves nothing. That is not
  # hypothetical: this test first used amount_readout.ctl, whose S2 matches the
  # observed compartment, so `scaling_hint[["2"]]` exists, the branch is skipped,
  # and removing the guard changed nothing. Caught by mutation, not by review.
  dir <- tmp_ctl_dir()
  ctl <- .scale_shape_ctl(dir, "other_cmt_scaled.ctl",
                          "  S1 = V", "Y = A(2)*(1 + EPS(1))")
  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  expect_length(grep("declares scaling for compartment", result$warnings), 0L)
})

test_that("the unscaled-readout INFO does not fire when a y readout is emitted", {
  skip_if_not_installed("nonmem2rx")
  # The INFO explains why no `obs_scale` was emitted. When the $ERROR block
  # produces a readout, `obs_scale` is discarded in favour of the `y` block
  # regardless of what the scaling decision was -- so the message would be
  # announcing a suppression that drove nothing, and its wording ("No [scaling]
  # block was emitted") is simply false: the file does get one.
  #
  # qss_tmdd.mod is the case that showed it. Its $ERROR divides by the PARAMETER
  # VC rather than by `scale1`, so the DV reads as unscaled, while the emitted
  # file carries `y = if (FLAG == 2) c_RTOT else CENT/VC`.
  ctl <- system.file("testmodels/nonmem/qss_tmdd.mod", package = "ferxtranslate")
  skip_if(ctl == "", "fixture not installed")
  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  body <- grep("^#", strsplit(result$ferx_text, "\n")[[1]],
               value = TRUE, invert = TRUE)
  expect_true(any(grepl("^  y( |[[])", body)))
  expect_length(grep("reads the compartment amount", result$warnings), 0L)
})

test_that("only a scale<n> symbol counts as the source having scaled the DV", {
  skip_if_not_installed("nonmem2rx")
  # nonmem2rx spells the $PK `S<n>` assignment `scale<n>`, and the predicate
  # matches that exactly. A user parameter whose name merely CONTAINS "scale"
  # must not be mistaken for it -- here `SCALEF`, an ordinary multiplier, in a
  # model whose $ERROR reads a bare A(2). A substring match would read this as
  # already scaled and emit obs_scale, dividing by V on top of a prediction
  # NONMEM never divided.
  dir <- tmp_ctl_dir()
  ctl <- .scale_shape_ctl(dir, "scalef.ctl",
                          "  SCALEF = 2.0\n  S2 = V",
                          "Y = A(2)*SCALEF*(1 + EPS(1))")
  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  body <- grep("^#", strsplit(result$ferx_text, "\n")[[1]],
               value = TRUE, invert = TRUE)
  expect_false(any(grepl("obs_scale", body, fixed = TRUE)))
  expect_length(grep("reads the compartment amount", result$warnings), 1L)
})

test_that("a $MODEL that cannot be reconciled by name reports the numbering as a gap", {
  skip_if_not_installed("nonmem2rx")
  # The decline path. A $MODEL compartment declared AFTER every DADT-bearing one
  # makes nonmem2rx drop it and lose the remaining names (issue #26), so the COMP
  # list and the state list cannot be matched up: 3 compartments against 2
  # placeholder-named states. No permutation is available, so d/dt order stands
  # -- and ferx will then read a CMT column against a numbering that is not the
  # source's, which is silent in the numbers. That has to be said out loud.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "trailing_gap.ctl")
  writeLines(c(
    "$PROBLEM a trailing $MODEL compartment with no DADT",
    "$INPUT ID TIME DV AMT EVID MDV CMT",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(DEPOT, DEFDOSE)",
    "  COMP=(CENTRAL, DEFOBS)",
    "  COMP=(DUMMY)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    "  S2 = V",
    "$DES",
    "  DADT(1) = -KA*A(1)",
    "  DADT(2) =  KA*A(1) - (CL/V)*A(2)",
    "$ERROR",
    "  Y = A(2)*(1 + EPS(1))",
    "$THETA (0,5) (0,50) (0,1)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  num <- grep("^ERROR .*numbers the compartments", result$warnings, value = TRUE)
  expect_length(num, 1L)
  expect_match(num, "DEPOT, CENTRAL, DUMMY", fixed = TRUE)
  expect_match(num, "including the dose", fixed = TRUE)
  # An action list entry, not just a console message.
  expect_length(grep("numbering could not be reconciled", result$unsupported), 1L)
  # Declined, so d/dt order is left exactly as it was -- no half-permutation.
  expect_length(grep("put in \\$MODEL compartment order", result$warnings), 0L)
})

test_that("an unnamed scaling_hint does not abort the translation", {
  skip_if_not_installed("nonmem2rx")
  # names() on an unnamed list is NULL, so as.integer() gives integer(0) and
  # `!is.na(integer(0))` is logical(0) -- which aborts the `if` rather than
  # declining. rxui_to_ir() is exported and scaling_hint is its argument, so
  # this is a caller-reachable crash, and the same trap the ui$central check
  # below it already records.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "unnamed_hint.ctl")
  writeLines(c(
    "$PROBLEM no compartment named anywhere",
    "$INPUT ID TIME DV AMT EVID MDV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(DEPOT)",
    "  COMP=(CENTRAL)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    "$DES",
    "  DADT(1) = -KA*A(1)",
    "  DADT(2) =  KA*A(1) - (CL/V)*A(2)",
    "$ERROR",
    "  IPRED = F",
    "  Y     = IPRED*(1 + EPS(1))",
    "$THETA (0,5) (0,50) (0,1)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)
  ui <- suppressMessages(nonmem2rx::nonmem2rx(ctl))

  expect_no_error(
    ir <- suppressWarnings(rxui_to_ir(ui, source_format = "nonmem",
                                      scaling_hint = list("V"))))
  # Declined rather than picked: an unnamed entry names no compartment.
  expect_length(grep("taken from \\$PK's S", ir$warnings), 0L)
  expect_length(grep("obs_cmt guessed", ir$unsupported), 1L)
})

test_that("a $MODEL with no DEFOBS is not reported as contradicting the DV expression", {
  skip_if_not_installed("nonmem2rx")
  # .extract_nm_defobs() now returns the COMP list even when no compartment is
  # marked DEFOBS, with name = NA. .same_cmt_name() reports NA as a
  # disagreement, so without an is.na() guard every model that resolves at tier
  # 1 and declares no DEFOBS would be accused of contradicting itself.
  dir <- tmp_ctl_dir()
  ctl <- file.path(dir, "explicit_no_defobs.ctl")
  writeLines(c(
    "$PROBLEM DV names A(2); $MODEL marks nothing",
    "$INPUT ID TIME DV AMT EVID MDV",
    "$DATA d.csv IGNORE=@",
    "$SUBROUTINE ADVAN6 TOL=6",
    "$MODEL",
    "  COMP=(DEPOT)",
    "  COMP=(CENTRAL)",
    "$PK",
    "  CL = THETA(1)*EXP(ETA(1))",
    "  V  = THETA(2)",
    "  KA = THETA(3)",
    "$DES",
    "  DADT(1) = -KA*A(1)",
    "  DADT(2) =  KA*A(1) - (CL/V)*A(2)",
    "$ERROR",
    "  Y = A(2)*(1 + EPS(1))",
    "$THETA (0,5) (0,50) (0,1)",
    "$OMEGA 0.09",
    "$SIGMA 0.04",
    "$EST METHOD=1"), ctl)

  result <- suppressWarnings(nm_to_ferx(ctl, validate = FALSE))
  expect_match(result$ferx_text, "obs_cmt=CENTRAL", fixed = TRUE)
  expect_length(grep("declares 'NA' as DEFOBS", result$warnings), 0L)
  expect_length(grep("but \\$MODEL declares", result$warnings), 0L)
})
