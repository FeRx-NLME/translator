# ferxtranslate 0.1.0.9000

## Breaking changes

* Thetas whose name matches an individual parameter are now renamed
  (`theta CL` becomes `theta TVCL`, with `CL = TVCL * exp(ETA_CL)` unchanged).
  In ferx a theta silently shadows an identically named individual parameter in
  every block where thetas are in scope, so the individual definition was
  written and never read. Code that indexes `fit$theta` by name may need
  updating; the emitted output of `pk_1cmt_oral.mod` and
  `pk_1cmt_oral_ampsim.ctl` changes, and `pk_1cmt_oral.mod` now fits with
  inter-individual variability on clearance for the first time
  (`inst/testdata/ode_1cpt_oral_concordance.csv` was regenerated accordingly).

* `to_ferx()` and its wrappers now validate the emitted `.ferx` with
  `ferx::ferx_model_validate()` and, by default, abort when the engine rejects
  it. Previously a translation that could not be parsed was returned without
  comment. Pass `strict = FALSE` to warn instead, or `validate = FALSE` to skip
  the check.

## New features

* The emitted `.ferx` records ODE state renames as `# renamed: state <from> ->
  <to>` provenance comments, and `ferx_ir` gains a `$state_renames` field. The
  `.ferx` is the artefact that gets shared, and a reader holding only that file
  could not otherwise map a sanitised state such as `c_RTOT` back to the
  `$MODEL` compartment or `A(n)` index it came from. These are provenance, not
  diagnostics: they do not count towards `# Warnings:`, and renamed thetas are
  deliberately excluded because `TVCL` is self-evidently derived from `CL`.

* `to_ferx(validate = TRUE, strict = TRUE)`. When the NONMEM `$DATA` file can be
  resolved relative to the control stream it is passed to the validator, which
  is what allows covariate references to be checked at all -- without a dataset
  every unknown identifier is read as a covariate and reported as valid.

* `ferx_translate_result` gains a `$validation` field holding `ok`,
  `diagnostics`, the `data_file` used, and the notes on how validation ran.

## Bug fixes

* A theta referenced from an ODE is now carried in by an individual parameter.
  ferx resolves a `d/dt` right-hand side against states, individual parameters and
  covariates only -- a theta is not in scope there -- so a `$DES` line naming a
  theta produced an `[odes]` block ferx could not resolve (issue #6, defect 2).
  Measured against ferx 0.2.0 and 0.3.0, a theta IS readable from
  `[individual_parameters]`, `[scaling] y` and `obs_scale`, and is NOT readable
  from a `d/dt` right-hand side, an ODE-block intermediate, an `init()` expression
  or a pk macro argument; only the second group gets a carrier, so a theta read
  from `[individual_parameters]` is left exactly as it was.

  Where the source already supplies the assignment nothing is added: `KTP =
  THETA(3)` arrives as the alias `KTP <- KTP`, which the alias filter drops as a
  self-assignment but which theta de-shadowing turns into the surviving `KTP <-
  TVKTP` -- the same carrier, spelled by the model. `inst/testmodels/nonmem/ode_theta_ref.ctl`
  covers the shape that has nothing to convert.

  A carrier is named after the source (`KTP = TVKTP`) where that name is free,
  which is the form ferx's own examples use and keeps the emitted `[odes]` diffable
  against the source `$DES`. Where it is taken the name is derived from the theta's
  *emitted* name instead: `TVCL_ODE = TVCL`. Not `CL_ODE`, which reads as "the CL
  used in the ODE" -- that is the individual value, and distinguishing it from the
  theta is the entire point. Not `CL_1` either: `_1` already means "state
  disambiguated" (`central_1`) and is `.free_theta_name()`'s last resort, so a
  third meaning on that spelling leaves a reader unable to tell a renamed
  compartment from a renamed theta from a carrier, and `CL_1`/`V_1` are plausible
  model variables in their own right. A numbered form is still possible if both
  spellings are taken, and now warns: the number is positional, so it moves when
  another carrier is added before it.

* A carrier reuses an existing individual parameter only when that parameter is a
  *pure alias* of the theta. Matching on the name alone gets the arithmetic wrong:
  `frac <- central/cl` in a model that later writes `cl <- cl*exp(eta.cl)` reads
  the theta, so pointing that reference at the individual parameter `CL` silently
  substitutes the IIV-applied value. Both forms parse and both fit. Such a
  reference now gets its own carrier (`CL_1 = TVCL`).

* The `[odes]` scope is applied to the emitted right-hand side rather than during
  the walk, because that is where every path converges. An ODE-block intermediate
  that touches a state is inlined instead of emitted, and its text was normalised
  for a different context: `ki <- KTP*CENT` reached `[odes]` as `TVKTP * CENT`, a
  bare theta, with the carrier defined and unreferenced beside it. Resolving the
  `d/dt` line alone also made the result depend on statement order -- the alias
  written above the `d/dt` bound in time, the identical line below it did not.

* A theta is recognised as ODE-referenced by its *current emitted* name, or by its
  source name when nothing else declares that symbol -- not by source name alone.
  In any de-shadowed model the source name is the individual parameter's name, so
  `d/dt(ABS) = -KA * ABS` beside `theta TVKA` and `KA = TVKA * exp(ETA_KA)` reads
  the parameter, and matching on `KA` invented a second carrier for a reference
  that was already correct. The test is recomputed each round rather than
  accumulated, for the same reason: a theta matches on its source name only until
  de-shadowing renames it.

* nonmem2rx's `rxmissingvars` placeholders are no longer emitted. When a NONMEM
  model names a theta by its `$THETA` label (`FLUX = KTP*A(1)` for `(0,0.2) ;
  KTP`), nonmem2rx does not bind the symbol: it leaves a free `KTP` and records
  the theta in `rxmissingvars1 <- t.KTP`. That produced the meaningless individual
  parameter `RXMISSINGVARS1 = TVKTP` while the reference itself stayed dangling.
  The placeholder is dropped and the reference bound by the carrier.

  Unlike the covariate case noted above, `[odes]` does not need a dataset for the
  engine to object: it reports `[odes]: RHS references undefined name(s)` either
  way, because thetas and etas are out of scope there and a bare name can only be
  a state, a parameter or a covariate.

* Every name the emitted `[odes]` block references is now checked against what the
  model declares, and an unresolved one is an `ERROR` with an `$unsupported` entry.
  ferx does report this itself, but only where the engine runs -- and the
  identifier-legality check beside this one tests the *grammar*, not whether
  anything declares the name, so `KTP` passed it. Two defects reached the corpus
  through that gap (issue #6 defects 2 and 4).

  The check is possible because the set of names a ferx ODE right-hand side may
  reference is *closed*: declared states, individual parameters, ODE-block
  intermediates, and `TIME`/`TAFD`/`TAD`/`MACHEPS`. Not thetas, etas or sigmas --
  and not covariates either. ferx's own message gives the remedy, which this one
  repeats: pre-compute the covariate-dependent term in `[individual_parameters]`
  and reference it from the ODE by that name.

  No covariate list is consulted, and that is deliberate rather than an omission.
  `.covariate_names()` *defines* a covariate as a symbol nothing else binds, and
  rxode2's `ui$allCovs` does much the same, so both classify a name the translator
  failed to bind as a legitimate covariate -- measured, both call the unbound `CF`
  in the QSS TMDD model and the unbound `KTP` in `ode_theta_ref.ctl` covariates.
  An earlier version of this check consulted them, passed the entire test suite,
  and could not fire on either defect it was written for.

  Two reports, because the remedies differ: a theta, eta or sigma needs a carrier;
  anything else resolves to nothing at all. The first makes an eta referenced from
  an ODE reportable for the first time -- ferx rejects it as `E_PARSE` and it was
  previously emitted without comment.

  Scoped to `[odes]` on purpose. In `[individual_parameters]` thetas *are* in
  scope, so the set of legitimate names is far larger and a leftover carries much
  less information.

  Measured against ferx 0.3.0 (ferx-r tag `v0.3.0`), which is the version shipped
  to users. 0.2.0 was not re-measured; no bundled model references a covariate from
  `[odes]`, so nothing that translated before is affected either way. Note that the
  `engine` CI job still pins `ferx-r@731adc9` = 0.2.0, so it validates against a
  build no user runs -- tracked separately from this change.

* A theta renamed only because it is ODE-referenced is no longer reported as
  shadowing an individual parameter. At the point of the rename nothing of that
  name is one -- the trigger predicts a carrier that is created afterwards -- so
  the message described a collision the model never had.

* A theta whose name is also an ODE state name is not treated as an ODE reference:
  in `[odes]` the state wins. Without the exclusion the emitted ODE still came out
  right, but the theta was renamed for a reference that was never to it.


* Emitted identifiers are now sanitised to the ferx grammar
  (`[A-Za-z_][A-Za-z0-9_]*`). ODE state names never went through normalisation,
  so a compartment that `nonmem2rx` had renamed to `c.RTOT` -- which it does
  whenever a `$MODEL` compartment collides with a variable name -- made the whole
  file unparseable (issue #6, defect 1). The rename is applied to the
  declaration and to every reference at once (`obs_cmt=`, `states=[...]`, the
  `d/dt` target, and the state inlined into other ODE right-hand sides), because
  renaming the declaration alone leaves the references pointing at nothing.
  Only names that must change do change: an already-legal state keeps its name
  and its case.

* Covariate references are deliberately left alone, and an illegal one is now
  reported as an `ERROR` rather than renamed. ferx matches covariates to data
  columns by exact name, case included, so sanitising one would move the failure
  from `E_PARSE` to `E_MISSING_COVARIATE` at fit time instead of fixing it. The
  dataset column has to be renamed. A regression test pins the existing
  case-preserving behaviour, which nothing previously covered.

* The ODE state rename is resolved through a dedicated state map rather than the
  general name map. The name map holds every parameter alias and grows as the
  model block is walked, which made the emitted state name something nobody
  decided: a state was renamed to a *parameter's* name whenever its raw spelling
  happened to be an `iniDf` key (`d/dt(CENT)` became `d/dt(VC)`), and an
  unrelated assignment written above rather than below a `d/dt` line changed
  whether the state was renamed at all. `obs_cmt` is resolved the same way; it
  previously read a different map from the `d/dt` target, so the two could
  disagree and emit an `obs_cmt` naming no declared state.

* The state sanitiser reserves every assignment target, covariate and reserved
  ODE name, not just the random-effect names. A sanitised state landing on an
  assignment target was absorbed into the auxiliary-variable set, dropped from
  `[individual_parameters]`, and its references re-resolved to the state --
  turning `d/dt(A.B) = -A_B * A.B` into `d/dt(A_B) = -A_B * A_B`, the amount
  squared, with the rate constant and its IIV gone. ferx validates that output
  clean, so nothing downstream would have caught it.

* The identifier-legality check runs on every symbol emitted verbatim, not on
  the covariate set. Covariates are classified by *normalised* name, so an
  illegal raw symbol whose normalised form matched a known name was filtered out
  before the check saw it, while still being emitted verbatim.

* Smaller fixes from the same review: `ui$central` is length-guarded (a
  zero-length or multi-element value aborted the translation); `.ferx_ident()`
  maps `NA` to a legal identifier rather than passing it through; `.ddt_state()`
  takes the first element, so a non-symbol `d/dt()` argument cannot invent a
  second state; and the two hand-inlined copies of the old normalisation rule
  (in `[scaling]` matching and in `.extract_sigmas()`) now call `.norm()`.

* `inst/testmodels/nonmem/dotted_state.ctl` joins the corpus, and the sweep in
  `test-integration.R` now tokenises the *emitted text* rather than listing IR
  declaration fields. An illegal name parses where it is declared and fails only
  where it is referenced, so a declaration-only check inspected the half that
  works -- it passed cleanly on a file the engine rejects.

* The concordance data generator (`data-raw/generate_concordance_data.R`)
  produced empty datasets against `ferx` 0.2.0: its templates wrote `DV = "."`,
  which the engine reads as "no observation" and skips. It now uses a numeric
  placeholder, restores `.` on dose rows when writing, and asserts that one
  value came back per observation row. Reported upstream as
  FeRx-NLME/ferx-r#283 and FeRx-NLME/ferx-core#957.

* The translation gap report now validates every bundled model and fails on
  engine errors, instead of only listing features the translator already knew it
  could not express.

## Internal

* `ferx` and `amp.sim` are declared in `Suggests`; both were already used by the
  test suite without being declared.

* The CI engine pin moves from `ferx-r@54f25d4` (0.1.5) to `ferx-r@731adc9`
  (0.2.0), the build the concordance suite is now verified against. No reference
  *value* changed: the theta assertions were renamed to follow the `TV` prefix,
  two omega assertions were added, and `ode_1cpt_oral_concordance.csv` was
  regenerated because its source model changed. The other three datasets
  regenerate byte-identical against this pin. The 2-cpt omega references still
  carry their original 0.1.x provenance and were re-verified, not re-fitted.

* The `engine` CI job now runs the full test suite rather than
  `filter = "concordance"`. Four of the validation tests require `ferx`, so a
  filtered run left the output-validation feature executed by neither job.
