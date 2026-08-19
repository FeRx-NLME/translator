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

* Every ODE state -- renamed or not -- is pinned on the maps the model block is
  parsed with, so neither an `iniDf` key that happens to share a state's name nor
  an assignment walked earlier can rebind one. The general name map holds every
  parameter alias and grows as the walk proceeds, which made the emitted state
  name something nobody decided, on the declaration *and* on every reference:
  `d/dt(CENT) = -K * VC` when a theta keyed `CENT` was labelled `VC` -- the
  compartment amount replaced by a fixed theta, with no diagnostic and a file the
  engine accepts -- and an ODE whose right-hand side changed depending on whether
  `central <- 0` stood above or below the `d/dt` line. `obs_cmt` is resolved the
  same way; it previously read a different map from the `d/dt` target, so the two
  could disagree and emit an `obs_cmt` naming no declared state.

* A state whose source name is also a model parameter's is now resolved by
  scope -- the state inside `[odes]`, where thetas and etas are out of scope in
  ferx too, and the parameter everywhere else -- and reported as an `ERROR`,
  because the source cannot say which reading it meant. It was previously
  refused outright, which left the state sharing a name with the parameter: the
  assignment referencing it was absorbed as an ODE intermediate, dropped from
  `[individual_parameters]` (leaving the block empty), and self-inlined to the
  depth cap, emitting the same `exp(ETA_X)` factor fifteen times.

* A symbol referenced *before* its defining assignment is now emitted under the
  same name as one referenced after it. The alias was installed only once the
  walk had passed the assignment, so `d/dt(central) = -cl*central*f.rac` written
  above `f.rac <- 0.5` emitted the illegal `f.rac` verbatim beside a declared
  `F_RAC` -- two spellings of one variable, an unparseable file, and no warning,
  since the legality check treats every assignment target as legal by
  construction regardless of position. Dotted local names are idiomatic nlmixr2.

* A theta named after a ferx solver builtin (`TIME`, `T`, `TAFD`, `TAD`,
  `MACHEPS`) is renamed, and an individual parameter with such a name is
  reported. ferx-core checks its reserved list against states, individual
  parameters and ODE intermediates alike; only states were covered. The theta
  half was silent rather than loud: ferx resolves the bare name to the value the
  solver injects, so `KA = TIME * exp(ETA1)` made `KA` read the integrator clock
  and left the estimated theta unreferenced, reported only as a
  `W_UNUSED_PARAM`.

* A NONMEM sigma whose source name needs sanitising still reaches
  `[error_model]`. `nonmem2rx` keeps sigma out of `iniDf`, so nothing bound its
  source spelling to its emitted name and the two agreed only while both were
  plain uppercase; once the declaration was sanitised, the error assignment
  matched no sigma and the block was omitted entirely -- a model with no residual
  error, and an `ERROR` misreporting the sigma as an illegal covariate.

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

* The CI engine pin moves to `ferx-r@7889719` (0.3.0, tag `v0.3.0`) -- the version
  users install. It had stayed on `ferx-r@731adc9` (0.2.0) after 0.3.0 shipped,
  which meant the `engine` job proved only that a build nobody runs accepted the
  emitted `.ferx`. That is the one job standing between a silent engine regression
  and a merge, so pointing it at a stale engine defeats its purpose.

  No re-baselining was needed: the full suite passes against 0.3.0 with the
  concordance references in `test-concordance.R` and the bundled datasets in
  `inst/testdata/` unchanged. `CLAUDE.md` says both are tied to the pin and must be
  re-baselined alongside it; that stays true as guidance -- this bump happening to
  need none is not evidence the next one won't.

  The pin is a SHA rather than the tag on purpose: a tag can be moved.

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
