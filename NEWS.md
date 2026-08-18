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

* `inst/testmodels/nonmem/dotted_state.ctl` joins the corpus, and the sweep in
  `test-integration.R` now asserts that every emitted identifier is legal.
  Without such a model in the corpus that assertion would hold vacuously.

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
