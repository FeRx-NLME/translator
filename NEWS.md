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
