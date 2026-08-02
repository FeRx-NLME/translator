# ferxtranslate (development version)

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

* `ferx` and `amp.sim` are declared in `Suggests`; both were already used by the
  test suite without being declared.

* The CI engine pin moves from `ferx-r@54f25d4` (0.1.5) to `ferx-r@731adc9`
  (0.2.0), matching the re-baselined concordance references and datasets.
