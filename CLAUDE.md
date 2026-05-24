# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

`ferxtranslate` is an R package that translates pharmacometric models written
in NONMEM, nlmixr2, or Monolix into ferx `.ferx` format, so users can move to
ferx without rewriting their models by hand. mrgsolve support is planned for
v0.2.

See `plans/v0.1-implementation.md` for the full design, architecture, and
day-by-day build plan.

## Sibling repositories

| Repo | Local path | Purpose |
|---|---|---|
| ferx-r | `../ferx-r` | R package for ferx; source of `.ferx` format spec and example models |
| ferx-core | `../ferx-core` | Rust engine; defines what `.ferx` files can express |

When in doubt about what ferx supports, read `../ferx-r/inst/examples/models/`
and `../ferx-core/src/parser/model_parser.rs`. The ferx format is the target —
never invent syntax that is not already accepted by the parser.

## Build and check commands

```bash
# Install the package locally
R CMD INSTALL .

# Run all tests
Rscript -e 'devtools::test()'

# Rebuild documentation (do this before every PR)
Rscript -e 'roxygen2::roxygenize()'

# Full CRAN check
R CMD check --as-cran .

# Quick check during development
Rscript -e 'devtools::check(cran = FALSE)'
```

## Known gotchas

**nonmem2rx drops `S2=V`** — NONMEM `$PK` scaling assignments (`S1`, `S2`, etc.) are
silently omitted from `ui$lstExpr`. Without them, an ODE model predicts amounts but
data are concentrations, so IPRED >> DV and estimation diverges silently.
Always parse the raw `.ctl`/`.mod` file for scaling via `.extract_nm_scaling()` (in
`R/utils.R`); do not rely on the rxode2 UI object for this.

**Fixed-effect PK params are absent from `indiv_params`** — When a PK param has no
ETA (e.g. `V = THETA(3)`), nonmem2rx does not emit it as an assignment in `lstExpr`.
The linCmt pk macro arg lookup then silently drops `v=V`, producing IPRED=0.
The passthrough logic in `rxui_to_ir.R` handles this; do not remove it.

**Snapshot acceptance after output changes** — If a code change affects the `.ferx`
text of any bundled test model, the integration snapshots in
`tests/testthat/_snaps/integration.md` will fail. Run
`testthat::snapshot_review("integration")` to inspect the diff before accepting.
Only accept if the new output is deliberately correct.

**amp.sim package** — `amp.sim` (GitHub: LeidenAdvancedPKPD/amp.sim) is used for
the external NONMEM reference benchmark in `test-concordance.R`. It is a `Suggests`
dependency. Install with `remotes::install_github("LeidenAdvancedPKPD/amp.sim")`.
`NM.theoph.02B.csv` is NOT bundled in amp.sim — the concordance dataset is
pre-simulated and stored in `inst/testdata/ampsim_1cpt_oral_concordance.csv`.
If the amp.sim reference estimates ever change, re-run
`data-raw/generate_concordance_data.R` to regenerate the dataset.

## Four-tier test structure

Every new function or behaviour needs a test. Put the test in the lowest tier
that covers it. Do not write tests at the end — write them as you build.

**Tier 1 — Unit tests** (`tests/testthat/test-*.R`, inline, no file I/O)

Test the smallest unit in isolation. Use inline R objects (hand-built
`ferx_ir` lists, inline `rxode2::rxode2()` model functions). Must run in
milliseconds. These run on every PR and block merge if they fail.

```r
test_that("diagonal omega emits correctly", {
  ir <- new_ferx_ir(omegas = list(list(type="diagonal", name="ETA_CL", value=0.07)))
  expect_match(emit_ferx(ir), "omega ETA_CL ~ 0.07")
})
```

**Tier 2 — Integration tests** (`tests/testthat/test-integration-*.R`)

Call the full pipeline (`to_ferx()` / `nm_to_ferx()` / etc.) but use inline
model definitions or small bundled `.ctl` files from `inst/testmodels/`.
Must complete in under 10 seconds. These also run on every PR.

```r
test_that("1-cpt oral NONMEM model round-trips to ferx", {
  ctl  <- system.file("testmodels/nonmem/1cpt_oral.ctl", package = "ferxtranslate")
  result <- nm_to_ferx(ctl)
  expect_snapshot(result$ferx_text)
  expect_length(result$unsupported, 0)
})
```

**Tier 3 — Reference snapshot tests** (`tests/testthat/test-snapshots-*.R`)

Compare `emit_ferx()` output against committed reference `.ferx` files in
`inst/testmodels/reference/`. These are the ground-truth correctness tests.
Run with `devtools::test()` locally; also run in CI. If a snapshot changes,
review the diff carefully before accepting — it means the translation output
changed for a real model.

Reference snapshots live in `tests/testthat/_snaps/`. Accept updated snapshots
with `testthat::snapshot_accept()` only after manually verifying the new output
is correct.

**Tier 4 — Numerical concordance tests** (`tests/testthat/test-concordance.R`)

Translate a bundled model, fit pre-simulated data with `ferx_fit()`, and assert
that estimated parameters are within tolerance of the known true values. These
are the only tests that can catch silent semantic errors: wrong ODE sign, swapped
parameter, missing scaling, wrong sigma interpretation.

Gated with `skip_if_not_installed("ferx")` and `skip_on_ci()` — run locally only.
Require the `ferx` binary in PATH and take ~2 minutes.

```r
# Run concordance tests locally
Rscript -e 'devtools::test(filter="concordance")'
```

Current test suite and tolerances:

| Test | Model | True params | Tolerance |
|---|---|---|---|
| linCmt 1-cpt oral: TVCL/TVV | `1cpt_oral.ctl` | TVCL=0.134, TVV=8.1 | 15% |
| linCmt 2-cpt IV: 4 thetas | `2cpt_iv.ctl` | CL=5, V1=20, Q=8, V2=60 | 10% |
| linCmt 2-cpt IV: omegas | `2cpt_iv.ctl` | om_CL/V1=0.10, om_Q=0.08 | 10% |
| amp.sim linCmt benchmark | `pk_1cmt_oral_ampsim.ctl` | KA=0.0825, CL=2.676, V=1.588 | 10% |
| ODE 1-cpt oral with S2=V | `pk_1cmt_oral.mod` | KA=0.1, CL=2.0, V=1.0 | 15% |

Datasets live in `inst/testdata/`. Regenerate with `data-raw/generate_concordance_data.R`
if model files or theta initials change. Commit the regenerated CSVs.

**Translation gap report** — a dedicated test translates every model in
`inst/testmodels/nonmem/` and prints a `model -> gap` table of any
`$unsupported` features. It always passes; the table is the signal. Do NOT
add hardcoded skip skeletons for individual gaps — the report is generic and
picks them up automatically. Add a model file to extend coverage; the gap
disappears from the report when the translator or ferx-core gains support.

## Documentation rules

**Run `roxygen2::roxygenize()` and commit the updated `man/` files before
opening or updating any PR.** The `.Rd` files are checked into the repo.

Every exported function must have:
- `@title` (implicit from first sentence)
- `@param` for every argument
- `@return` describing the `ferx_translate_result` fields or other return value
- `@examples` block — use `\dontrun{}` for anything requiring file paths or
  heavy packages

**No non-ASCII characters anywhere in `R/*.R` files.** This causes `R CMD check
--as-cran` failures and broken help pages. Violations to avoid:

| Avoid | Use instead |
|---|---|
| `\uXXXX` escape sequences in comments | literal ASCII |
| em-dash `—`, en-dash `–` | `-` or `:` |
| ellipsis `...` (Unicode) | `...` (three dots) or `etc.` |
| box-drawing characters | `-- Section --` |

Check before every PR:

```r
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  n <- sum(chartr(rawToChar(as.raw(128:255)), strrep("x", 128),
                  readLines(f, warn = FALSE)) !=
           readLines(f, warn = FALSE))
  if (n > 0) message(f, ": ", n, " non-ASCII lines")
}
```

Or from the shell:

```bash
python3 -c "
import glob
for f in glob.glob('R/*.R'):
    if any(b > 127 for b in open(f,'rb').read()):
        print(f)
"
```

## Pull requests

Before opening a PR:

1. Read `.github/PULL_REQUEST_TEMPLATE.md` and fill every section.
2. Run `roxygen2::roxygenize()` and commit `man/`.
3. Run the non-ASCII check above.
4. Run `R CMD check --as-cran .` — zero ERRORs, zero WARNINGs.
5. Run `devtools::test()` — all tiers must pass.

One PR per logical change. Do not bundle unrelated fixes.

If a change here requires a corresponding change to the ferx format (a new
`.ferx` keyword, a new `[section]`, a new pk macro), open a ferx-core PR first
and link it. This package is a consumer of the ferx format — it must never
generate syntax that the current ferx parser rejects.

## Warning system

Every translation warning must be:

- Prefixed with `INFO`, `WARN`, or `ERROR` in the stored string
- Emitted via `cli::cli_warn()` or `cli::cli_inform()` at translation time
  so the user sees it immediately in the console
- Stored in `ir$warnings` for programmatic inspection via `result$warnings`
- Placed as a `# WARNING: ...` comment at the exact location in the `.ferx`
  output where the unsupported feature would have appeared

See `plans/v0.1-implementation.md` Section 9 for the full catalogue of
translatable, lossy, and untranslatable features.

## What ferx does not yet support

Features listed in `plans/v0.1-implementation.md` under "ferx feature roadmap"
are things we detect in source models but cannot emit because ferx-core does
not support them yet. When one of these is encountered:

- Do NOT silently drop it.
- Add an `ERROR`-level warning and a `# WARNING:` comment in the output.
- Add it to `result$unsupported`.

This gives users a clear action list and gives the ferx-core team a concrete
signal of what to prioritise.
