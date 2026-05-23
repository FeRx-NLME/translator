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

## Three-tier test structure

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
