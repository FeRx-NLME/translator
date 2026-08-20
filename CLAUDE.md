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

## Working on a shared machine

Several sessions may be working this repo at once, and they share one machine,
one system-wide R library and one CI runner pool. Every item below is something
that actually went wrong here, not a precaution.

**A broken package install leaves a half-populated directory and a live
`00LOCK-*`. Stop touching it. Do not clean up.** This is the counter-intuitive
one: `R CMD INSTALL` moves the previous version into `00LOCK-<pkg>/` and restores
it if the build fails, so the lock IS the backup. Tidying away a stale-looking
lock, or moving the old package back yourself, destroys the thing R is about to
restore. An empty `library/ferx` next to a `00LOCK-ferx-r` is a recoverable state
and usually an install still running - check `pgrep -fl cargo` and the lock's
mtime before concluding it is abandoned.

**Never kill by a pattern that is not scoped to your own build.**
`pkill -f "cargo build --release"` matches every session's build, not yours.
`pkill -f "rustc --crate-name"` matches every Rust compile on the machine. Run
`pgrep -fl <pattern>` first and read the list; kill a PID you captured yourself,
or a pattern containing your own build directory. A ferx build costs 15-30
minutes and blocks the engine test tier, so a stray kill is expensive for whoever
started it - and you will not be told.

**Announce before installing into the system-wide R library.** `.libPaths()[1]`
is `/Library/Frameworks/R.framework/...`, shared by every session. Two concurrent
installs of `ferx` or a shared dependency clobber each other, and the symptom
(the previous paragraph) does not look like a collision.

**Pushing to a branch cancels that branch's in-flight CI run.** `check.yml` has a
`concurrency:` block keyed on `github.ref`, so a rapid series of pushes leaves only
the last commit tested. That is deliberate - four concurrent runs once piled up on
one branch here in fifteen minutes, three of them testing commits already moved
past - but it means that if you want a specific commit's engine result before
merging, you have to let that run finish rather than pushing on top of it.

Runs on `main` are never cancelled:

```yaml
cancel-in-progress: ${{ github.ref != 'refs/heads/main' }}
```

With no required status checks, a `main` commit's engine run is the only thing
verifying what actually landed, so those must always complete. A superseded branch
run is merely wasted; a cancelled `main` run leaves a commit permanently
unverified on the one branch where that matters.

What that rule protects is the RESULT, not the run. If the `engine` job on a
`main` run has already completed successfully and only `fast` is wedged, the
verification is banked and cancelling the run throws nothing away -- do it, and
free the runner. `fast` duplicates a check you can run locally; `engine` does not,
which is the whole reason main runs are exempt from cancellation. Read
`cancel-in-progress` as "never lose a main commit's engine result", not as "never
touch a main run": the second reading makes you hesitate at exactly the moment
cancelling is correct. This is not hypothetical -- a wedged `fast` sat on a main
run for 111 minutes whose `engine` had passed in five.

**A hung CI job looks identical to a slow one in `gh pr checks`.** It prints
`pending` and nothing else, which is worth an hour if you let it be. The
step-level view is the one that answers the question:

```bash
gh api repos/FeRx-NLME/translator/actions/runs/<run-id>/jobs \
  --jq '.jobs[] | "== \(.name) \(.status)", (.steps[] | "   \(.name): \(.status) \(.started_at) -> \(.completed_at // "-")")'
```

The two jobs are each other's control: they run the same `setup-r-dependencies@v2`
against different package sets, so a step that takes 63 seconds in one and 65
minutes in the other is hung, not slow. Observed on both jobs, on the same commit,
with the same YAML -- so it is not specific to either dependency set, and a
dependency that succeeds in eight minutes on one run and hangs for fifty-three on
the next is not the cause. The remedy is a re-run, not a workflow change: cancel,
`gh run rerun <run-id>`, and the same step completes in about a minute. Suspect a
transient in the action before you suspect the diff.

## Known gotchas

**A theta that shares a name with an individual parameter shadows it** - ferx
resolves an identifier as theta first, then eta, then individual parameter, in
every block where thetas are in scope (`[individual_parameters]`, `[scaling]`,
`[initial_conditions]`). `CL = CL * exp(ETA_CL)` followed by `K20 = CL/V` reads
the *theta* in the second line, so the first is dead and the IIV silently
vanishes. ferx emits no diagnostic. `.deshadow_theta_names()` in `rxui_to_ir.R`
renames the theta (`CL` -> `TVCL`) and is the single owner of theta naming; it
also uniquifies duplicate `$THETA` labels, which ferx would otherwise resolve to
the first. Do not rename thetas anywhere else, and do not compare theta names
against individual-parameter names by prediction - the parser is run twice for
exactly this reason.

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

**A schedule-only workflow is not registered until it has run once** - GitHub does
not index a workflow whose only triggers are `schedule` and `workflow_dispatch`
until an actual run exists for it. Until then it is absent from `gh workflow list`
and from the Actions API, and `gh workflow run` answers
`404: not found on the default branch` - even though the file is on the default
branch, is valid YAML, and Actions is enabled. Merging one and walking away leaves
something that looks installed and does nothing.

It is a chicken-and-egg, not a defect: it cannot be dispatched because it is not
registered, and it is not registered because it has never run. To break it, push
the workflow to a throwaway branch with a temporary `push:` trigger; that first run
registers the file path permanently, after which `workflow_dispatch` works on the
default branch and the temporary trigger can be dropped. `engine-pin-drift.yml` was
registered exactly this way.

Worth being deliberate about, because the workflows most likely to be
schedule-only are watchdogs - the ones nobody watches, whose whole purpose is to
notice something on your behalf. A silently unregistered watchdog is worse than
none: it reads as coverage.

**After a merge, parse-check the test files before trusting the suite** - a
conflict boundary can cut through a `test_that` block and leave it without its
closing `})` while `R/` auto-merges cleanly and is genuinely correct. This
happened on the #15 merge, to `test-rxui_to_ir.R`.

Measured on this package, `testthat::test_local()` then halts:

```
ir: Error in parse(...) : test-ir.R:181:0: unexpected end of input
Execution halted
```

Loud, but it *aborts the run*, so what reaches a summary line is a short or
missing result rather than a failure. A clean auto-merge plus a suite that says
little reads exactly like success.

The pre-flight is cheap and names the file faster than a suite run does:

```bash
for f in tests/testthat/*.R R/*.R; do Rscript -e "invisible(parse('$f'))"; done
```

Then compare the test count against the previous run. `FAIL 0` on a suite that
quietly shed several hundred tests has not gotten better, and the count is the
only thing that notices. Same family as the `main` run-cancellation rule above:
the failure mode is an absence of signal rather than a signal.

## Four-tier test structure

Every new function or behaviour needs a test. Put the test in the lowest tier
that covers it. Do not write tests at the end — write them as you build.

**A fixture picked to show a behaviour must be able to show that behaviour being
wrong.** Before adding a test, break the code it guards and confirm the test
fails. If it still passes, the fixture cannot distinguish the two cases and the
test proves nothing - it will sit in the suite reading as coverage.

Three have shipped that way, all three through review:

- `F1 = 1` cannot show bioavailability being dropped. Measured: `F1 = 1` and no
  `F` at all give predictions identical to every printed digit, so a concordance
  test built on it passes whether the feature works or not.
- Bare `init(X) = TIME` cannot show drop-versus-substitute. It is the one
  spelling where dropping the statement and substituting `TIME := 0` give the
  same answer; `init(X) = K + TIME` is the case that separates them.
- An omega guard starting at its own simulation truth cannot show an eta failing
  to reach the ODE - "recovered the truth" and "never moved" are the same
  observation. See the Tier 4 notes below.

None of these were caught by review. All three would have been caught by one
question: what would this test do if the code under it were broken? That is not
a maxim, it is a command - run it.

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

Gated with `skip_if_not_installed("ferx")`. They require the `ferx` package
(the Rust engine) and take ~2 minutes. They run:

- **Locally** whenever `ferx` is installed - `Rscript -e 'devtools::test(filter="concordance")'`.
- **In CI** only in the dedicated `engine` job (`.github/workflows/check.yml`),
  which installs a pinned `ferx` and sets `FERXTRANSLATE_ENGINE_TESTS=true`. The
  fast PR job skips them. This is what makes a green check actually mean the
  engine accepted and fit the emitted `.ferx` - do not re-add a blanket
  `skip_on_ci()`, or the only tier that exercises the engine stops running in CI.

The `engine` job pins `ferx` (currently `ferx-r@9c97c13`). Pin the SHA, not the
tag: a tag can be moved, and this job is the only thing between a silent engine
regression and a merge.

Keep the pin on the version users actually install, and note that "the version
users install" means MAIN, not the latest tag - both READMEs say
`pak::pak("FeRx-NLME/ferx-r")`, which takes the default branch. The rule has now
been broken in both directions. It sat on `ferx-r@731adc9` (0.2.0) after 0.3.0
shipped, so every green engine job proved only that a build nobody runs accepted
the output; and pinning the `v0.3.0` tag today would do the same thing again,
because ferx-r main has moved past it.

Do not identify the pinned engine by its version string. `9c97c13` still reports
0.3.0 in DESCRIPTION even though ferx-core's parser changed underneath it
(`E_DOSE_ATTR_DOUBLE_USE`, ferx-core#993/#1003), so "ferx 0.3.0" names two
different engines. Cite the ferx-r SHA.

The pin is TRANSITIVE and reading only ferx-r's `Cargo.toml` will mislead you:
it takes ferx-core as `branch = "main"`, which looks unpinned, while the
committed `src/rust/Cargo.lock` carries the exact ferx-core rev. So there is no
single ferx-core SHA this repo can pin, and picking up a ferx-core change is a
two-step cross-repo bump - ferx-r's lock first (via its
`tools/update-ferx-core-lock.sh`; never a bare `cargo update`, which its local
`[patch]` silently unpins), then this pin to that ferx-r commit.

The reference omegas in `test-concordance.R` AND the bundled datasets in
`inst/testdata/` (regenerated by `data-raw/generate_concordance_data.R` against
that engine) are tied to that build - re-baseline both in the same commit as a
bump IF they move. They did not move for 0.2.0 -> 0.3.0: the full suite passed
against 0.3.0 with no change to either. Check rather than assume, in both
directions - a bump that needs no re-baseline is not evidence the next one won't.

```r
# Run concordance tests locally (ferx installed)
Rscript -e 'devtools::test(filter="concordance")'
```

Current test suite and tolerances:

| Test | Model | Reference | Tolerance |
|---|---|---|---|
| linCmt 1-cpt oral: TVCL/TVV | `1cpt_oral.ctl` | truth: TVCL=0.134, TVV=8.1 | 15% |
| linCmt 2-cpt IV: 4 thetas | `2cpt_iv.ctl` | truth: CL=5, V1=20, Q=8, V2=60 | 10% |
| linCmt 2-cpt IV: omegas | `2cpt_iv.ctl` | reference fit (ML), not truth | 10% |
| amp.sim linCmt benchmark | `pk_1cmt_oral_ampsim.ctl` | NONMEM: KA=0.0825, CL=2.676, V=1.588 | 10% |
| ODE 1-cpt oral with S2=V | `pk_1cmt_oral.mod` | truth: TVKA=0.1, TVCL=2.0, V=1.0 | 15% |
| ODE 1-cpt oral omegas | `pk_1cmt_oral.mod` | truth: ETA_KA=0.01, ETA_CL=0.02 | 35% |

Structural thetas are well identified, so they are asserted against the nominal
simulation truth. Omega variances from 50 subjects carry ~20% sampling SE, so
the 2-cpt IV omega test asserts against the **reference ML fit** (a deterministic
regression check on eta wiring), not the nominal truth - see the comment in
`test-concordance.R`. Re-baseline those reference omegas if the engine or the
dataset changes.

The ODE omega test is the exception and asserts nominal truth at a wider 35%,
because it is a regression guard rather than a precision check: it exists to
prove `ETA_CL` still reaches the ODE through the derived `K20 = CL/V`. Note
that a model's `$OMEGA` initials are often the simulation truth, which makes
"recovered the truth" and "never moved off its start" indistinguishable - that
assertion was a tautology until the test was changed to start `ETA_CL` away
from the truth and additionally assert that it moved. This is one instance of
the discriminating-fixture rule at the top of this section, which applies to
every tier, not only to omega guards - scoping it to omega guards is why it
went on to bite twice more.

Theta names carry a `TV` prefix wherever the source names a theta after the
parameter it defines (see the de-shadowing note under Known gotchas), so
concordance references index `fit$theta[["TVCL"]]`, not `fit$theta[["CL"]]`.

Datasets live in `inst/testdata/`. Regenerate with `data-raw/generate_concordance_data.R`
if model files or theta initials change. Commit the regenerated CSVs.

**Corpus sweeps** — two tests translate every model in `inst/testmodels/nonmem/`.
Do NOT add hardcoded skip skeletons for individual gaps; both are generic and
pick up a new model file automatically.

- `test-integration.R` (engine-free, runs on every PR): asserts nothing crashes
  and that no theta shadows an individual parameter, and prints a
  `model -> gap` table of `$unsupported`. Gaps are reported, not failed — the
  table is the signal, and a gap disappears when the translator or ferx-core
  gains support. Crashes and shadowing DO fail.
- `test-concordance.R` (needs `ferx`): validates each emitted `.ferx` with the
  engine. Engine warnings are reported; engine **errors fail**.

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
