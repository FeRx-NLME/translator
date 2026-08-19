# ferx-core / ferx-r output fixes

Prerequisite fixes identified in `ferx-output-gaps.md`. All are small, targeted
changes in ferx-core and/or ferx-r. Ordered by urgency — items 1 and 2 must
land before `ferx_gof()` is written; items 3-6 must land before `ferx_runlog()`.

**Status as of 2026-06-02:** All six fixes are in flight across two open PRs:
- ferx-core PR#172 covers Fixes 1-6 (Rust side)
- ferx-r PR#107 covers Fixes 1-6 (R side) + implements ferx_runlog()
Neither PR has merged yet. Fix 4 is only partially implemented (see below).

---

## Fix 1 -- sdtab uses sequential index instead of original ID [ferx-core]

**File:** `ferx-core/src/io/output.rs`, `sdtab()` function, line ~488

**Current:**
```rust
ids.push(si as f64 + 1.0);   // sequential index
```

**Fix (in PR#172):**
```rust
ids.push(sr.id.parse::<f64>().unwrap_or(si as f64 + 1.0));
```

For non-numeric IDs (should not occur with valid NONMEM data), the fallback now
also emits a diagnostic to `FitResult.warnings` listing the affected subject
count and an example ID. This replaces the silent fallback that was flagged
in the ultra-review.

Tests added in PR#172: `sdtab_id_column_uses_subject_id_not_loop_index` and
`sdtab_id_column_falls_back_for_non_numeric_ids`.

**Impact:** Breaks any downstream join from `$sdtab` to original data by ID for
datasets with non-consecutive numeric IDs. Blocks `ferx_gof(fit, data)` ETA vs
covariate plots. Must fix before `ferx_gof()` ships.

**Status: RESOLVED in PR#172 / PR#107.**

---

## Fix 2 -- CMT absent from sdtab for per-CMT models [ferx-core]

**File:** `ferx-core/src/io/output.rs`, `sdtab()` function

**Current:** No CMT column is built. `obs_cmts: Vec<usize>` is available on
`Subject` but never used in `sdtab()`.

**Fix (in PR#172):** Conditionally include CMT when `any(obs_cmts != 1)` -- same
conditional-inclusion pattern already used for CENS and OCC. Column order in
sdtab: `ID, TIME, DV, [CENS], [OCC], [CMT], PRED, ...`

Tests added in PR#172: `sdtab_cmt_column_present_for_multi_cmt` and
`sdtab_cmt_column_absent_for_single_cmt`.

**Impact:** Without CMT, `ferx_gof()` cannot stratify residual plots by endpoint
for joint PK-PD models. Plots mixing CMT=1 (PK) and CMT=2 (PD) residuals are
actively misleading.

**Status: RESOLVED in PR#172 / PR#107.**

---

## Fix 3 -- Initial estimates not stored in fit result [ferx-core + ferx-r]

**What:** Add `theta_init`, `omega_init`, `sigma_init` to `FitResult` -- the
parameter values at the start of estimation after any NCA warm-start.

**Implementation (in PR#172/107):**
- Captured from `stage_params` *after* the NCA block (not from `init_params`
  before it -- the first commit had this wrong; second commit corrected it).
  This ensures `theta_init` reflects what the optimiser actually started from.
- Wired through `FitWire` with `#[serde(default)]` for backward-compatible
  loading of older `.fitrx` bundles (missing keys deserialise to empty vecs).
- `omega_init` backward-compat fallback in `wire_to_fit_result` is `omega.clone()`
  (not `DMatrix::zeros`) to avoid breaking Cholesky-based consumers.
- R side: `fit$theta_init`, `fit$omega_init` (matrix + `fit$omega_init_dim`),
  `fit$sigma_init`.

**Needed for:** `ferx_runlog()` (the `.lst` equivalent). Also useful for
diagnosing convergence -- how far did the optimiser move from the starting point?

**Status: RESOLVED in PR#172 / PR#107.** No new tests for this fix (PR#107
notes "pure output formatting; underlying fields tested on ferx-core side").

---

## Fix 4 -- Data summary not in fit result

**Status: MOSTLY RESOLVED across pre-existing code + PR#172 / PR#107.**

Earlier analysis was wrong about the scope of this gap. Actual state:

| Field | Status |
|---|---|
| `fit$n_subjects` | ALREADY PRESENT (pre-existing in lib.rs line 1320) |
| `fit$n_obs` | ALREADY PRESENT (pre-existing in lib.rs line 1321) |
| `fit$obs_time_range` | ADDED in PR#172 / PR#107 |
| `fit$covariate_names` | ADDED in PR#107 (lib.rs line 1490) |
| `fit$n_doses` | Still missing; not in FitResult |
| `fit$obs_per_subj` (min/median/max) | Still missing; not in FitResult |

**Consequence for ferx_runlog():** The [DATA SUMMARY] section in runlog.R
shows N subjects, N observations, and time range -- all available. What it
cannot show: n_doses and per-subject obs distribution. This is secondary
data that most NONMEM .lst files show but that is not required for the
primary regulatory purpose.

**Consequence for ferx_compare():** `fit$n_obs` is already present, so the
comparison table CAN surface the N that BIC was computed with. The
"ferx_compare() depends on Fix 4" dependency is weaker than stated -- the
only remaining prerequisite is the AIC/BIC DOF fix for FIX thetas.

**Remaining gap:** n_doses and obs_per_subj. Estimated 1-2 days in ferx-core
(datareader.rs counts) + 0.5 day in ferx-r. Low priority -- not blocking
any current planned function.

---

## Fix 5 -- Model text not stored in fit result

**What:** Store the full `.ferx` text in the fit result so the fitted model is
self-contained and traceable without the original file.

**Implementation (in PR#172/107):** Moved to ferx-core, not ferx-r only as
originally planned. `FitResult.model_text: Option<String>` is populated via a
second `read_to_string` in all three file-based entry points
(`run_model_with_data_inits`, `run_model_simulate`, `fit_from_files`).
`load_fit()` sets it from the zip's `model.ferx` entry. `save_fit()` uses
`result.model_text.as_deref().unwrap_or(model_source)` so old callers are
unaffected.

R side: `fit$model_text` (NULL for in-memory fits, character for file-based).
`ferx_runlog()` gracefully omits the [MODEL] section when it is NULL.

**Needed for:** `ferx_runlog()`. Also for `ferxtranslate` round-trip inspection
(translate -> fit -> confirm what was actually run).

**Status: RESOLVED in PR#172 / PR#107.**

---

## Fix 6 -- Final gradient vector not exposed [ferx-core + ferx-r]

**What:** At convergence of gradient-based methods (FOCE, FOCEI, GN), the final
gradient vector should be exposed. A non-zero final gradient indicates the
optimiser stopped before true convergence (MAXITER hit or tolerance too loose).

**Implementation (in PR#172/107):**
- `Arc<Mutex<Option<Vec<f64>>>>` shared state before NLopt objective closures
  (primary and SLSQP fallback). Both closures share one Arc so the final value
  is from whichever run found the better OFV.
- Stored inside `if let Some(g) = grad { … }` gated on `best_seen` (the global
  accumulator) -- NOT on the closure-local `state.best_ofv`, which the first
  commit had wrong and the second commit corrected. The gate ensures the stored
  gradient corresponds to the best-OFV point.
- `None` for BOBYQA, built-in BFGS, GN, SAEM, and trust-region.
- R side: `fit$final_gradient` (NULL or named numeric vector).
- `ferx_runlog()` prints `max|gradient|` and flags individual components
  exceeding `gradient_tol` (default 0.01). Wraps at 8 components per line.

**Needed for:** `ferx_runlog()` (NONMEM prints the final gradient). Also a
direct convergence quality diagnostic.

**Status: RESOLVED in PR#172 / PR#107.**

---

## Summary table

| Fix | Repo | Status | Remaining gap |
|---|---|---|---|
| 1. sdtab original ID | ferx-core | RESOLVED (PR#172) | none; non-numeric warning added |
| 2. CMT in sdtab | ferx-core | RESOLVED (PR#172) | none |
| 3. Initial estimates | ferx-core + ferx-r | RESOLVED (PR#172/107) | no R-side tests |
| 4. Data summary | ferx-core + ferx-r | PARTIAL (PR#172/107) | n_subjects, n_obs, n_doses, covariates still missing |
| 5. Model text | ferx-core (moved) | RESOLVED (PR#172/107) | none |
| 6. Final gradient | ferx-core + ferx-r | RESOLVED (PR#172/107) | no R-side tests |

**Remaining work after PRs merge:**

Fix 4 is the only partially unresolved item, and the remaining gap (n_doses,
obs_per_subj) is secondary -- not blocking any planned function.

PR#107 also exposed several fit result fields beyond the original fix plan:
`optimizer_label`, `outer_maxiter`, `outer_gtol`, `bloq_method_label`,
`inits_from_nca`, `n_starts`, `multi_start_seed`, `saem_seed`,
`sir_seed_used`, `is_seed`, `covariate_names`. These power the ESTIMATION
SETTINGS section of ferx_runlog().

**Test coverage:** ferx-r PR#107 includes `tests/testthat/test-runlog.R`
(360 lines, 32 test cases across Tiers 1-3). The PR body was wrong when it
said "no new automated tests." Tests cover: return type, all mandatory
section headers, parameter table content, SE/%RSE columns, INITIAL column,
covariance step, estimation settings, diagnostics, gradient tolerance,
wrapping, graceful degradation with NULL fields, and pure ASCII output.

**What ferx_runlog() does NOT yet cover:**
- Per-iteration log (trace). The .lst shows each iteration's OFV and thetas.
  ferx has $trace_path (CSV) but runlog doesn't include it. Would make the
  runlog very long; likely should be a separate ferx_trace() display.
- n_doses and obs_per_subj in [DATA SUMMARY].
- File output parameter (verbose=TRUE/FALSE controls console print; no
  `output="path.lst"` argument to write to disk).
