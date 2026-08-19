# ferx output gaps: ferx_gof prerequisites, missing fields, and the .lst equivalent

Three related topics in one document:

1. Confirmed answers to the two open ferx_gof() questions
2. Output fields that are missing or wrong in the current ferx_fit result
3. What NONMEM's `.lst` file contains and what the ferx equivalent would be

---

## 1. Confirmed answers to the ferx_gof() open questions

### 1a. Does $sdtab include a CMT column for per-CMT models?

**Confirmed absent -- and fixed in ferx-core PR#172.**

The original Rust code built exactly:

```
ID, TIME, DV, [CENS if any BLOQ rows], [OCC if any IOV], PRED, IPRED,
CWRES, IWRES, EBE_OFV, N_OBS
```

CMT was absent. PR#172 adds CMT when `any(obs_cmts != 1)`, following the same
conditional-inclusion pattern as CENS and OCC. New column order:

```
ID, TIME, DV, [CENS], [OCC], [CMT], PRED, IPRED, CWRES, IWRES, EBE_OFV, N_OBS
```

**Consequence for ferx_gof():** Resolved once PR#172 merges. The "Live question
#1" in ferx-output-advantage.md section 3g is answered: CMT will be present for
multi-endpoint models.

**Remaining ferx_gof() constraint:** Until PR#172 merges, ferx_gof() must accept
`data` as a required argument for per-CMT models, or detect `ErrorSpec::PerCmt`
and refuse to plot without data.

### 1b. Does $eta_normality row order match $ebe_etas column order?

**Yes, and it is name-based — safe to rely on.**

`.ferx_compute_eta_normality()` in `fit.R` (line 1704):

```r
eta_cols <- setdiff(names(ebe_etas), "ID")
do.call(rbind, lapply(eta_cols, function(col) {
  data.frame(eta = col, W = ..., p_val = ..., flag = ...)
}))
```

`eta_normality$eta` contains the exact column names from `ebe_etas` in the same
order as `names(ebe_etas)[-1]`. Matching is by name, not fragile position.

**Correct usage in ferx_gof():**

```r
# Safe — match by name:
for (eta_col in names(ebe_etas)[-1]) {
  norm_row <- fit$eta_normality[fit$eta_normality$eta == eta_col, ]
  # annotate histogram of ebe_etas[[eta_col]] with norm_row$W, norm_row$p_val
}
# NOT safe — do not assume row i corresponds to column i+1 by position alone
```

---

## 2. Missing or wrong fields in $sdtab

### 2a. [FIXED in PR#172] Original subject ID not preserved

`sdtab()` used `ids.push(si as f64 + 1.0)` -- the sequential loop index. Fixed
in PR#172 to use `sr.id.parse::<f64>().unwrap_or(si as f64 + 1.0)`.

For non-numeric IDs (should not occur with valid NONMEM data), the fallback now
emits a diagnostic to `FitResult.warnings`. The silent regression path from the
ultra-review is addressed.

The corresponding R-side workaround in persist.R (which remapped sequential
sdtab IDs via ebe_etas) has been removed in PR#107. This confirms that
`ebe_etas$ID` was already using the original subject ID -- the bug was isolated
to `sdtab()` only. The ebe_etas merge issue from ultra-review finding #3 is
therefore less severe than feared: ebe_etas was not broken.

### 2b. CMT absent from sdtab (see 1a above)

### 2c. No MDV flag in sdtab

Observations with `MDV=1` in the original data are excluded from estimation but
some analyses (e.g., visual overlays with predicted profiles at all timepoints)
need to know which rows were MDV. Currently sdtab contains only the observation
rows that contributed to the likelihood. Rows where `MDV=1` are dropped.

For standard GOF this is fine. For producing individual profile plots that include
predicted concentrations at scheduled timepoints where no sample was taken, MDV
rows are needed.

**Impact:** Low for GOF; medium for profile plots and VPC pre-processing.

### 2d. No individual predicted profile over time

`$sdtab` contains PRED and IPRED only at observation timepoints. For producing
the standard "individual fits" plot (observed DV vs time with individual predicted
curve), you need IPRED at a fine time grid — not just at observed timepoints.

This is a simulation task (`ferx_simulate(fit, data, n_sim=1)` at fine times
gives the curve), but it's not part of the fit result itself. Currently the
user must run a separate simulation step to get smooth individual curves.

**Proposal:** `ferx_ipred_profile(fit, data, times=NULL)` that returns a data
frame of IPRED at user-specified or auto-generated dense timepoints per subject.
Enables the individual fits plot without a full simulation pass.

---

## 3. Missing fields in ferx_fit that would be useful

### 3a. [FIXED in PR#172/107] Initial estimates not stored

`fit$theta_init`, `fit$omega_init`, `fit$sigma_init` now populated from
`stage_params` *after* the NCA warm-start block (so init values reflect what
the optimiser actually started from, not the raw .ferx parameters before NCA).
Backward-compatible: missing keys in older `.fitrx` bundles deserialise to
empty vecs; omega_init fallback is `omega.clone()` (not zeros).

### 3b. [FIXED in PR#172/107] Final gradient vector not exposed

`fit$final_gradient` now populated at the best-OFV point inside the NLopt
closures, gated on the global `best_seen` accumulator. NULL for
BOBYQA/BFGS/GN/SAEM. `ferx_runlog()` prints max|gradient| and flags
components exceeding `gradient_tol`.

### 3c. [MOSTLY RESOLVED] Data summary not in fit result

Earlier analysis was wrong: `fit$n_subjects` and `fit$n_obs` were already
present in the fit result before these PRs (lib.rs lines 1320-1321 pre-existing).
`fit$obs_time_range` was added in PR#172. `fit$covariate_names` was added in
PR#107 (lib.rs line 1490).

`ferx_runlog()` [DATA SUMMARY] therefore shows N subjects, N observations,
and time range -- all three of the primary data summary lines that open every
NONMEM .lst file.

Still missing: `n_doses` and `obs_per_subj` (min/median/max per subject).
These are secondary and not blocking any current planned function.

### 3d. [FIXED in PR#172/107] Model text not stored in fit result

`fit$model_text` (character or NULL) now populated from `FitResult.model_text`.
Added to ferx-core (not ferx-r only as originally planned), populated via a
second `read_to_string` in all three file-based entry points. `load_fit()`
reads it from the zip's `model.ferx` entry. NULL for in-memory fits.
`ferx_runlog()` omits the [MODEL] section gracefully when NULL.

---

## 4. The .lst equivalent in ferx

### What NONMEM's .lst contains

A NONMEM `.lst` file is the primary output artifact — the regulatory run
documentation. It contains, in order:

```
1.  Run identification: NONMEM version, licence, date/time
2.  $PROBLEM description (the model label)
3.  Full input control stream verbatim (the complete .ctl/.mod text)
4.  Data summary: N records, N individuals, N observations, key column stats
5.  $THETA initial values and bounds
6.  Per-iteration estimation log:
      ITERATION  OFV  THETA(1) THETA(2) ... (wide fixed-width table)
7.  Minimisation status line:
      "MINIMIZATION SUCCESSFUL" or specific failure + reason
8.  NPARAMETERS (number of free parameters)
9.  Final THETA table with SE and %RSE
10. Final OMEGA lower triangular with SE
11. Final SIGMA with SE
12. Covariance step status and computation notes
13. Correlation matrix of all estimates
14. Eigenvalues of the correlation matrix
15. Condition number
16. Computation summary: CPU time, N function evaluations
```

### What ferx currently has as equivalents

| .lst item | ferx equivalent | Gap? |
|---|---|---|
| NONMEM version | `$ferx_version` | Version field exists |
| Date/time | not stored | **Missing** |
| $PROBLEM / model label | `$model_name` | Present |
| Full model text | not in fit result | **Missing** (see 3d) |
| Data summary stats | not in fit result | **Missing** (see 3c) |
| Initial estimates | not stored | **Missing** (see 3a) |
| Iteration log | `$trace_path` CSV | Opt-in only; always-on needed |
| Convergence status | `$converged` | Present |
| N parameters | `$n_parameters` (struct) | Check if R-exposed |
| Final THETA + SE | `$theta`, `$se_theta` | Present |
| Final OMEGA + SE | `$omega`, `$se_omega` | Present |
| Final SIGMA + SE | `$sigma`, `$se_sigma` | Present |
| Correlation matrix | `ferx_cor_matrix()` | Present |
| Eigenvalues | `$eigenvalues` | Present |
| Condition number | `$condition_number` | Present |
| CPU / wall time | `$wall_time_secs` | Present |
| N function evals | from trace or fit | Check exposure |
| Final gradient | not exposed | **Missing** (see 3b) |

### ferx_runlog() — proposed function

A `ferx_runlog(fit, output=NULL)` function in ferx-r that assembles the above
into a plain-text document formatted similarly to NONMEM's `.lst`:

```
================================================================
  ferx run log
  Model:   warfarin.ferx
  Data:    warfarin.csv
  Method:  FOCEI
  ferx:    0.1.5
  Date:    2026-06-02 14:23:11
================================================================

[MODEL]
# One-compartment oral PK model (warfarin)
...complete .ferx text...

[DATA SUMMARY]
  Subjects:          32
  Observations:      251  (min 3 / median 8 / max 14 per subject)
  Doses:             32

[INITIAL ESTIMATES]
  THETA:  TVCL=0.134  TVV=8.1  TVKA=1.0
  OMEGA:  ETA_CL=0.07  ETA_V=0.02  ETA_KA=0.40
  SIGMA:  PROP_ERR=0.01

[ESTIMATION]
  Iterations:    47
  OFV final:    -842.317
  Status:       MINIMIZATION SUCCESSFUL
  Gradient method: Enzyme AD

[FINAL ESTIMATES]
  THETA:
    TVCL   0.1338  SE=0.00821  RSE=6.1%
    TVV    8.093   SE=0.412    RSE=5.1%
    TVKA   1.021   SE=0.147    RSE=14.4%
  OMEGA (variance):
    ETA_CL 0.0694  SE=0.0189   RSE=27.2%
    ...
  SIGMA:
    PROP_ERR 0.00987 SE=0.000821 RSE=8.3%

[COVARIANCE STEP]
  Status:         computed
  Condition number: 124.3
  Eigenvalues:    0.142  0.387  0.841  1.203  ...

[DIAGNOSTICS]
  ETA shrinkage:  ETA_CL=12.3%  ETA_V=8.7%  ETA_KA=31.4%
  EPS shrinkage:  8.2%
  DW statistic:   1.98 (no autocorrelation)
  ETA normality:
    ETA_CL: W=0.973  p=0.612
    ETA_V:  W=0.981  p=0.843
    ETA_KA: W=0.891  p=0.021 [!]

[RUNTIME]
  Wall time:   12.4 seconds
  Gradient:    Enzyme AD (inner) / Enzyme AD (outer)
================================================================
```

**Status: ferx_runlog() is IMPLEMENTED in ferx-r PR#107 (R/runlog.R, 408 lines).**

`ferx_runlog(fit, gradient_tol = 0.01, verbose = TRUE)` produces sections:
- Run header (model name, method, ferx version, date)
- Model file (fit$model_text; omitted for in-memory fits)
- Data summary (N subjects, N obs, time range -- all available)
- Parameter estimates (INITIAL/FINAL/SE/%RSE; handles FIXED parameters)
- Estimation settings (optimizer, maxiter, gtol, BLOQ method, NCA warm-start,
  seeds, covariate names) -- this section has NO NONMEM equivalent; it is a
  genuine advantage of ferx_runlog() over NONMEM .lst
- Objective function (OFV, AIC, BIC, convergence flag)
- Covariance step (condition number with flag if >1000, eigenvalues)
- Diagnostics (ETA/EPS shrinkage, Durbin-Watson, Shapiro-Wilk per ETA with [!])
- Final gradient (max|gradient|, per-component check, convergence warning)
- Runtime (wall time, gradient method inner/outer)

**Tests:** `tests/testthat/test-runlog.R` (360 lines, 32 test cases). The PR
body incorrectly said "no new automated tests" -- tests exist and are
comprehensive. Covers: return type, all section headers, parameter table
layout, SE/%RSE, INITIAL column, covariance, estimation settings, diagnostics,
gradient tolerance, wrapping at 8 components, graceful degradation with NULL
fields, and pure ASCII output.

**Non-ASCII fix:** First commit used em-dashes in runlog.R -- CLAUDE.md
violation. Fixed in second commit.

**What ferx_runlog() still lacks vs NONMEM .lst:**
- Per-iteration trace (NONMEM .lst shows each iteration's OFV + thetas;
  ferx has $trace_path CSV but runlog does not embed it)
- n_doses and obs_per_subj in [DATA SUMMARY]
- File output: no `output="path.lst"` argument; verbose=TRUE prints to console,
  verbose=FALSE returns string silently

---

## 5. Current status and remaining work

| Item | Status |
|---|---|
| sdtab ID bug | FIXED (PR#172) |
| CMT in sdtab | FIXED (PR#172) |
| Initial estimates | FIXED (PR#172/107) |
| Model text in fit | FIXED (PR#172/107) |
| Final gradient | FIXED (PR#172/107) |
| Data summary (n_subjects, n_obs, etc.) | PARTIAL -- obs_time_range only |
| ferx_runlog() | IMPLEMENTED (PR#107) but incomplete [DATA SUMMARY] |
| ferx_gof() | Not yet built; prerequisites now met once PRs merge |

**Remaining prerequisite for a complete ferx_runlog():**
Complete Fix 4 -- add n_subjects, n_obs, n_doses, obs_per_subj to FitResult.
Estimated 4-5 days across both repos. No PR open.

**Next build priority:** ferx_gof() -- sdtab ID and CMT prerequisites are now
met. The "Live questions" in ferx-output-advantage.md section 3g are answered:
CMT will be in sdtab for multi-endpoint models; eta_normality is name-matched.
ferx_gof() can be scoped and built without further ferx-core changes.
