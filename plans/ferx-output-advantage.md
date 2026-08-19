# ferx output: what NONMEM has, what ferx has, what neither has

This plan maps the full output space of a pharmacometric estimation run across
three categories:

1. What NONMEM produces that ferx already matches or exceeds
2. What ferx already produces that NONMEM does NOT
3. What neither produces today — genuine gaps that ferx could close first

The goal is to identify where ferx can make a credible superiority claim and where
there is net-new capability worth building.

---

## 1. What NONMEM produces — and ferx's status

| NONMEM output | NONMEM file | ferx status |
|---|---|---|
| OFV (-2 log-likelihood) | `.ext` last row | `$ofv` — identical concept |
| Parameter estimates (THETA, OMEGA, SIGMA) | `.ext` last row | `$theta`, `$omega`, `$sigma` |
| Standard errors | `.ext` SE row | `$se_theta`, `$se_omega`, `$se_sigma` |
| %RSE per parameter | derived | computable from `se / estimate * 100`; `ferx_estimates()` shows this |
| Wald 95% CI | derived | computable from SE; `ferx_estimates()` |
| Full parameter covariance matrix | `.cov` | `$cov_matrix` |
| Parameter correlation matrix | `.cor` | `ferx_cor_matrix()` |
| Inverse covariance matrix | `.coi` | not exposed separately; rarely needed |
| Eigenvalues of correlation matrix | `.lst` body | `$eigenvalues` |
| Condition number | `.lst` body | `$condition_number`; auto-flagged if > 1000 |
| Convergence status | `.lst` header | `$converged`, `$covariance_status` |
| Number of iterations | `.lst` | `$n_iterations` (outer iterations; field confirmed in fit.R) |
| Individual ETA estimates (EBE) | `.phi` | `$ebe_etas` |
| Individual OFV contributions | `.phi` | `$sdtab$EBE_OFV` |
| Population predictions (PRED) | `$TABLE` | `$sdtab$PRED` |
| Individual predictions (IPRED) | `$TABLE` | `$sdtab$IPRED` |
| Conditional weighted residuals (CWRES) | `$TABLE` | `$sdtab$CWRES` |
| Individual weighted residuals (IWRES) | `$TABLE` | `$sdtab$IWRES` |
| ETA shrinkage | `.lst` | `$shrinkage_eta` |
| EPS shrinkage | `.lst` | `$shrinkage_eps` |
| Iteration-by-iteration parameters | `.ext` body | `$trace_path` (CSV file) |
| Wall time | not reported | `$wall_time_secs` |

**Result:** ferx matches every meaningful NONMEM output. The only gap is
presentation — NONMEM users expect flat files with specific column names.
`ferx_estimates()` and `$sdtab` cover this for programmatic use; a `write_nm_tables()`
convenience function that writes NONMEM-style `.tab`/`.sdtab`/`.patab` flat files
is the missing UX layer.

One deliberate non-match: WRES (first-order weighted residuals) appear in NONMEM
$TABLE outputs. WRES is a FO approximation residual with known limitations. CWRES
(conditional weighted residuals) is strictly superior. ferx has CWRES and IWRES;
WRES is not planned.

---

## 2. What ferx already has that NONMEM does NOT

These are existing ferx_fit outputs with no NONMEM equivalent.

### 2a. AIC and BIC (automatic)

`$aic` and `$bic` are returned on every fit.

NONMEM reports only OFV. Practitioners compute AIC/BIC manually as:
`AIC = OFV + 2 * n_params`, `BIC = OFV + log(N) * n_params`. This arithmetic
is trivial but every NONMEM shop reinvents it via post-processing scripts or
tools like Xpose/PsN. ferx returns them directly.

**Value:** Immediate model comparison without external tools. Enables a
`ferx_compare(fit1, fit2)` function that shows delta-OFV, delta-AIC, delta-BIC in
one call.

### 2b. ETA normality test (Shapiro-Wilk, automatic)

`$eta_normality`: a data frame with columns `eta`, `W`, `p_val`, `flag`. A `[!]`
flag appears when `p < 0.05`.

NONMEM does not test ETA normality. Practitioners diagnose this visually via
histograms or Q-Q plots in Xpose, which requires exporting the `.phi` file and
loading into R. ferx runs the test as part of every fit and flags violations
automatically.

**Value:** Catches misspecified random-effect distributions (e.g., a bimodal ETA
distribution suggesting a mixture population) without any external step. This is
a genuine model misspecification alarm that NONMEM does not have.

### 2c. Durbin-Watson autocorrelation of IWRES

`$dw_statistic`: pooled Durbin-Watson statistic. Values near 2.0 = no
autocorrelation; < 1.5 = positive (missing structural dynamics); > 2.5 = negative
(over-parameterisation). `$iwres_lag1_r`: lag-1 Pearson correlation.

NONMEM has no equivalent. Practitioners detect autocorrelation visually via
IWRES vs TIME plots in Xpose/ggPMX. ferx quantifies it automatically.

**Value:** A single number that flags models where the residuals are not
exchangeable in time — the most common structural misspecification signal. Can
be added to a model comparison table as a goodness-of-fit metric.

### 2d. SIR confidence intervals (profile-likelihood equivalent)

`$sir_ci_theta`, `$sir_ci_omega`, `$sir_ci_sigma`: 95% CIs derived via sampling
importance resampling, which approximates the profile-likelihood confidence
intervals. These are asymptotically more correct than Wald intervals (which assume
the likelihood surface is quadratic around the mode).

NONMEM gives only Wald SEs. Bootstrap CIs require PsN bootstrap (hours of compute
time on a cluster, external tool). SIR in ferx takes seconds on top of the fit.

**Value:** Statistically superior uncertainty quantification at negligible
additional cost. Especially important for OMEGA parameters where Wald intervals
systematically underperform.

### 2e. Importance sampling marginal log-likelihood

`$importance_sampling$minus2_log_likelihood`: a Monte-Carlo estimate of the true
marginal log-likelihood, not the FOCE approximation. The `mc_standard_error` field
quantifies its precision.

NONMEM's OFV is an FOCE approximation to -2 log L. The true marginal likelihood
requires integrating out ETAs, which NONMEM only approximates. ferx's `imp` stage
estimates this directly.

**Value:** Gold-standard model comparison metric. When `imp` is run after FOCEI,
the importance-sampling OFV is a more valid basis for likelihood ratio tests and
information criteria, especially for models with high shrinkage or non-normal ETAs.

### 2f. Exact gradient (AD precision)

`$gradient_used = "ad"` signals that Enzyme AD was used throughout.

NONMEM uses central finite differences for all gradient computations. The gradient
error is O(h^2) where h is the step size, and step-size selection is heuristic.
ferx's AD gradient is exact to floating-point precision. This means:
- Covariance matrix estimates are more accurate (the R and S matrices in the
  sandwich estimator use the exact gradient)
- Convergence is to a tighter tolerance
- Standard errors are more trustworthy, especially for poorly identified parameters

**Value:** Directly translates to better SE estimates. For models near
non-identifiability (high condition number), this difference matters most.

### 2g. ETA lognormality metadata

`$eta_log_transformed`: logical vector indicating which ETAs are lognormally
parameterised. This drives `$omega_param_corr` (parameter-level correlation for
block omega using the bivariate lognormal formula, not the eta-level approximation).

NONMEM does not track this; it reports omega on the eta scale only. When omega is
block and ETAs are lognormal, the reported correlation is not the correlation of
the underlying PK parameters — ferx corrects for this.

### 2h. Gradient method transparency

`$gradient_method_inner`, `$gradient_method_outer`: exact report of what method
was used. NONMEM never tells you; it always uses FD silently.

**Value:** Reproducibility and trust. When ferx uses AD, the user knows their
gradient is exact. When it falls back to FD (Form B/C ODE models), the user knows
why and can reason about the quality of the result.

### 2i. Estimation settings in the run log (ferx_runlog only)

`ferx_runlog()` (PR#107) produces an ESTIMATION SETTINGS section that NONMEM's
`.lst` does not have:

- Optimizer name (LBFGS, BOBYQA, trust-region, etc.)
- Maximum iterations and gradient tolerance used
- BLOQ handling method
- Whether NCA warm-start was applied and which algorithm
- Multi-start count and random seeds for every stochastic stage (SAEM, SIR,
  importance sampling, multi-start)
- Covariate column names in the data

NONMEM's `.lst` gives no information about optimizer settings, seeds, or
warm-start strategy. Exact reproducibility of a NONMEM run requires the
control stream and the PsN run script; ferx embeds this information in the
standard run log.

**Value:** A ferx run is self-documenting. Any future re-run can recover the
exact settings from the `.lst` equivalent without a separate configuration
file. This is a direct advantage for regulatory submissions and GxP audit trails.

---

## 3. What neither NONMEM nor ferx has — genuine gaps ferx could close first

These are capabilities that practitioners currently get from external tools (PsN,
vpc R package, Xpose, nlme tools) or do not have at all.

### 3a. Prediction-corrected VPC (built-in)

A visual predictive check (VPC) simulates from the fitted model, bins by
prediction-corrected observations, and overlays observed vs simulated percentiles.
It is the standard tool for model adequacy evaluation.

Neither NONMEM nor ferx produces VPCs. Practitioners run PsN vpc (NONMEM), then
`vpc::vpc()` in R. ferx has `ferx_simulate()` — the raw material is there.

**What to build:** `ferx_vpc(fit, data, n_sim=1000L, bins=10L)` that:
1. Simulates `n_sim` replicates via `ferx_simulate()`
2. Applies prediction correction (divide simulated and observed DV by PRED)
3. Bins by PRED, computes 5th/50th/95th quantiles of simulated and observed
4. Returns a ggplot2 plot object (or a data frame for custom plotting)

**Effort:** 2-3 weeks. The simulation infrastructure is already there.
**Payoff:** This is the single most-used diagnostic in NCA/NLME modelling. Having
it built-in, native, fast (ferx_simulate is Rust), is a major differentiator.

### 3b. NPDE (Normalized Prediction Distribution Errors, built-in)

NPDE is a simulation-based residual that does not depend on the FOCE approximation.
It is theoretically superior to CWRES for detecting model misspecification.

Currently available via the `npde` R package, but requires: running the model,
simulating, computing NPDEs manually. Neither NONMEM nor ferx computes NPDEs
natively.

**What to build:** `ferx_npde(fit, data, n_sim=1000L)` that:
1. Simulates `n_sim` replicates
2. For each observed DV, computes the empirical quantile in the simulation
   distribution (the prediction distribution function, PDF)
3. Applies the probit transform to get NPDE ~ N(0,1)
4. Returns a data frame with NPDE and a diagnostic plot

**Effort:** 1-2 weeks once VPC simulation infrastructure exists.
**Payoff:** NPDE is the simulation-based residual recommended in FDA guidance for
confirmatory modelling. Built-in NPDE would be unique among NLME engines.

### 3c. Automatic model comparison table

When a user runs multiple models (base, + covariate 1, + covariate 2, etc.),
they need a comparison table with:
- Model label
- OFV, delta-OFV vs base, p-value (chi-squared test)
- AIC, BIC
- Number of parameters
- Convergence status
- Condition number
- ETA shrinkage summary (max shrinkage)
- DW statistic

**What to build:** `ferx_compare(fit_list, labels=NULL)` that:
1. Takes a named list of `ferx_fit` objects
2. Produces the comparison table as a data frame
3. Optionally prints in NONMEM-style or knitr-ready format

**Effort:** 1 week. All the inputs already exist in `ferx_fit` fields.
**Payoff:** Replaces PsN `runrecord`, Xpose model comparison, and manual Excel
tables that are the current standard workflow.

### 3d. Per-subject influence diagnostics (case-deletion)

Identifies subjects that disproportionately influence parameter estimates. In
NONMEM this requires manually re-running the model N times (once per subject
excluded). External tools like Perl scripts automate this but it remains slow.

With ferx's speed, case-deletion could be automated:

**What to build:** `ferx_influence(fit, data, model)` that:
1. Re-fits the model N times, each time dropping one subject
2. Returns a data frame: subject ID, delta-OFV, delta-theta, delta-omega per
   parameter
3. Flags subjects whose removal changes any parameter estimate by > X%

**Effort:** 3-5 weeks (design + parallelisation). Requires ferx to be fast enough
that N re-fits is affordable. AD-based ferx on analytical PK models fits in
seconds — 100 subjects = 100 fits = a few minutes.
**Payoff:** Influence diagnostics are standard in regulatory submissions but
rarely done because the compute cost is prohibitive with NONMEM. ferx speed
makes them routine.

### 3e. Bootstrap parameter distribution (subject-level resampling)

Resample subjects with replacement, re-fit, collect the distribution of parameter
estimates. This is the frequentist gold standard for uncertainty quantification.

PsN bootstrap runs this but it requires a NONMEM license and typically runs
overnight on a cluster. ferx's SIR approximates this well for well-identified
models, but true bootstrap makes a different (weaker) assumption.

**What to build:** `ferx_bootstrap(fit, data, model, n_boot=200L, ncores=NULL)` that:
1. Resamples subjects with replacement
2. Re-fits on each bootstrap dataset (parallelised across cores)
3. Returns the full bootstrap distribution as a data frame
4. Computes 95% bootstrap CIs

**Effort:** 4-6 weeks. Requires parallel execution infrastructure (could use
`parallel` or `future` package). The main design question is how to handle
datasets where a subject appears multiple times.
**Payoff:** True bootstrap is the benchmark that SIR is compared to. Offering
both is a strong story.

### 3f. Covariate screening table (forward/backward)

The standard covariate model building workflow is:
1. Run base model
2. For each candidate covariate × parameter combination, add the covariate and
   check delta-OFV
3. Accept or reject based on a significance threshold
4. Backward elimination

With the `ferx_ir` modifier API (Phase 5 of the direct-run plan), this can be
automated in R:

**What to build:** `ferx_covariate_screen(base_ir, data, covariates, params,
models=c("power","proportional","linear"), alpha_forward=0.05, alpha_backward=0.01)` that:
1. Runs the base model
2. For each covariate × param × model type: uses `ir_add_covariate()` to generate
   a modified ferx_ir, fits it, records delta-OFV
3. Returns the full covariate screening table
4. Optionally runs backward elimination

**Effort:** 3-4 weeks (builds on Phase 5 ferx_ir modifier API).
**Payoff:** Replaces PsN SCM (stepwise covariate modelling), which requires a
NONMEM license and separate Perl-based configuration files. A native R-based
covariate screening tool with ferx speed would be genuinely novel.

### 3g. ferx_gof() — standard goodness-of-fit panel

The most-used diagnostic after fitting a model is a 4-6 panel GOF page: DV vs
PRED, DV vs IPRED, CWRES vs TIME, CWRES vs PRED, and ETA distributions. Every
pharmacometric report contains one. Currently the user exports `$sdtab` and
`$ebe_etas` and builds plots manually or via Xpose.

**What ferx already provides for free:**

All data needed for the core GOF plots lives in the fit object with no
recomputation:

| Plot | Source |
|---|---|
| DV vs PRED / DV vs IPRED | `$sdtab` (PRED, IPRED, DV, CENS) |
| CWRES vs TIME / CWRES vs PRED | `$sdtab` (CWRES, TIME, PRED) |
| CWRES distribution + Q-Q | `$sdtab$CWRES` |
| ETA histograms + N(0,ω) density | `$ebe_etas` + `$omega` diagonal |
| ETA Q-Q plots | `$ebe_etas` |
| ETA vs covariate | `$ebe_etas` merged with original `data` by ID |

**What makes ferx_gof() different from Xpose:**

The tests are already computed. `ferx_gof()` can annotate directly:
- CWRES vs TIME gets the **DW statistic** and its interpretation (from `$dw_statistic`)
  as a subtitle, coloured green/orange/red
- ETA histograms get **Shapiro-Wilk W and p-value** from `$eta_normality`,
  with a flag indicator when `p < 0.05`
- ETA histograms get **shrinkage %** from `$shrinkage_eta` — high shrinkage (>30%)
  means the ETA is driven by the prior, not the data; this warrants annotation
- ETA vs covariate scatter plots get **eta_corr_test p-value** if available from
  the fit, auto-flagging significant correlations

No external tool does this integration. Xpose computes shrinkage from scratch;
ferx already has it.

**Interface:**

```r
ferx_gof(fit,
         data     = NULL,      # needed for ETA vs covariate plots only
         log_dv   = TRUE,      # log-log scale for DV vs PRED/IPRED (typical PK)
         annotate = TRUE,      # overlay DW, Shapiro-Wilk, shrinkage on plots
         type     = c("standard", "residuals", "etas", "covariates", "all"))
```

Returns a named list of ggplot2 objects with class `"ferx_gof"`. The `print`
method lays them out via patchwork if available, otherwise prints sequentially.
Individual plots accessible by name (`gof$dv_pred`, `gof$cwres_time`, etc.).

**Standard panel layout (type = "standard"):**

```
| DV vs PRED          | DV vs IPRED         |
| CWRES vs TIME       | CWRES vs PRED       |
| CWRES histogram     | ETA_k histograms... |
```

**Implementation notes:**

- ggplot2 and patchwork go into ferx-r `Suggests` (not `Imports`); checked at
  runtime. `ferx_plot_trace()` already uses **base graphics** — `ferx_gof()` is
  the first ggplot2 surface in ferx-r. This is a deliberate decision: ggplot2
  is opt-in, not required for the core engine.
- BLOQ rows (`CENS = 1` in `$sdtab`) are shown as downward-pointing triangles
  in DV vs PRED/IPRED plots and flagged in CWRES plots.
- Multi-CMT models (per-CMT error models) get plots faceted by CMT column
  automatically if `CMT` appears in `$sdtab`.
- IOV models: `$sdtab$OCC` is present; CWRES vs TIME stratified by occasion.
- For ETA vs covariate: ferx merges `$ebe_etas` with `data` by ID, then produces
  one scatter + loess per ETA × covariate combination. If `$eta_names` is
  available (it is), x-axis labels are correct without any user effort.

**Answered questions (verified from source code and PRs):**

1. CMT column in $sdtab: **Yes, for multi-CMT models** -- fixed in ferx-core
   PR#172. Column appears when `any(obs_cmts != 1)`, following the CENS/OCC
   conditional-inclusion pattern. Column order: ID, TIME, DV, [CENS], [OCC],
   [CMT], PRED, IPRED, CWRES, IWRES, EBE_OFV, N_OBS. ferx_gof() can facet by
   CMT from $sdtab directly once PR#172 merges.

2. $eta_normality row order: **Name-based, safe to rely on.** Confirmed from
   fit.R line 1704: iterates `setdiff(names(ebe_etas), "ID")` and stores the
   column name as `eta`. Match by `eta_normality$eta == eta_col`; never by
   positional index.

3. Shrinkage annotation threshold: **Default 30%.** This is the FDA guidance
   threshold below which ETAs are considered informative. Expose as
   `shrinkage_threshold` argument. Flag >50% in red (unreliable EBE), 30-50%
   in orange (caution). This can be hard-coded for MVP and made configurable.

**Effort:** 2-3 weeks for the standard panel (residuals + ETA histograms +
annotations). ETA vs covariate plots add 1 week. Full multi-CMT / IOV / BLOQ
handling adds another week. Total: 3-5 weeks for production quality. This is
where the plan-review correction applies — do not promise 1 week.

**Scope boundary:** `ferx_gof()` lives in **ferx-r**, not ferxtranslate.
It depends only on `fit$sdtab`, `fit$ebe_etas`, and the fit metadata — all
native ferx-r output. No translation layer involved.

---

## Priority ranking

| Item | Effort | Strategic value | Depends on | Repo | Status |
|---|---|---|---|---|---|
| `ferx_runlog()` .lst equivalent | -- | High -- regulatory traceability | Fixes 1-6 (ferx-core PR#172) | ferx-r | IMPLEMENTED in PR#107 (incomplete data summary) |
| `ferx_compare()` model comparison table | 1 week | High -- replaces manual workflow | Fix 4 complete (n_subjects, n_obs); AIC/BIC DOF fix | ferx-r | Not built |
| `ferx_gof()` standard GOF panel | 3-5 weeks | Very high -- first plot after every fit | PR#172 merged (CMT, sdtab ID) | ferx-r | Not built |
| `write_nm_tables()` flat file exporter | 1 week | Medium -- UX for NONMEM migrants | nothing | ferx-r | Not built |
| `ferx_vpc()` prediction-corrected VPC | 6-8 weeks MVP; 3 months full | Very high -- most-used diagnostic | ferx_simulate (exists) | ferx-r | Not built |
| `ferx_npde()` simulation-based residuals | 3-5 weeks | High -- regulatory gold standard | ferx_vpc simulation path | ferx-r | Not built |
| `ferx_covariate_screen()` | 3-4 weeks | Very high -- unique capability | ferx_ir modifier API (Phase 5) | ferxtranslate | Not built |
| `ferx_influence()` case-deletion | 3-5 weeks (linCmt only) | Medium -- niche but high-value | ferx speed; warm-start from fit | ferx-r | Not built |
| `ferx_bootstrap()` true bootstrap | 4-6 weeks | Medium -- SIR covers most of this | ID renaming logic; parallel | ferx-r | Not built |

**Notes on effort corrections from plan-review.md:**
- `ferx_vpc()` originally listed as 2-3 weeks. Revised: 6-8 weeks for an MVP
  (single endpoint, no IOV, no BLOQ stratification). Full production quality
  with binning, IOV, BLOQ, confidence bands on percentile lines is 3 months.
- `ferx_npde()` originally 1-2 weeks. Revised: 3-5 weeks because the
  decorrelation step (Brendel 2006) is not just empirical quantiles. The plan
  body description of the algorithm is wrong -- see ultra-review finding #1.
- `ferx_influence()` caveated: feasible only for analytical PK models.
  ODE models at 5-30 min per fit make case-deletion computationally intractable.
- `ferx_compare()` has one remaining prerequisite: the AIC/BIC FIX-theta
  DOF issue (plan-review §2.1) must be resolved or documented as a known
  limitation. `fit$n_obs` is already present (pre-existing field), so the
  N used for BIC can be surfaced in the comparison table without further
  ferx-core changes.

**Recommended build order:** `ferx_gof()` → `ferx_compare()` → `ferx_vpc()` →
`ferx_npde()` → `ferx_covariate_screen()`.

`ferx_gof()` is now unblocked (PR#172 fixes CMT and sdtab ID). Start here.
`ferx_compare()` depends on Fix 4 completion and the AIC/BIC DOF resolution;
do not build until those are settled. `ferx_vpc()` follows once the simulation
path is exercised by GOF. Together the first three make ferx sufficient for a
complete model-building report without any external tools.
