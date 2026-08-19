# Critical review of the three ferx plans

Covers: `direct-run-feasibility.md`, `ferx-output-advantage.md`,
`ferx-strategic-roadmap.md`. Each issue is rated by severity:
**[BLOCKER]** = will stop an implementation cold, **[MAJOR]** = significant
rework needed, **[MINOR]** = manageable but must not be ignored.

---

## Part 1: direct-run-feasibility.md

### 1.1 [BLOCKER] Headerless NONMEM data files

The plan says "Phase 1 data_file must be a named-header NONMEM-format CSV."
NONMEM does NOT require column headers. The entire $INPUT block exists precisely
because NONMEM data files are often positional and headerless. In real-world
clinical trial datasets, a significant fraction have no header row; $INPUT lists
the column names in order.

Phase 1 therefore fails silently for any dataset without a header — the ferx
data reader would read the first data row as column names, producing garbled
columns. The plan treats this as a known documentation limitation; it is actually
a blocking failure mode for a large class of real-world NONMEM runs.

**Fix:** Either (a) make Phase 1 explicitly require header detection and error
on headerless files before any fitting occurs, or (b) make data preprocessing
(reading $INPUT to add the header) a prerequisite to Phase 1, not Phase 2.

### 1.2 [BLOCKER] Silent misparse from nonmem2rx — no detection strategy

The plan notes: "nonmem2rx errors or silently returns a malformed UI object, and
the translator either crashes or produces wrong `.ferx`." The hard crash case is
handled by `tryCatch`. The silent misparse case — nonmem2rx parses without error
but extracts wrong expressions — is explicitly acknowledged but no mitigation is
proposed.

This is the most dangerous failure mode because a user gets a fit result and
wrong parameter estimates with no warning. The plan says "add it to the risk
register" but proposes nothing actionable.

**Fix:** Add a structural sanity check on the rxUi before translation: at minimum,
verify that the number of thetas in `iniDf` matches the number of distinct theta
indices referenced in `lstExpr`. A mismatch is a strong signal of silent
misparse. Also verify that the sigma count matches the error model. These checks
catch the most common nonmem2rx silent failures.

### 1.3 [MAJOR] ir_add_covariate() is not "30 lines of R"

Phase 5 describes `ir_add_covariate(ir, param="CL", covariate="WT", model="power",
reference=70)` as a simple helper. The function must:

1. Find the row in `ir$indiv_params` where `lhs == "CL"`
2. Parse the current `rhs` string (e.g., `"TVCL * exp(ETA_CL)"`) into an expression
3. Insert `* (WT/70)^0.75` at the correct location — multiplicatively before the
   exponential term
4. Re-serialise back to a string

Step 2 is non-trivial: `rhs` is stored as a **character string** in the IR, not
as an R expression object. The plan treats string manipulation as straightforward
but consider: what if rhs is `"TVCL * F_CRCL * exp(ETA_CL)"` (already has a
covariate)? Or `"(TVCL + DELTA_CL) * exp(ETA_CL)"` (additive theta structure)?
Or for linear: `CL = TVCL * exp(ETA_CL)` and WT goes in the TVcl part: `CL =
TVCL * (1 + BETA_WT * (WT - 70)) * exp(ETA_CL)`.

Getting the insertion point right for general expressions requires either a
mini-expression parser or storing `rhs` as an R language object (not a string)
throughout the IR. The plan skips this fundamental design decision.

**Fix:** Decide now: does `ferx_ir$indiv_params$rhs` become an R expression
object instead of a character string? If yes, the modification API is natural
but every emitter and consumer of that field changes. If no, document the
limitation explicitly: `ir_add_covariate()` only works for the standard
`TVCL * exp(ETA)` pattern.

### 1.4 [MAJOR] nlmixr2 model modification — unverified API assumption

Phase 8 Layer claims: "nlmixr2 already has update methods: `ini(model, ...)`,
`model(model, ...)` replaces model block expressions." This was not tested. The
R session was blocked by disk space before verification. The plan treats this as
proven but it is an assumption about an external API.

Specifically: does `rxode2::model(ui, cl <- tvcl * (WT/70)^0.75 * exp(eta.cl))`
create a modified rxUi with the new expression replacing the old one, or does it
append it? Does the modified rxUi still have the correct `$lstExpr` structure
that `rxui_to_ir()` expects?

**Fix:** Before Phase 5 is designed around this, verify the nlmixr2 model update
API with a concrete test. If it doesn't work as assumed, the "programmatic modify
→ translate → run" loop breaks at the first step.

### 1.5 [MAJOR] The "target audience" question is never answered

The plan ends with: "If the target audience is nlmixr2 users: 1→5→3→2→4→6. If
NONMEM users: 1→2→3→4→5→6." The document never commits. This is not a neutral
omission — without a target audience decision, the prioritisation is undefined and
the roadmap cannot be actioned. The two audiences have different data formats,
different tooling expectations, and different translation gaps that matter to them.

**Fix:** Decide the primary audience before implementation begins, state it
explicitly, and let that drive the phase ordering. An ambiguous roadmap produces
an ambiguous product.

### 1.6 [MINOR] WRES omission carries hidden regulatory risk

The output comparison table says "WRES is not planned" with the justification
that CWRES is strictly superior. This is technically true but FDA submission
guidelines and some older internal validation SOPs at CROs explicitly reference
WRES. Saying "not planned" without acknowledging the potential regulatory
objection is a happy path.

**Fix:** Note that WRES may be requested in specific regulatory contexts. Add a
low-priority ticket for WRES computation (it is trivially computable given CWRES
and the relationship between the FO and conditional approximations).

### 1.7 [MINOR] FO/LAPLACIAN hard-stop — what does the user see?

The plan says "hard-stop on FO" but doesn't specify the error message or what
information the user needs to proceed. A user running `nm_run("model.ctl",
"data.csv")` on a model with `METHOD=0` gets an error and no path forward.

**Fix:** The error message must:
1. Name the unsupported NONMEM option explicitly
2. Explain the statistical reason it cannot be silently upgraded
3. Suggest the nearest ferx equivalent (`focei` for `METHOD=1 INTERACTION`)
4. Point to documentation on expected numerical differences

### 1.8 [MINOR] Phase ordering rationale buried at the end

The phase ordering rationale is in the last section, after the recommended
starting point. A reader who stops before the end (common for long documents)
misses it. More critically: the rationale says "decide based on audience" but
the preceding "recommended starting point" section gives advice without knowing
the audience.

---

## Part 2: ferx-output-advantage.md

### 2.1 [MAJOR] AIC/BIC degrees of freedom — regulatory ambiguity

Section 2a presents AIC/BIC as a clear ferx advantage. Missing: how are the
degrees of freedom counted?

For a diagonal omega with 3 ETAs: n_params = n_theta + 3 (omega diagonal) + 1
(sigma) = n_theta + 4. Unambiguous.

For a block omega with 2 ETAs: the omega block has 3 parameters (2 variances + 1
covariance). Still unambiguous.

The ambiguity arises for fixed thetas: do FIX parameters count toward degrees of
freedom? In NONMEM they do not. In information criteria they should not. If ferx
includes FIX parameters in n_params, AIC/BIC are wrong for models with fixed
thetas. This matters because fixed thetas are common (e.g., allometric exponents
fixed at 0.75) and the model comparison will have wrong deltas.

**Fix:** Document explicitly how n_params is counted: only freely estimated
parameters (no FIX thetas, no Cholesky elements constrained to zero). Verify
this in the ferx-core output. If it's wrong, it needs to be corrected before
`ferx_compare()` is built on top.

### 2.2 [MAJOR] SIR CIs — "always better than Wald" is wrong

Section 2d: "SIR confidence intervals are asymptotically more correct than Wald
intervals." Then the strategic value says "superior uncertainty quantification."
This is too strong.

SIR requires the importance weights to have finite variance (the proposal must
have heavier tails than the target). The Gaussian proposal used in ferx's SIR is
fine when the posterior is unimodal and roughly Gaussian — which is exactly when
Wald intervals are also fine. When the posterior is multimodal or has heavy tails
(poorly identified models, high shrinkage), SIR with a Gaussian proposal can have
degenerate weights (low ESS → wide variance of the ESS estimate → unreliable CIs).

In the cases where Wald is worst (non-identifiability, high shrinkage), SIR with
a Gaussian proposal may also fail. The plan should say:
- SIR is better when the posterior is unimodal and `$sir_ess` is high (> 200)
- SIR degrades when ESS is low; the ESS must always be checked
- Profile likelihood or full Bayesian HMC is the correct answer when SIR ESS is low

### 2.3 [MAJOR] VPC complexity severely underestimated

Section 3a: `ferx_vpc()` described as 2-3 weeks with a 4-step algorithm. The
algorithm is superficially correct but misses:

**Binning:** "Bins by PRED" — auto-binning for irregular schedules is hard.
Equal-width bins ignore data density; quantile bins can have empty cells. The vpc
R package uses a stratified binning approach with visual bin boundaries. The
plan gives no thought to binning strategy.

**Multiple DVIDs/endpoints:** Joint PK-PD models have multiple DVs (CMT=1 is
PK, CMT=2 is PD). A VPC must be stratified by CMT/DVID. Not mentioned.

**IOV:** When models have IOV (kappa), the simulation must include the kappa
random effects. ferx_simulate handles this but the VPC must correctly stratify
by occasion. Not mentioned.

**BLOQ handling in simulation:** When observations are below LLOQ, simulated
values below LLOQ must be censored consistently. M3/M4 method-consistent VPCs
require specific treatment. Not mentioned.

**Stratification:** VPCs for covariate models require stratification by covariate
quartiles. Not mentioned.

**Confidence bands on the percentile lines:** The 95% PI on the 5th/50th/95th
percentiles (from simulating the simulation distribution) requires a nested
simulation. Not mentioned.

These aren't optional polish items — the vpc R package's equivalent code is
2000+ lines. A 2-3 week estimate is off by a factor of 3-4x for a production-
quality implementation. A minimum viable VPC (no stratification, no IOV, single
DVID) is feasible in 2-3 weeks but must be clearly scoped as such.

**Fix:** Either scope the MVP VPC explicitly (single endpoint, no IOV, no BLOQ,
simple bins), or revise the estimate to 6-8 weeks.

### 2.4 [MAJOR] NPDE algorithm described incorrectly

Section 3b describes NPDE as:
1. Simulate n_sim replicates
2. For each DV, compute the empirical quantile in the simulation distribution
3. Apply probit transform

This is the **simulation-based quantile** (PBPK) approach but NOT the correct
NPDE algorithm. The correct NPDE algorithm (Brendel et al. 2006; Mentré and
Escolano 2006) involves:

1. Simulate n_sim replicates **from the conditional distribution given the
   individual's covariate design** (not the marginal)
2. Compute the **decorrelated** residuals using the mean vector and variance
   matrix of the predicted distribution
3. Apply the probit to get N(0,1) residuals

The decorrelation step (step 2) accounts for within-subject correlation — without
it, the resulting residuals are not independent and the probit transform does not
yield N(0,1). This is why the `npde` R package is ~3000 lines, not 100.

**Fix:** Correct the algorithm description. Acknowledge that NPDE requires
within-subject simulation from the conditional distribution, not the marginal
unconditional simulation that `ferx_simulate()` currently does. Revise effort
to 3-5 weeks.

### 2.5 [MAJOR] ferx_influence() timing claim is only true for linCmt models

Section 3d: "100 subjects = 100 fits = a few minutes." This is true for
analytical PK models (linCmt → pk macro, milliseconds per fit). For ODE models,
a single fit can take 5-30 minutes even with ferx. 100 subjects × 10 min = 16+
hours. The plan doesn't distinguish model types.

Also: the re-fit after case-deletion starts from default initials each time. If
the original converged theta is not used as a warm start, many re-fits will not
converge. Warm starting from the original fit avoids this but is an implementation
detail that must be in the design.

**Fix:** Scope `ferx_influence()` explicitly for analytical PK models. Add a
warning for ODE models. Specify that re-fits use the original theta as warm start.

### 2.6 [MAJOR] ferx_bootstrap() — the ID collision problem is not resolved

The plan notes: "The main design question is how to handle datasets where a
subject appears multiple times." It then says "4-6 weeks" without resolving it.

This IS the hard problem for bootstrap: when subject 23 appears 3 times in the
bootstrap resample, their ID appears 3 times in the bootstrapped dataset, but
ferx groups by ID (so all their rows merge into one subject). The standard
solution is to reassign IDs in the bootstrapped dataset (ID_new = bootstrap_draw
index). But this requires knowing which subjects were resampled how many times
and renaming them consistently across all their rows. This is not a footnote — it
is the entire data preparation step.

**Fix:** Specify the ID renaming procedure explicitly. This is actually 1-2 weeks
of careful data manipulation code.

### 2.7 [MINOR] ferx_covariate_screen() — forward stepwise is known to be biased

Section 3f presents stepwise covariate screening (SCM) as the method to implement.
This is the standard approach in PMx (PsN SCM) but it is statistically well-known
to:
- Overfit when the number of candidates is large
- Produce inflated Type I error due to multiple testing
- Give different results depending on the order of evaluation

The plan doesn't acknowledge these limitations. If ferx builds SCM and people use
it uncritically, results will be statistically invalid for confirmatory purposes.

**Fix:** Add a note that forward SCM is an exploratory tool, not a confirmatory
one. Consider implementing LASSO-penalised covariate selection as a more
statistically principled alternative.

### 2.8 [MINOR] run log / regulatory submission artifact — missing

The output comparison table lists every NONMEM file (`.ext`, `.cov`, `.cor`,
`.phi`, `.tab`) and declares ferx has an equivalent for each. Missing: the `.lst`
run log, which contains:
- The complete input control stream (the model code as submitted)
- Summary statistics of the data
- Full estimation run log including each iteration's objective function
- Final gradient vector
- NONMEM version and FORTRAN compiler information

The `.lst` is a regulatory artifact required in every pharmacometric report and
NDA submission. Reviewers want to see it. ferx has no equivalent. The plan
should call this out as a gap requiring a `ferx_runlog()` function that generates
a GxP-compliant run summary document.

**Update (2026-06-02):** `ferx_runlog()` is now implemented in ferx-r PR#107.
The [DATA SUMMARY] section is currently incomplete (shows obs_time_range only;
n_subjects, n_obs, n_doses are missing because Fix 4 is partial). The function
exists but is not yet regulatory-complete for submissions.

---

## Part 3: ferx-strategic-roadmap.md

### 3.1 [BLOCKER] Non-Gaussian likelihoods require data format changes

Section 1 (count data, ordinal, binary) describes new error model types but
never addresses the data format prerequisite:

- **Poisson/NB**: DV must be integer. The current data reader reads DV as f64.
  For count data, non-integer DV values (e.g., 2.5 observations) should be an
  error. The lgamma computation for Poisson is defined for non-integer values
  (the Stirling extension) but produces meaningless likelihood for fractional
  counts.

- **Ordinal**: DV must be one of {0, 1, 2, ..., K-1} for K categories. The
  current per-observation likelihood assumes continuous DV. For ordinal, the
  likelihood is `P(DV = k) = logistic(SCORE - TH_{k}) - logistic(SCORE - TH_{k+1})`.
  The concept of CWRES has no meaning for ordinal outcomes.

- **TTE**: DV is an event indicator (0=censored, 1=event) and TIME is the
  observed event or censoring time. This is a completely different row structure
  from PK concentration data (where DV is concentration, TIME is sample time).

None of these can be shoehorned into the current `{ID, TIME, DV, EVID, AMT, CMT,
...}` data format. Each new likelihood type requires:
1. Data format specification for that endpoint type
2. Data reader changes in ferx-core
3. Parser changes for the new error model syntax
4. Gradient changes for the new likelihood
5. New diagnostic outputs (no CWRES/IWRES for non-Gaussian; what replaces them?)

The plan treats this as "3-4 weeks for Poisson." In reality Poisson is 6-10
weeks because the data format and diagnostics infrastructure must be built first.

### 3.2 [BLOCKER] TTE data structure is fundamentally incompatible with current format

Section 2 presents TTE as an extension of the existing ODE machinery. The data
structure issue is more severe than for count data:

NONMEM TTE data uses either:
- **Single event**: one row per subject with DV=0 (censored) or DV=1 (event),
  TIME=event time. EVID=0 for observations, EVID=1 for doses.
- **Repeated TTE**: multiple rows per subject, each with (start_time, end_time,
  event_indicator).

The survival function integral `S(t) = exp(-integral_0^t h(u) du)` must be
computed over the full observation window per subject. For a subject followed for
365 days with drug concentrations varying over time (because they received doses),
the hazard function depends on the PK trajectory, which requires integrating the
PK ODE simultaneously with the survival ODE.

This is a completely different simulation structure from the current event-driven
ODE engine. The "ferx already has ODE infrastructure" claim is misleading — the
ODE engine computes compartment states at observation times, not the cumulative
hazard integral. The survival integral requires either:
- A separate ODE for the cumulative hazard (dH/dt = h(t)), run alongside the PK
- Or Gaussian quadrature over the hazard at each observation interval

This is 3-4 months of ferx-core work, not 6-8 weeks.

### 3.3 [BLOCKER] GPU acceleration incompatible with Enzyme AD

Section 13a proposes GPU acceleration using CUDA or wgpu. Critical omission: the
current AD framework (Enzyme) operates on LLVM IR for CPU targets. Enzyme does
NOT support CUDA or GPU codegen. There is ongoing Enzyme research for GPU but it
is not production-ready.

This means GPU acceleration of the AD-differentiated likelihood is NOT achievable
by extending the current Enzyme approach. To run on GPU, one of:
1. Implement hand-coded forward-mode AD for GPU kernels (enormous effort)
2. Use a different AD framework that supports GPU (e.g., ad-hoc WGSL, cuDNN custom ops)
3. Run the likelihood on GPU without AD (use finite differences for the GPU path)
4. Use GPU only for ODE solving (RK45 on GPU) and keep the gradient on CPU

Option 3 (FD on GPU) is the most pragmatic but gives up the key advantage (exact
AD). Option 4 (GPU ODE + CPU gradient) is architecturally complex. The plan
presents GPU as a straightforward parallelism extension; it is actually an
architectural constraint that requires a fundamental decision about the AD
framework.

**Fix:** Acknowledge the Enzyme-GPU incompatibility explicitly. Reframe GPU
acceleration as "GPU-accelerated ODE solving with CPU-side AD" and estimate
6-12 months rather than 3-6 months.

### 3.4 [BLOCKER] Regulatory validation not mentioned anywhere

Three plan documents totalling 1693 lines make zero mention of the regulatory
validation required for pharmacometric software in regulatory submissions.

For ferx to be used in NDA/BLA/MAA submissions:
- The software must be validated per ICH E9, FDA 21 CFR Part 11, and EMA
  annex 11 requirements
- A formal Validation Plan and Validation Report are required
- Each new estimation method or model type requires separate validation
- Validation typically requires: theoretical justification, numerical comparison
  against a reference (NONMEM or literature), test cases with known answers,
  edge case analysis

This is not a minor administrative step. It is a multi-month process per feature
and must be planned for from the start. A feature that works correctly but has no
validation documentation cannot be used in a regulatory submission. Industry
adoption will be blocked until at least a preliminary validation package exists.

**Fix:** Add a dedicated section on regulatory validation strategy. Every feature
in the roadmap must have a corresponding validation approach. The concordance
tests are the beginning of a validation strategy but far from sufficient for
regulatory use.

### 3.5 [MAJOR] Full Bayesian HMC — the dimension problem

Section 3 presents Option B (SAEM + HMC polish) as "6-8 weeks" and the natural
path. The fundamental problem:

HMC on the joint posterior `p(theta, {eta_i} | data)` for N subjects with d ETAs
each has dimension `n_theta + N * d`. For a 2-cpt oral model with 5 thetas, 4
ETAs, and 200 subjects: dimension = 5 + 800 = 805. Standard HMC (identity mass
matrix) is known to mix poorly in >100 dimensions — each leapfrog step must make
a half-step, full-step, half-step for 805 parameters.

The correct approach for high-dimensional NLME HMC is the **partially collapsed
Gibbs** strategy or **marginalisation**: integrate out the ETAs analytically (or
via Laplace approximation), leaving only the n_theta-dimensional theta posterior.
This is what Stan does for its hierarchical models (`integrate_ode_rk45` + Laplace
marginalisation). It requires computing the Laplace approximation of the
individual posteriors — which FOCE already does — and using that as a surrogate
for the marginal.

The plan doesn't address dimensionality. Option B ("polish the SAEM mode with
HMC") doesn't specify what space is being sampled or how the 800-dimensional
individual effect space is handled.

**Fix:** Specify the marginalisation strategy explicitly. The 6-8 week estimate
is only valid if HMC operates in theta-space only (marginalising over eta via the
Laplace approximation). Sampling the full joint is a research project.

### 3.6 [MAJOR] OED / FIM — constraint structure for design optimisation ignored

Section 4 proposes BOBYQA for optimising sampling times. BOBYQA handles box
constraints `lower[i] <= x[i] <= upper[i]` for each dimension independently.
Sampling time optimisation has structured constraints:

- Times must be strictly increasing within a subject: `0 <= t_1 < t_2 < ... < t_n`
- Times must be within the study window: `0 <= t_i <= T_max`
- Subjects at the same dose group must have the same nominal sampling schedule

These are linear inequality constraints, not simple box constraints. BOBYQA
cannot enforce `t_1 < t_2 < t_3` as written. The plan is wrong about the
optimizer choice for this specific subproblem.

**Fix:** Either use SLSQP (which handles linear inequality constraints) with
explicit monotonicity constraints, or reparameterise the sampling times as
successive differences `delta_i = t_{i+1} - t_i > 0` (converting ordering
constraints to box constraints), then optimise in the delta space.

### 3.7 [MAJOR] MIPD — a medical device, not an R function

Section 5 presents `ferx_mipd()` as a near-term R function. The regulatory
reality:

Any software that makes or supports clinical dosing decisions in the US falls
under FDA Software as Medical Device (SaMD) guidance. Depending on the risk class:
- Class I: general wellness, low risk
- Class II: moderate risk (510(k) required, e.g., TDM software providing dosing
  suggestions)
- Class III: high risk (PMA required, e.g., closed-loop systems)

An R function that "returns the dose that achieves the target" for clinical use
is likely Class II SaMD. This requires:
- A 510(k) premarket notification or De Novo classification
- Quality Management System (21 CFR Part 820)
- Software development lifecycle documentation (IEC 62304)
- Clinical performance evaluation

This is not a software engineering problem; it is a regulatory strategy problem.
The plan presents it as "4-6 weeks" which is the coding time only. The path to
clinical use is 12-18 months of regulatory process, assuming the regulatory
classification is correct.

**Fix:** Separate `ferx_mipd()` as a research/development tool (acceptable
without SaMD classification if not used for clinical decisions directly) from
a clinical-grade MIPD system (requires full regulatory clearance). The R package
can be the research tool; the clinical product is a different scope.

### 3.8 [MAJOR] Structural identifiability claim is backwards

Section 6b: "algebraic structural identifiability is computationally expensive
for complex ODE models but tractable for 1-3 cpt analytical models."

This is backwards. Analytical 1-3 cpt models have known identifiability
properties (all are globally identifiable from rich data — this has been proven
algebraically decades ago and is in pharmacokinetics textbooks). The hard, open
problem is structural identifiability for **complex ODE models** where the
analytical proof is intractable.

DAISY and similar tools were developed specifically for ODE systems where
observability and identifiability must be checked symbolically. For analytical
PK models, you don't need DAISY — you need a lookup table.

**Fix:** Correct the description. Section 6b should focus on ODE models (where
the capability is novel and useful) and note that analytical PK model
identifiability is already known.

### 3.9 [MAJOR] NPML complexity underestimated — not just an E-step change

Section 7 presents NPML as "replace the Gaussian prior with a discrete mixture
during the E-step" in 6-10 weeks. The full scope:

1. The NPML support points live in eta-space. Their positions must be optimised
   (gradient-based or EM). This is a nested optimisation problem.
2. The number of support points K must be chosen (AIC/BIC on the NPML marginal
   likelihood, or the 2-step NPEM approach).
3. Starting values for support points are critical — bad starts produce degenerate
   solutions where all mass concentrates on one point.
4. The standard algorithm (Mallet 1986) alternates between estimating the masses
   (linear programming problem) and the support points (gradient ascent). The LP
   step requires a linear programming solver (ferx currently has no LP dependency).
5. Simulation from the estimated NPML distribution (for VPC) requires sampling
   from a discrete distribution over the support points — the standard
   `ferx_simulate()` infrastructure doesn't support this.
6. Diagnostics: standard CWRES/IWRES are valid; ETA shrinkage is not defined for
   NPML (there are no ETAs, only posterior class probabilities).

This is closer to 3-6 months than 6-10 weeks.

### 3.10 [MAJOR] Mixture models — parameter space constraints ignored

Section 8 proposes mixture model syntax:
```
mixture(K=2, probs=[PI1, PI2])
theta TVCL_class1(0.1, ...)
theta TVCL_class2(0.5, ...)
```

The constraints:
- `PI1 + PI2 = 1` (probabilities sum to 1) — one parameter is redundant;
  only K-1 mixing probabilities are free
- `0 < PI1 < 1`, `0 < PI2 < 1` — box constraints needed
- `TVCL_class1 != TVCL_class2` — mixture is non-identifiable if the class
  means are swapped (label switching problem)

The label switching problem is fundamental: if you start the optimiser with
class1 and class2 swapped, the EM converges to the same solution with labels
exchanged. Point estimates are therefore undefined without an identification
constraint (e.g., `TVCL_class1 < TVCL_class2`).

NONMEM handles this partly via the user specifying ordered initial values; partly
it doesn't handle it well (NONMEM $MIX has known convergence issues). The plan
doesn't acknowledge label switching at all.

**Fix:** The mixture model design must include:
1. K-1 free probabilities (parameterised as softmax to enforce sum-to-one)
2. Optional ordering constraints on class-specific thetas
3. Warning when classes are near-identical (potential unidentifiability)
4. Multiple random starts (label switching makes single-start convergence
   unreliable)

### 3.11 [MAJOR] Joint PK-TTE models — shared gradient path not addressed

Section 9: "the joint likelihood is `log L_joint = log L_PK + log L_TTE`, both
conditioned on the same individual parameters. FOCE handles this naturally."

"Naturally" is too optimistic. The challenge:

The FOCE E-step minimises the individual objective for each subject. For PK only,
this is the Laplace approximation of `log p(y_PK | eta, theta)` with respect to
eta. For joint PK-TTE, it becomes `log p(y_PK | eta, theta) + log p(TTE | eta,
theta)` — both terms contribute to the eta gradient.

The gradient of `log p(TTE | eta, theta)` with respect to eta requires
differentiating through the survival integral, which itself depends on the PK ODE
trajectory (which depends on eta). This is a second-order AD problem: the hazard
function depends on concentration, which comes from the ODE, which is
differentiated with respect to eta. The current Enzyme AD setup would need to
differentiate through the ODE solver at both the eta and theta levels
simultaneously.

This is not handled by simply adding a second likelihood term. It requires careful
implementation of the AD chain through the combined ODE system.

### 3.12 [MINOR] "No other NLME tool has native MIPD support" — incorrect

Section 5: "No other NLME estimation engine has native MIPD support."

This is not accurate. MwPharm++ (Mediware/Morf Medics), InsightRx, DoseMeRx,
and Tucuxi are all commercial MIPD platforms in clinical use. They have their own
MAP estimation engines, though these are based on older algorithms. Additionally,
the R package `PKPDmap` and the `BayesianPK` Python package provide MAP Bayesian
dosing. The claim should be: "no open-source NLME estimation engine has native
MIPD support as a first-class feature."

### 3.13 [MINOR] Stan comparison claim "better suited than any current NLME tool"

Section 3: "ferx's exact AD makes it better suited for this than any current
NLME tool." Stan also uses exact AD (via templated C++ with its own autodiff
library). The claim that ferx is better suited than Stan for Bayesian NLME needs
a specific argument (Rust vs C++ performance? NLME-specific priors? Integration
with PK structures?). As written it is unsupported.

### 3.14 [MINOR] Effort estimates are systematically 2-3x optimistic

Across all sections, the effort estimates follow a pattern of optimism:

| Feature claimed | Stated effort | More realistic |
|---|---|---|
| Poisson/NB likelihood | 3-4 weeks | 6-10 weeks (data format + diagnostics) |
| Ordinal logistic | 4-6 weeks | 3-4 months (threshold constraints, diagnostics) |
| VPC | 2-3 weeks | 6-8 weeks (binning, IOV, BLOQ, multi-DVID) |
| NPDE | 1-2 weeks | 3-5 weeks (decorrelation step) |
| TTE parametric | 6-8 weeks | 3-4 months (data format, survival ODE) |
| NUTS full Bayesian | 3-6 months | 6-12 months |
| NPML | 6-10 weeks | 3-6 months |
| GPU acceleration | 3-6 months | Not feasible with Enzyme; rearchitecture needed |

The optimism is systematic, not random. Every estimate is for the "happy path
where the existing infrastructure just works." The realistic estimate must include
data format changes, new diagnostics, test coverage, and edge cases.

---

## Cross-cutting gaps — missing from all three plans

### X.1 [BLOCKER] No mention of where features live (ferxtranslate vs ferx-r vs ferx-core)

The plans mix features from three separate repositories without specifying which
repo owns each:

- `ferxtranslate`: translation layer (R), owned here
- `ferx-r`: R interface to ferx (R + Rust via extendr), different repo
- `ferx-core`: Rust engine, different repo

`ferx_vpc()`, `ferx_compare()`, `ferx_mipd()`, `ferx_covariate_screen()`,
`ferx_bootstrap()`, `ferx_influence()` — these would logically live in ferx-r,
not ferxtranslate. The ferx_ir modifier API (`ir_add_covariate()`, etc.) lives
in ferxtranslate. The new likelihoods (Poisson, ordinal, TTE) live in ferx-core.

Building something in the wrong repo creates maintenance, dependency, and API
coherence problems later. This decision must be made before Phase 1 begins.

### X.2 [MAJOR] No testing strategy for new statistical methods

Every new statistical method (Poisson likelihood, TTE, Bayesian, NPML) requires
a simulation study to verify it works correctly:
- Generate data from the assumed model
- Fit with the new method
- Verify parameter recovery, coverage of CI, type I error of LRT

These validation simulation studies are not mentioned anywhere. They are essential
before claiming any method is correct. For regulatory acceptance, they must be
documented and published (or at minimum available as reproducible vignettes).

The concordance tests (test-concordance.R) are the right idea but they test only
a handful of models at fixed true values. A proper validation study sweeps across
parameter values, sample sizes, and model complexity.

### X.3 [RESOLVED] .ferx format versioning — cargo feature flags, not separate crates

The concern was: adding TTE, count, mixture syntax breaks existing `.ferx` files
and old parsers. The resolution is already in the codebase — cargo feature flags.

`ferx-core/Cargo.toml` already has:
```toml
[features]
default  = ["autodiff"]
nn       = []     # neural covariate model — already shipped this way
```

The correct extension:
```toml
tte      = []     # time-to-event likelihoods
count-pd = []     # Poisson, NB, ordinal, binary
mixture  = []     # finite mixture / latent class
full     = ["nn", "tte", "count-pd", "mixture"]
```

**A separate `ferx-core-tte` crate is the wrong level.** Joint PK-TTE models
need the full PK engine, so `ferx-core-tte` would be a superset (effectively a
fork). `ferx-r` can only link one native binary. Feature flags solve this
natively in Rust.

**The format versioning concern is smaller than it looked.** Adding new syntax
is additive: old parsers see an unknown token and error cleanly; new parsers read
old files fine. An explicit format version number is only needed if the *meaning*
of existing syntax changes — nothing in the roadmap proposes that. The parser
must produce a readable error ("unknown error model type 'tte': rebuild ferx with
--features tte") rather than a cryptic crash, which is a parser robustness issue,
not a versioning issue.

### X.4 [MAJOR] GOF plots -- the most basic diagnostic, completely missing

The plans discuss VPC, NPDE, PSIS-LOO, bootstrap CIs -- all relatively advanced
diagnostics. The most basic goodness-of-fit plots used in every pharmacometric
report are not mentioned:

- PRED vs DV (population fit)
- IPRED vs DV (individual fit)
- CWRES vs TIME (residual trend)
- CWRES vs PRED (heteroscedasticity)
- ETA vs covariate (covariate discovery)
- Histogram and Q-Q plot of ETA distribution

These are 1-2 weeks of ggplot2 code in ferx-r. `ferx_gof(fit, data)` returning a
patchwork of these plots would be the single highest-frequency use case after
fitting. More frequently used than VPC.

**Update (2026-06-02):** `ferx_gof()` is now described in detail in
ferx-output-advantage.md section 3g. Prerequisites (sdtab ID, CMT column) are
fixed in ferx-core PR#172. The "live questions" about CMT and eta_normality
ordering are answered. ferx_gof() is the recommended next build target.

### X.5 [MAJOR] Run management / project organisation not considered

In industry, a modelling project produces 20-200 model runs. Currently:
- NONMEM uses NMproject or manual directory naming (run001, run002, etc.)
- There is no equivalent in the ferx ecosystem

A `ferx_project` object that tracks runs, stores fit results, and generates a
run record (like PsN runrecord) would be the workflow infrastructure that makes
ferx usable for multi-run projects. Without it, users manage this manually.

This is not in any of the three plans but is a prerequisite for the covariate
screening and influence diagnostic tools (which themselves generate many runs).

### X.6 [MINOR] Population generation for simulation — missing prerequisite for VPC

`ferx_vpc()` requires simulating virtual subjects. For VPC of a dataset, you
simulate from the same design (same covariate values, same dosing) as the
observed data. This is straightforward.

For clinical trial simulation (not yet in the plans), you need to generate a
synthetic population with realistic covariate distributions (WT, AGE, CRCL, etc.).
This requires a population generator — not mentioned anywhere. The `ferx_vpc()`
function implicitly solves the simpler case; the harder case is completely absent.

### X.7 [MINOR] Covariate forward plots and forest plots

These are among the most-requested visualisations in PMx:

- **Forest plot**: shows parameter estimates ± CI across covariate subgroups
  (WT quartiles, age groups, renal function). Standard in every regulatory submission.
- **ETA vs covariate plot**: shows correlation between individual ETAs and
  covariates; the primary tool for identifying missing covariates.

`ferx_r` already has `eta_corr_test` (Shapiro-Wilk per ETA). The ETA vs
covariate scatter plot is a trivial 1-week addition. The forest plot requires
simulation from the fitted model at extreme covariate values — 2-3 weeks. Neither
is in any plan.

### X.8 [MINOR] Drug-drug interaction (DDI) model templates

DDI models (competitive inhibition, mechanism-based inactivation, enzyme
induction) follow standard structural forms. While ferx's ODE engine handles
them, having standard template models (`ferx_ddi_model(type="competitive_inhibition")`)
would make ferx the go-to tool for DDI assessment (a large regulatory use case).
Not in any plan.

---

## Summary: the top five blockers to address before implementing anything

1. **Regulatory validation strategy** — without it, industry adoption is blocked
   regardless of how good the software is. Must be designed first.

2. **Data format extensibility design** — count, ordinal, TTE, and mixture data
   all require format changes. Design the extensible format once, not piecemeal
   per feature.

3. **Repository assignment** — which features live in which repo. Cannot be
   retrofitted.

4. **ferx_ir rhs storage decision** — string vs R expression object. Every
   modifier function depends on this. Cannot be changed after the modifier API is built.

5. **Target audience commitment** — without it, phase ordering is undefined and
   the roadmap is unactionable.

## Summary: the top five realistic effort corrections

| Item | Plan says | Reality |
|---|---|---|
| VPC | 2-3 weeks | 6-8 weeks for an MVP; 3 months for production quality |
| Poisson/NB likelihoods | 3-4 weeks | 6-10 weeks (includes data format + diagnostics) |
| TTE parametric | 6-8 weeks | 3-4 months |
| NUTS full Bayesian | 3-6 months | 6-12 months |
| GPU acceleration | 3-6 months | Requires rearchitecting AD; timeline undefined |
