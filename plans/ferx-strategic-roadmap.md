# ferx strategic roadmap: becoming the definitive NLME engine

This document evaluates what additions — statistical, computational, workflow,
and ecosystem — would make ferx the default choice for pharmacometric modelling.
It covers three horizons: near-term (months), medium-term (6-18 months), and
long-term (research territory). Every item is evaluated for effort, genuine
novelty, and strategic leverage.

The audience is clinical pharmacology and pharmacometrics. The current
competition is NONMEM (industry standard, proprietary, slow, text-centric),
Monolix (good UX, no AD, estimation accuracy comparable to NONMEM), nlmixr2
(open, flexible, uses R ecosystem), and Stan (gold-standard Bayesian, not
pharmacometrics-specific).

---

## What ferx already has (baseline)

**Estimation:** FOCE, FOCEI, Gauss-Newton (BHHH), GN-Hybrid, SAEM,
Importance Sampling (IMP), SIR. HMC exists in the codebase as an E-step
proposal inside SAEM (leapfrog integrator is already implemented).

**Error models:** Additive, Proportional, Combined — all Gaussian. LTBS
(log-transform-both-sides). M3 (Tobit likelihood for BLOQ). Per-CMT dispatch.

**Structural:** 1/2/3-cpt analytical PK (oral, IV, infusion). Arbitrary ODE.
SDE with Extended Kalman Filter likelihood. Neural network covariate model
(DCM, [covariate_nn] block).

**Random effects:** Diagonal and block omega. IOV via kappa. Lognormal and
additive parameterisations. Logit-normal for bioavailability F.

**Diagnostics:** AIC, BIC, condition number, eigenvalues, ETA Shapiro-Wilk
normality test, Durbin-Watson IWRES autocorrelation, shrinkage, EBE OFV.

**Output:** Full covariance matrix, SIR CIs, IMP marginal log-likelihood,
individual estimates, sdtab with CWRES/IWRES. `.fitrx` bundle format.

**Key infrastructure:** Exact automatic differentiation (Enzyme). Rayon CPU
parallelism. Multi-start. NCA-based initial estimates.

---

## 1. Non-Gaussian observation likelihoods

### Delivery model: cargo feature flags, not separate crates

New likelihood types are gated behind cargo feature flags — the same pattern
already used for `nn` (deep compartment model). `ferx-core/Cargo.toml` gains:

```toml
count-pd = []   # Poisson, NB, ordinal, binary
```

`ferx-r` builds ferx-core with `--features count-pd` for the relevant release.
Old `.ferx` files are unaffected — the parser is additive. A binary built without
`count-pd` seeing `DV ~ poisson(...)` produces a clean "unknown error model type"
error, not a crash.

### Current state
The `ErrorModel` enum has three variants: `Additive`, `Proportional`,
`Combined`. All assume the residual is Gaussian.

### What to add

**1a. Count data: Poisson and negative binomial**

Pharmacodynamic endpoints measured as counts — seizure frequency, pain event
counts, hospitalisation counts — require integer-valued likelihoods. NONMEM
handles these via custom PRED blocks; no standard count-data support exists in
any NLME engine as a first-class feature.

Ferx error model syntax could be:

```
[error_model]
  DV ~ poisson(lambda=IPRED)
  DV ~ negative_binomial(mu=IPRED, phi=DISP)
```

Where `DISP` (overdispersion) is estimated alongside structural parameters.
The likelihood is exact: `log P(y | lambda)` with no Gaussian approximation.

AD compatibility: Poisson log-likelihood is a sum of log-factorials (handled
via lgamma) — fully AD-differentiable. NB likewise.

**1b. Ordered categorical: proportional odds**

Sedation scales, pain scores, adverse event grades, and many PD endpoints are
ordinal. The proportional odds model assigns a latent continuous score and
cuts it at thresholds:

```
[error_model]
  DV ~ ordered_logistic(eta=SCORE, thresholds=[TH1, TH2, TH3])
```

Where `SCORE = IPRED` is the individual predicted latent score and `TH1..THK`
are estimated threshold parameters.

NONMEM implements this via custom $ERROR code (FORTRAN). No NLME engine has
it as a first-class syntax. This is a genuine differentiator.

**1c. Binary / Bernoulli**

Binary PD endpoints (response/no-response, toxicity yes/no). Logistic mixed
model with IIV on the linear predictor:

```
[error_model]
  DV ~ bernoulli(p=inv_logit(IPRED))
```

**1d. Zero-inflated variants**

Zero-inflated Poisson or negative binomial for overdispersed count data with
excess zeros. Requires a mixture weight parameter (probability of structural
zero). Common in adverse event counts.

**Effort:** 1a (Poisson/NB) is 3-4 weeks in ferx-core. 1b (ordinal) is 4-6
weeks. 1c is 1-2 weeks. 1d follows from 1a+1b, 2-3 additional weeks.

**Strategic value:** Immediately opens ferx to the entire PD modelling space
that NONMEM handles only via laborious custom code. Any pharmacometrician
fitting a count or categorical endpoint will choose the tool that has native
support.

---

## 2. Time-to-event (TTE) and survival models

Feature flag: `tte = []` in ferx-core. Joint PK-TTE models need the full PK
engine — a separate `ferx-core-tte` crate would be a fork. Feature flags keep
it in one binary, opt-in at build time.

### Why this matters
Oncology pharmacometrics revolves around overall survival, progression-free
survival, and time to response. Joint PK-TTE models (where PK drives the
hazard function) are the state of the art. NONMEM handles TTE via `$DES`
(numerically integrating the survival function as an ODE), which works but is
cumbersome. No NLME engine has TTE as a first-class citizen.

### What to add

**2a. Parametric TTE**

Parametric hazard functions as native error model types:

```
[error_model]
  DV ~ tte(hazard=BASE_HAZARD * exp(-BETA * CONC), dist=weibull)
  DV ~ tte(hazard=BASE_HAZARD * (1 + EMAX * CONC / (EC50 + CONC)), dist=exponential)
```

Ferx already has the ODE infrastructure (survival function S(t) = exp(-integral
of h(t) dt) can be computed via the ODE solver). The innovation is wrapping
this in a clean syntax with automatic event-time log-likelihood contribution.

Supported baseline distributions: exponential, Weibull, Gompertz, log-normal,
log-logistic.

**2b. Repeated TTE (recurrent events)**

When the endpoint repeats (e.g., seizures, hospitalisations), the likelihood
becomes a product of inter-event time contributions. The Andersen-Gill model
for recurrent events is the standard approach.

**2c. Joint PK-TTE models**

The PK ODE runs simultaneously with the survival ODE. Individual PK parameters
(from IIV) drive the hazard. The joint likelihood is:

```
L_total = L_PK * L_TTE = prod_i [ p(conc_i | theta_i) * p(TTE_i | theta_i) ]
```

Both contributions use the same theta/eta. This is the current gold standard
for benefit-risk modelling but requires manual implementation in NONMEM.

**Effort:** 2a (parametric TTE) is 6-8 weeks in ferx-core. 2b is additional
4-6 weeks. 2c (joint) is additional 4-6 weeks after 2a is working.

**Strategic value:** The highest-value addition for the oncology market. Every
oncology PK/PD submission that includes survival modelling would prefer ferx
if TTE is native.

---

## 3. Full Bayesian estimation (HMC/NUTS)

### Why this matters
The HMC leapfrog integrator is already implemented in ferx-core (`hmc.rs`) as
an E-step proposal for SAEM. The infrastructure for exact AD gradients (Enzyme)
is also there. Stan's dominance in Bayesian statistics is built on exactly this
combination: HMC + AD. The gap is using HMC as a full population-level sampler
rather than a within-subject E-step tool.

### What full Bayesian would give

- Proper uncertainty quantification: posterior distributions over all
  parameters, not just point estimates with Wald SEs
- Naturally handles non-identifiability (posterior remains broad rather than
  crashing the covariance step)
- Informative priors for sparse data populations (pediatrics, rare disease)
- MCMC convergence diagnostics (R-hat, effective sample size per parameter)
- No covariance step that can fail — uncertainty is sampled, not approximated

### Architecture options

**Option A: Full HMC population sampler (No-U-Turn Sampler, NUTS)**

NUTS is HMC with adaptive step size and path length. It requires:
- A joint log-posterior function: `log p(theta, {eta_i}) | data`
- AD gradient of the above with respect to all parameters
- NUTS adaptation (dual averaging for step size, online mass matrix estimation)

The per-subject NLL (already computed by FOCE) becomes a likelihood term.
Priors go on theta (Normal, half-Normal, log-Normal, uniform).

The marginalised sampler (marginalises over eta analytically) would be the
most efficient but is complex. The joint sampler (samples theta AND eta
simultaneously) is straightforward.

**Option B: SAEM + HMC E-step → post-hoc HMC refinement**

Run SAEM to convergence for point estimates, then run a short HMC chain
initialised at the SAEM mode to characterise the posterior. This is much
cheaper than pure MCMC from scratch and produces calibrated posterior intervals.

Moreoever: SAEM already uses HMC proposals — so the extension is making the
post-convergence phase a full posterior sampler rather than stopping at the mode.

**Option C: Variational Bayes (ADVI)**

Approximate the posterior with a factored Gaussian (mean-field VI). Optimise
the ELBO (evidence lower bound) using AD gradients. Much faster than MCMC;
useful as a quick uncertainty estimate or as an initializer for HMC. Used in
Stan's `vb()` mode and PyMC.

**Recommended path:** Option B as the near-term deliverable (weeks), Option A
as the full implementation (months). Option C as a speed-accuracy tradeoff.

**Effort:** Option B (SAEM + HMC polish) 6-8 weeks. Option A (full NUTS) 3-6
months.

**Strategic value:** Very high. Bayesian NLME is the future of regulatory
modelling (FDA guidance increasingly accepts Bayesian approaches for pediatric
extrapolation, rare disease, and adaptive trials). ferx's exact AD makes it
better suited for this than any current NLME tool.

---

## 4. Optimal experimental design (OED / PopED)

### Why this matters
Before running a clinical trial, pharmacometricians compute the Fisher
Information Matrix (FIM) to evaluate how precisely parameters can be estimated
given a proposed design (sampling times, doses, number of subjects). This drives
protocol optimisation.

Currently: PopED (R package) and PFIM (R package) do this. Both use approximate
FIM computation. Neither uses exact AD.

### ferx's advantage

The expected FIM under a population model is:

```
FIM(design) = sum_i E_eta[ (d log L_i / d theta)^2 ]
```

This is exactly the Hessian of the expected log-likelihood with respect to
theta, averaged over subjects. With ferx's exact AD, this can be computed
analytically rather than via finite differences. The result is a more accurate
FIM and faster computation.

### What to build

`ferx_fim(model, design, n_sim=1000L)`:
- Takes a ferx_ir + a proposed design (times, doses, n_subjects, dose regimen)
- Computes the expected FIM via Monte Carlo (simulate ETAs, compute exact Hessian)
- Returns: expected SE per parameter, D-optimal criterion (det(FIM)^(1/p)),
  power to detect a covariate effect, and optimal design suggestions

`ferx_opt_design(model, design_space, criterion="D")`:
- Optimises the design (sampling times, doses) to maximise the chosen criterion
- Uses an outer optimiser (BOBYQA for sampling time optimisation — already in NLopt)

**Effort:** FIM computation 3-4 weeks. Design optimisation 4-6 additional weeks.

**Strategic value:** Very high. OED is standard practice before phase I trials
and for protocol amendments. Currently requires separate tools (PopED) with
manual workflow. Built-in OED would be unique among estimation engines.

---

## 5. MAP Bayesian estimation and MIPD (Model-Informed Precision Dosing)

### Why this matters
In clinical practice, a population PK model is used to individualise dosing.
A patient's drug concentration is measured at 1-3 timepoints; a MAP Bayesian
estimator uses the population model as a prior to estimate that patient's
individual PK parameters; the estimated parameters are used to compute an
optimal dose for the next interval.

This is called Model-Informed Precision Dosing (MIPD) and is the primary
clinical translation of population PK modelling. NONMEM can compute MAP
estimates but is not designed for clinical use. Dedicated TDM software
(MwPharm, InsightRx, DoseMeRx) does this but requires a pre-loaded model.

### What ferx has and what's missing

ferx EBE computation is already MAP estimation. The `individual_estimates`
output is each subject's MAP estimate. The missing pieces are:
- An API designed for clinical use (one subject, iterative updates as new
  observations arrive)
- Dose recommendation given an individual PK profile and a target (e.g.,
  AUC24h = 400 mg·h/L)
- Uncertainty propagation from MAP estimates to dose recommendation

### What to build

`ferx_map(model, fit, individual_data)`:
- Given a population fit and a set of individual observations (times + DVs),
  compute that individual's MAP ETA estimates
- Returns individual parameter estimates with posterior uncertainty

`ferx_mipd(model, fit, individual_data, target_metric, target_value, dosing_interval)`:
- Computes the MAP estimates
- Simulates the individual PK forward in time under candidate doses
- Returns the dose (with credible interval) that achieves the target

**Effort:** `ferx_map` 2-3 weeks (thin wrapper around existing EBE machinery).
`ferx_mipd` 4-6 weeks (requires forward simulation + dose search).

**Strategic value:** Very high. This is the bridge between drug development and
clinical practice. No other NLME estimation engine has native MIPD support.
ferx's speed (milliseconds per individual fit) makes it deployable in real-time
clinical systems.

---

## 6. Structural identifiability analysis

### The problem
A model can be statistically identifiable in theory (enough data to estimate
parameters) but structurally non-identifiable from a given dataset (e.g., a
3-cpt model with only central compartment observations). Fitting proceeds,
produces estimates, but they are meaningless. NONMEM does not warn about this.

### What to build

**6a. Practical identifiability: simulation-based recovery test**

`ferx_identifiability(ir, data_design, n_sim=200L)`:
1. Generate `n_sim` datasets by simulating from the model at the prior/initial
   parameter values
2. Fit each dataset
3. Report: coverage of 95% CI (should be ~95%), bias, and precision of each
   parameter
4. Flag parameters that cannot be recovered (high bias, huge variance)

This is related to the simulation-based power analysis workflow, but focused
on parameter recovery rather than hypothesis testing.

**6b. Structural identifiability: differential algebra (long-term)**

Algebraic structural identifiability analysis (like the DAISY algorithm) tests
whether parameters are globally/locally identifiable from the model equations
alone, before seeing any data. This is computationally expensive for complex
ODE models but tractable for 1-3 cpt analytical models.

**Effort:** 6a is 3-4 weeks. 6b is a research project (months).

**Strategic value:** High. Catching non-identifiability before fitting saves
enormous time in drug development. A warning like "parameter TVQ is likely not
identifiable from your current study design — consider adding peripheral samples"
at the fit stage would be novel and valuable.

---

## 7. Nonparametric and semiparametric methods

### The case for nonparametric

Standard NLME assumes ETA ~ N(0, Omega). This assumption fails when:
- A population has a mixture of metabolizer phenotypes (bimodal ETA)
- Covariate effects are misspecified (the "ETA" absorbs non-Gaussian structure)
- The sample is enriched (e.g., patients who responded vs didn't)

### 7a. Nonparametric maximum likelihood (NPML / Mallet)

Estimate the discrete empirical distribution of ETAs (support points + masses)
rather than assuming Gaussianity. Each subject is assigned to one of K support
points. K is chosen by parsimony criteria.

This is what NONMEM NONPARAMETRIC does (via the NPEM algorithm). Also
available in nlmixr2 via the npde approach.

Implementation in ferx: replace the Gaussian prior on eta with a discrete
mixture during the E-step. Support points are optimised; masses via EM. The
FOCE likelihood structure (subject-level E-step) is compatible.

**7b. Kernel density ETA distribution**

Replace the Gaussian prior with a kernel density estimate updated iteratively
(like SAEM but nonparametric). Less principled than 7a but easier to implement.

**Effort:** 7a is 6-10 weeks (significant change to the E-step). 7b is 3-4
weeks.

**Strategic value:** Medium-high. NPML is particularly valuable for paediatric
modelling (mixed metaboliser populations) and oncology (responder/non-responder
mixtures). Would close the NONMEM NONPARAMETRIC gap.

---

## 8. Mixture models (latent class)

### Current state
`[parameters]` has no `mixture` keyword. The ferx-core parser does not support
it. This is the largest single capability gap vs NONMEM.

### What it enables
A mixture model assigns each subject probabilistically to one of K classes,
each with different typical values:

```
[parameters]
  mixture(K=2, probs=[PI1, PI2])
  theta TVCL_class1(0.1, ...)
  theta TVCL_class2(0.5, ...)
  ...
[individual_parameters]
  CL = mix(TVCL_class1, TVCL_class2) * exp(ETA_CL)
```

The EM structure of FOCE/SAEM is natural for mixtures: the E-step computes
posterior class probabilities for each subject; the M-step updates class-
specific parameters.

**Applications:** poor/extensive metaboliser distinction, responder/non-responder,
cancer responders.

**Effort:** 8-12 weeks (parser changes + E-step changes + output).

**Strategic value:** High. Currently listed as a ferx-core gap. NONMEM supports
this ($MIX). Not having it is a concrete deficiency for any mixture population
analysis.

---

## 9. Joint longitudinal and survival models

### Why this matters
The most sophisticated model in pharmacometrics integrates:
- A continuous longitudinal PK/PD trajectory (the "longitudinal" part)
- A time-to-event endpoint influenced by the trajectory (the "survival" part)

Example: tumour size over time (longitudinal) + overall survival as a function
of rate of tumour shrinkage (survival). This is the Ribba/Bruno joint model
framework used in oncology regulatory submissions.

### Current state
Neither NONMEM nor any other NLME engine supports this natively. NONMEM
approximates it via a joint ODE (survival ODE + PK ODE in $DES) and combined
likelihood tables. The JM R package does this for simpler cases (lme4 +
survial), not for NLME PK/PD.

### What to build
A joint model combines:
- The existing FOCE likelihood for the longitudinal component
- The TTE likelihood (Section 2) for the survival component
- Shared random effects (same ETA drives both)

The joint likelihood is:
```
log L_joint = log L_longitudinal + log L_TTE
```

Both are conditioned on the same individual parameters. FOCE handles this
naturally — the E-step objective includes both terms.

**Effort:** After TTE is implemented (Section 2), the joint model is 4-6
additional weeks.

**Strategic value:** Very high for oncology. Would be a first-in-class
capability among dedicated NLME engines.

---

## 10. Optimizers — what's actually missing

### Current optimizer landscape in ferx
SLSQP (gradient), COBYLA (derivative-free), MMA (gradient), BFGS (gradient),
trust-region Newton (gradient), BOBYQA (derivative-free). NLopt underlies most
of these. The outer optimizer choices are already good.

### What would genuinely help

**10a. L-BFGS-B** (limited-memory BFGS with box constraints)
- More memory-efficient than BFGS for large parameter vectors
- The standard in ML optimization (used in PyTorch, scipy.optimize)
- Would help for models with many thetas (PBPK, population covariate models)
- Effort: 1-2 weeks if done via NLopt or a direct implementation

**10b. CMA-ES** (Covariance Matrix Adaptation Evolution Strategy)
- Derivative-free global optimizer; state-of-the-art for moderately
  high-dimensional problems
- Superior to BOBYQA for multimodal likelihood surfaces
- Used in some pharmacometric global search applications
- Effort: 2-3 weeks (CMA-ES has clean open-source implementations to wrap)
- Strategic value: Medium. Multi-start covers most use cases but CMA-ES is
  more systematic and adapts its search distribution.

**10c. Stochastic gradient (Adam) for SAEM**
- The SAEM stochastic approximation step could use Adam-style adaptive
  learning rates instead of the Robbins-Monro schedule
- Reduces sensitivity to step-size tuning in SAEM
- Used in deep learning SAEM variants (saemix with Adam)
- Effort: 2-3 weeks
- Strategic value: Low-medium. Would improve SAEM stability in difficult models.

**10d. Natural gradient for population-level step**
- Precondition the outer gradient with the Fisher information matrix
- Equivalent to Newton's method in the expectation-parameter space
- Known to improve SAEM convergence in practice
- Closely related to the existing Gauss-Newton method
- Effort: 3-4 weeks
- Strategic value: Medium.

---

## 11. Diagnostics borrowed from other fields

### 11a. PSIS-LOO (Pareto-smoothed importance sampling leave-one-out CV)

From Bayesian statistics (Vehtari, Gelman, Gabry 2017). A computationally
efficient approximation to leave-one-subject-out cross-validation:
- Uses importance weights computed from the existing posterior samples
- Each subject's LOO contribution is estimated without re-fitting
- Returns an LOO information criterion comparable across models
- The Pareto k̂ diagnostic flags subjects that are "surprising" given the model

For frequentist ferx (FOCE), the IMP stage already produces importance weights
per subject. PSIS-LOO is a natural extension.

**Effort:** 2-3 weeks.
**Strategic value:** High. LOO-CV is the theoretically preferred model
comparison criterion (better than AIC/BIC in finite samples). The Pareto k̂
per subject is a novel influential-subject diagnostic.

### 11b. Conformal prediction intervals for model predictions

Conformal prediction (Venn-Abers, split conformal) produces prediction
intervals with guaranteed marginal coverage without assuming the correct model.
Applied to ferx's `sdtab` predictions: for a new subject, what is the
guaranteed-coverage interval for their next observation?

This is distinct from VPC (which tests the model) — conformal PI is used for
clinical dosing decisions.

**Effort:** 2-3 weeks.
**Strategic value:** Medium. Increasingly asked for in regulatory submissions
as model-based dosing becomes more common.

### 11c. Posterior predictive p-values (Bayesian goodness-of-fit)

For each test statistic (max CWRES, % observations > 2 SD, etc.), compute the
fraction of simulated datasets where the test statistic exceeds the observed
value. This is the Bayesian analogue of a formal hypothesis test for model
misspecification.

Used in Stan via `posterior_predictive_check()`. Not available in any NLME
engine. Requires ferx_simulate (already exists).

**Effort:** 1-2 weeks once VPC simulation infrastructure is in place.
**Strategic value:** Medium. Gives a formal p-value for model adequacy,
complementing the visual VPC.

### 11d. Generalised additive models (GAMs) for covariate relationships

Instead of pre-specifying the covariate functional form (power, linear, etc.),
fit a smooth spline and let the data determine the shape. Useful for
discovering non-linear covariate effects before deciding on a parameterisation.

Implementation: after fitting the base model, extract ETA-BLUPs and residuals,
then fit mgcv::gam() to each ETA ~ covariate. Return the fitted smooth and the
predicted optimal functional form.

This is a post-hoc R-level tool (not engine changes): `ferx_covariate_gam(fit, data)`.

**Effort:** 2-3 weeks (pure R, uses mgcv).
**Strategic value:** High. Removes the "what functional form?" guesswork from
covariate model building. Would be unique among NLME tools.

### 11e. Functional data analysis of PK profiles

Treat each subject's full concentration-time profile as a functional object.
FPCA (functional principal component analysis) decomposes between-subject
variation into orthogonal "modes" — the PK analogue of the omega matrix but
without assuming lognormal distributions.

Useful for: identifying dominant modes of PK variability, detecting subjects
whose profiles are outliers, guiding random effect parameterisation.

Pure R implementation using the `fda` package on top of individual predictions
from ferx.

**Effort:** 2-3 weeks (R-level tool).
**Strategic value:** Medium. Research value; would be the first PK-specific FDA.

---

## 12. The PBPK direction (long-term)

Physiologically-based PK (PBPK) represents drug distribution as a network of
anatomical compartments with literature-fixed volumes and blood flows. Models
are typically pre-parameterised from in vitro data; clinical data is used to
adjust a small number of parameters (typically absorption rate or tissue binding).

**What would be needed:**
- A model library of PBPK structures (human adult, paediatric, pregnant, renal
  impaired) as ferx ODE definitions
- Prior distributions on PBPK parameters (the in vitro-to-in vivo extrapolation
  uncertainty)
- Bayesian estimation to update from sparse clinical data

ferx's ODE engine already handles arbitrary ODEs. The Bayesian infrastructure
(Section 3) is the prerequisite. The model library would be the main effort.

**Effort:** 12-24 months, requires domain expertise for model library.
**Strategic value:** High (long-term). PBPK is mandatory for paediatric
regulatory submissions. Currently owned by SIMCYP and PK-Sim (proprietary,
expensive). Open ferx-based PBPK would be a major disruptive move.

---

## 13. Speed and compute infrastructure

### 13a. GPU acceleration

Each subject's likelihood is computed independently. For SAEM (which needs many
Monte Carlo draws per subject per iteration), GPU parallelism is the natural
extension of ferx's Rayon CPU parallelism.

For ODE models with large subject counts (N > 500) and many parameters,
GPU-parallelised RK45 solvers would provide a 10-50x speedup over CPU-parallel
at the same thread count.

Framework options: CUDA via cuBLAS/cuDV for matrix operations, wgpu for
cross-platform GPU (Metal, Vulkan, CUDA).

**Effort:** 3-6 months (significant infrastructure).
**Strategic value:** High for large N, ODE models, SAEM. Would enable
analysis of electronic health record datasets (N = 10,000+) that are currently
intractable.

### 13b. Distributed fitting

Subject-level likelihood decomposition (FOCE subjects are independent given
theta) enables distributed computation: each node computes a subset of
subjects. Only the gradient sum needs to be communicated.

Standard pattern in distributed ML (parameter server). Would enable very large
datasets without GPU.

**Effort:** 2-3 months.
**Strategic value:** Medium. Most PMx datasets are small enough for CPU.
Real-world evidence (RWE) datasets are the use case.

### 13c. Python bindings (PyO3)

ferx-core is Rust. Exposing the core fit + simulate functions to Python via
PyO3 would enable:
- Use from Python pharmacometrics workflows
- Integration with PyTorch / JAX for hybrid ML-NLME models
- Scripting from Jupyter notebooks

**Effort:** 4-6 weeks for the core API.
**Strategic value:** High. The pharmacometrics Python ecosystem is growing (POPED-py,
nlme-ode). Being available in Python without a full re-implementation is a
strong position.

### 13d. Web API / server mode

An HTTP server wrapping ferx-core (already has an API layer in `api.rs`) would
enable:
- Language-agnostic clients (Julia, Python, R, JavaScript)
- Cloud deployment for hospital MIPD systems
- Multi-user server installations at CROs

**Effort:** 3-4 weeks for a basic REST API.
**Strategic value:** High for MIPD and clinical deployment.

---

## 14. Ecosystem integrations

### 14a. CDISC data standards (ADAM/ADPC)

Clinical trial data arrives in ADAM format (Analysis Data Model). ADPC
(Analysis Dataset for Pharmacokinetics) is the standard structure for PK
concentration data. Converting ADPC → ferx-ready CSV is a standard pre-analysis
step done by every modeller.

A `read_adpc(adpc_sas_data)` function (reading SAS7BDAT or XPT) that produces
a ferx-compatible data frame would reduce friction for clinical pharmacology.

**Effort:** 1-2 weeks (R-level, wraps haven/pharmaverseadam packages).
**Strategic value:** Medium — reduces data prep time.

### 14b. Xpose / ggPMX output compatibility

Xpose is the standard R package for PMx diagnostics (GOF plots, VPC, forest
plots). It reads NONMEM output. A `ferx_to_xpose(fit)` function that generates
an `xpdb` object compatible with Xpose from a `ferx_fit` result would
immediately give ferx all of Xpose's plotting infrastructure.

**Effort:** 2-3 weeks (mapping ferx_fit fields to xpdb structure).
**Strategic value:** High. Xpose is ubiquitous in industry. Compatibility
removes the "but my R scripts use Xpose" objection.

### 14c. NMproject / PsN workflow compatibility

NMproject (R package for NONMEM project management) and PsN (Perl toolkit) are
used by large CROs for model tracking, run management, and sensitivity analyses.
Compatibility at the project structure level (run-numbered directories, `.mod`
→ ferx translation, results in NMproject database) would ease institutional
adoption.

**Effort:** 4-6 weeks.
**Strategic value:** Medium. Relevant for enterprise/CRO adoption.

---

## Priority matrix

Organised by strategic leverage (impact on adoption) vs implementation effort.

### Tier 1: Near-term differentiators (< 3 months each)

| Item | Effort | Why |
|---|---|---|
| Non-Gaussian likelihoods: Poisson, NB | 3-4 weeks | Opens PD modelling; unique among NLME engines |
| `ferx_compare()` model table | 1 week | Immediate workflow gain; all data available |
| `ferx_vpc()` + `ferx_npde()` | 3-5 weeks | Most-used diagnostics; ferx_simulate exists |
| MAP Bayesian / `ferx_map()` | 2-3 weeks | Clinical translation; thin wrapper on EBEs |
| Ordered categorical likelihood | 4-6 weeks | Common PD endpoint; no competitor has native support |
| `ferx_covariate_gam()` | 2-3 weeks | Novel tool; pure R on top of fit |
| Python bindings (PyO3) | 4-6 weeks | Language reach |
| Xpose compatibility | 2-3 weeks | Removes toolchain objection |

### Tier 2: Medium-term structural capabilities (3-12 months)

| Item | Effort | Why |
|---|---|---|
| TTE / survival models (parametric) | 6-8 weeks | Oncology market; prerequisite for Tier 1+ |
| Bayesian HMC (SAEM+HMC polish first) | 6-8 weeks | Uncertainty QI; pediatric/rare disease |
| MIPD / `ferx_mipd()` | 4-6 weeks | Clinical deployment; unique |
| Optimal experimental design (FIM) | 6-10 weeks | Protocol design; ferx AD advantage |
| Mixture models | 8-12 weeks | NONMEM gap; metaboliser phenotyping |
| PSIS-LOO cross-validation | 2-3 weeks | Model comparison; borrows from Bayesian stats |
| `ferx_identifiability()` | 3-4 weeks | Prevents wasted runs |
| CMA-ES global optimizer | 2-3 weeks | Difficult multimodal models |

### Tier 3: Long-term research capabilities (12+ months)

| Item | Effort | Why |
|---|---|---|
| Full NUTS (population-level HMC) | 3-6 months | Gold-standard Bayes; major effort |
| Joint longitudinal + survival | 3-4 months (after TTE) | Oncology gold standard |
| Nonparametric NLME (NPML) | 6-10 weeks | Closes NONMEM gap |
| GPU acceleration | 3-6 months | Scale to N > 1000 |
| PBPK model library + Bayes | 12-24 months | Paediatric regulatory |
| Neural ODE structural models | 6-12 months | Mechanism discovery |
| Web API / server mode | 3-4 weeks | MIPD deployment |

---

## The single most leveraged addition

If only one thing were added in the next quarter, it should be **`ferx_vpc()`**
combined with **`ferx_compare()`**.

Reason: the VPC is the first plot every pharmacometrician produces after fitting
a model. Having it built-in, fast (Rust simulation), and requiring zero external
tools would eliminate the biggest friction point in the current ferx workflow.
`ferx_compare()` then makes model selection trivial. Together, a complete
model-building cycle becomes:

```r
base_fit   <- ferx_fit(base_model, data)
cov_fit    <- ferx_fit(cov_model, data)
ferx_compare(list(base=base_fit, cov=cov_fit))   # one line
ferx_vpc(cov_fit, data)                           # one plot
```

No NONMEM + PsN + Xpose pipeline. No R post-processing scripts. That is the
experience that makes someone switch tools and not go back.

---

## The moonshot: ferx as a clinical decision engine

The ultimate long-term position: ferx as the engine embedded in clinical dosing
software. The path:

1. **Now:** Fast, accurate NLME estimation. ✓
2. **Near-term:** MAP Bayesian + `ferx_mipd()` → individual dose optimisation
3. **Medium-term:** Web API → deployable in hospital pharmacy systems
4. **Long-term:** Real-time Bayesian updating as new patient observations arrive

No existing NLME tool (NONMEM, Monolix, nlmixr2) is architected for real-time
clinical deployment. ferx's speed (milliseconds per individual fit), Rust
reliability, and API-first design make this path realistic.

The regulatory modelling world and the clinical dosing world currently use
completely different software stacks. ferx could be the first tool to bridge
them.
