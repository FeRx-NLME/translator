# Direct-run feasibility: NONMEM / nlmixr2 as input language for ferx

The idea: users keep writing NONMEM control streams or nlmixr2 model functions,
submit them with their data, and ferx does the estimation. No `.ferx` syntax to
learn. Translation is internal and invisible.

This document is a pre-implementation feasibility audit. It is deliberately
pessimistic — every layer of the stack gets scrutinised for where it will break
before celebrating what already works.

---

## What already exists (the good news up front)

The translation pipeline (`ferxtranslate`) is already built and handles the core
of the problem:

```
NONMEM .ctl  →  nonmem2rx()  →  rxUi object
nlmixr2 fn   →  rxode2()     →  rxUi object
                                    ↓
                             rxui_to_ir()   →  ferx_ir
                                    ↓
                             emit_ferx()   →  .ferx text
                                    ↓
                             ferx_fit()   →  estimates
```

The data format is also already compatible: `ferx_fit()` takes NONMEM-format CSV
(ID, TIME, DV, EVID, AMT, CMT, RATE, MDV, SS, II, OCC, CENS). A NONMEM dataset
that already has named column headers works as-is.

A minimal `nm_run()` wrapper — translate, write temp `.ferx`, run `ferx_fit()` —
is roughly 30 lines of R and could be written today. That is the honest starting
point, not a limitation.

The hard work is everything below that surface.

---

## Layer 1: Data preprocessing — the first wall

### What NONMEM assumes about data

NONMEM control streams contain a `$DATA` statement that specifies:
- The path to the data file (relative to the run directory)
- `IGNORE=` rules: `IGNORE=@` (skip rows starting with `@`), `IGNORE=#`,
  `IGNORE=(C.EQN.1)`, `IGNORE=(ID.EQ.0)`, `IGNORE=FIRST` (skip first row), etc.

They contain an `$INPUT` statement that maps positional or named columns:
```
$INPUT ID TIME DV=CONC AMT EVID MDV WT CRCL DOSE=DROP
```
- `DV=CONC` means the column named CONC in the file becomes DV
- `DROP` silently discards a column
- The file may not have a header row at all — columns are positional

### What ferx_fit needs

ferx_fit needs a clean CSV with named columns matching what it expects
(ID, TIME, DV, EVID, AMT, CMT, ...). It does not apply IGNORE filters itself.

### The gap

The translator currently does not:
1. Parse `$DATA IGNORE=...` rules and filter the data file
2. Parse `$INPUT` column renaming and reorder/rename columns
3. Handle headerless NONMEM data files
4. Handle `IGNORE=(condition)` expressions

For the specific case where the user has a clean, NONMEM-standard named-header
CSV with no DROP columns and simple `IGNORE=@`, this is a non-problem. Many
academic/reference models are in this category (the test models all are).

For real clinical-trial NONMEM runs: the data file often has 30+ columns, several
DROPs, alias renames, and complex IGNORE conditions. Without a proper `$INPUT`/
`$DATA` parser in R, `nm_run()` cannot handle these.

**Verdict:** A `$DATA`/`$INPUT` preprocessing layer is required for production use.
It is buildable (parsing these blocks is straightforward), but is separate work
from the model translation and should not be mixed into the same PR.

---

## Layer 2: nonmem2rx parsing limits — what falls through before ferxtranslate even sees it

nonmem2rx is the first step; it must successfully parse the `.ctl` before the
translator sees it. nonmem2rx breaks or degrades silently on:

| NONMEM feature | Status |
|---|---|
| ADVAN1-4, ADVAN6 (ODE), ADVAN10 (Michaelis-Menten) | Mostly handled |
| `$PRED` block (user-coded `PRED`) | Not supported; nonmem2rx requires `$PK`/`$ERROR` |
| `$MIX` (finite mixture models) | Not supported |
| `NONPARAMETRIC` estimation | Not supported |
| `MTIME` (model event times) | Not supported |
| Multiple `$ERROR` endpoints tied to `DVID` | Partial; nonmem2rx may parse but translator does not map to per-CMT ferx |
| Complex `F1`, `R1`, `D1` modifiers in `$PK` | Partial; bioavailability F supported, R1/D1 infusion rate/duration less so |
| `ALAG1`, `ALAG2` absorption lags | Handled via `lagtime=` in ferx pk macros |
| `SS`/`II` steady-state dosing | Handled — ferx_fit has native SS/II column support |
| `A_0(n) = expr` initial conditions in `$PK` | nonmem2rx parses these; translator needs to check each case |
| `TABLE` output blocks | Irrelevant to model structure; not needed |
| Bayesian (`$PRIOR`) | Not supported |

**The critical failure mode:** nonmem2rx errors or silently returns a malformed UI
object, and the translator either crashes or produces wrong `.ferx`. In a silent
pipeline this is invisible. The `tryCatch` in `to_ferx()` catches hard crashes but
not malformed output that still parses.

**Mitigation:** The translation gap report (current concordance test) surfaces
`$unsupported` features, but it requires the model to have successfully parsed.
A failed nonmem2rx parse needs a separate, user-visible error path.

---

## Layer 3: Translation gaps in rxui_to_ir — what the translator cannot yet emit

Even when nonmem2rx succeeds, `rxui_to_ir()` has known gaps:

| Feature | Current status |
|---|---|
| Single-endpoint linCmt + ODE models | Working |
| 1/2/3-cpt pk macros | Working |
| Fixed thetas (FIX) | Working |
| Block omega | Working |
| IOV (kappa) | Working, diagonal only; warn on off-diagonal IOV omega |
| PKPD with multiple DVIDs / per-CMT error | **Not implemented** — ferx supports per-CMT error models but translator does not map DVID to CMT |
| Proportional / additive / combined error | Working |
| Logit-normal bioavailability (F parameter) | Partially working |
| Allometric scaling in `[individual_parameters]` | Working (pass-through expression) |
| Covariate effects (continuous, if/else) | Working |
| SDE / diffusion (ferx `[diffusion]` block) | **Not implemented** |
| Deep compartment model (`[covariate_nn]`) | **Not implemented** |
| Michaelis-Menten ODE models | Working via ODE path |
| Transit compartment absorption | Working via ODE path |
| MIXTURE / latent-class | **Not implemented in ferx-core** |

The multiple-DVID → per-CMT error gap is significant for PKPD models. A NONMEM
model with `CMT=1` PK and `CMT=2` PD observations uses DVID or CMT in `$ERROR`
to dispatch. ferx has a clean per-CMT `[error_model]` syntax. The translator does
not currently produce it.

---

## Layer 4: ferx-core feature gaps — things ferx cannot do regardless of translation

Some features are not gaps in the translator but in ferx-core itself:

| Feature | Status |
|---|---|
| FOCE / FOCEI | Supported |
| SAEM | Supported |
| Importance sampling | Supported |
| Gauss-Newton (BHHH) | Supported |
| IOV (kappa) | Supported |
| BLOQ M3 likelihood | Supported |
| Steady-state dosing | Supported |
| Multi-start | Supported (via settings) |
| MIXTURE models | **Not in ferx-core** |
| LAPLACIAN approximation | **Not in ferx-core** |
| FO (zeroth-order) estimation | **Not in ferx-core** |
| NONPARAMETRIC | **Not in ferx-core** |
| `$PRIOR` / Bayesian | **Not in ferx-core** |
| SIR (sampling importance resampling) | Supported |
| SDE (stochastic differential equations) | Supported |

When a user submits a NONMEM model that uses `METHOD=0` (FO), there is no ferx
equivalent. The run would have to silently upgrade to FOCE or refuse. Neither is
clean.

---

## Layer 5: Estimation method mapping

The method mapping for supported methods is clean and unambiguous:

| NONMEM `$ESTIMATION` | ferx method |
|---|---|
| `METHOD=1 INTERACTION` | `focei` |
| `METHOD=1` (no INTERACTION) | `foce` |
| `METHOD=SAEM` | `saem` |
| `METHOD=IMP` | `imp` |
| `METHOD=0` (FO) | **no equivalent — hard-stop** |
| `LAPLACIAN` | **no equivalent — hard-stop** |

For all mappable methods, numerical results will differ from NONMEM — and that is
the point. ferx uses exact automatic differentiation (Enzyme) where NONMEM uses
central finite differences. Converged estimates are more accurate; standard errors
are tighter because the curvature estimate is exact. This is an advantage to
communicate, not a problem to apologise for.

The only exception requiring a hard-stop is `METHOD=0` (FO). FO is not in
ferx-core and silently upgrading to FOCE changes the statistical model and the
interpretation of parameter estimates. The wrapper must stop with a clear message.

---

## Layer 6: Output format — ferx already has more than NONMEM

After reading `ferx_fit()` documentation in full, the output is already substantially
richer than NONMEM. The familiar NONMEM outputs are all covered:

| NONMEM output | ferx equivalent |
|---|---|
| OFV (`.ext`, final row) | `$ofv` |
| THETA/OMEGA/SIGMA estimates | `$theta`, `$omega`, `$sigma` |
| Standard errors, %RSE | `$se_theta`, `$se_omega`, `$se_sigma` |
| Covariance matrix (`.cov`) | `$cov_matrix` |
| Correlation matrix (`.cor`) | `ferx_cor_matrix()` |
| Condition number | `$condition_number` |
| Eigenvalues | `$eigenvalues` |
| Individual ETAs (`.phi`) | `$ebe_etas` |
| Individual predictions, CWRES, IWRES | `$sdtab` |
| Per-subject OFV | `$sdtab$EBE_OFV` |
| ETA shrinkage | `$shrinkage_eta` |
| EPS shrinkage | `$shrinkage_eps` |
| Iteration trace (`.ext` body) | `$trace_path` (CSV) |
| Convergence status | `$converged`, `$covariance_status` |
| Wall time | `$wall_time_secs` |

The only real gap for NONMEM users is visual familiarity — the output is an R list
rather than flat files with NONMEM-conventional column names. A formatter is a UX
layer, not a capability gap. See the separate output-advantage plan for detail.

The more important story is what ferx already has that NONMEM does NOT — also
covered in the separate plan.

---

## Layer 7: The debuggability problem — what happens when it goes wrong

In the current workflow, the user sees the `.ferx` file. When something looks
wrong, they can read it, compare it to their NONMEM model, and spot the
translation error. The `.ferx` file is a legible audit trail.

In a "silent pipeline" `nm_run()` workflow, the `.ferx` is internal. When a model
produces biologically implausible estimates, the user has no visibility into what
was actually run. They cannot tell if the problem is:
- A translation bug (wrong structural model)
- A ferx estimation issue (convergence to a local minimum)
- A data issue (wrong column mapping)
- An expected numerical difference from NONMEM

**Mitigation:** The `nm_run()` result must include the intermediate `.ferx` text,
accessible as `result$ferx_text` and printable/writable on demand. The translation
warnings and unsupported-features list must be prominently exposed, not buried.
The function should have a `verbose=TRUE` mode that prints the `.ferx` to console
before fitting.

A `show_ferx=TRUE` argument (default FALSE, but recommended in documentation)
that writes the intermediate `.ferx` to the working directory next to the input
`.ctl` would be valuable: it externalises the audit trail.

---

## Layer 8: The nlmixr2 R-object advantage — the underexplored story

This is the most important structural difference between the NONMEM path and the
nlmixr2 path, and it changes what ferx can offer.

### What an nlmixr2 model actually is

An nlmixr2 model function is not a configuration file or a string. When passed to
`rxode2()` or `nlmixr2()`, it produces an `rxUi` object — a structured R list
with the following key fields:

- `$iniDf`: a **data frame** of all parameters (thetas, omegas, sigmas, kappas)
  with columns `name`, `est`, `lower`, `upper`, `fix`, `neta1`, `neta2`,
  `condition`. This is directly programmable.
- `$lstExpr`: a **list of R language objects** — actual parse tree nodes (AST),
  not strings. `deparse(ui$lstExpr[[1]])` gives `"cl <- tvcl * exp(eta.cl)"`.
  These are first-class R objects you can inspect, traverse, substitute, or build
  from scratch using `quote()`, `bquote()`, `rlang::call2()`.
- `$theta`, `$omega`, `$sigma`: matrices and named vectors — standard R.
- `$covariates`: character vector of covariate names referenced in the model.

A NONMEM `.ctl` file is a text format. To manipulate it you must re-parse text.
An nlmixr2 model is an R-native object. Manipulation is R-native.

### What the R-object nature enables

**1. Programmatic model modification before translation**

nlmixr2 already has update methods: `ini(model, ...)` updates initial estimates,
`model(model, ...)` replaces model block expressions. The `$iniDf` data frame can
be mutated directly. A covariate model can be composed in R:

```r
base <- function() {
  ini({ tvcl <- 0.134; eta.cl ~ 0.07; prop.err <- 0.01 })
  model({ cl <- tvcl * exp(eta.cl); linCmt() ~ prop(prop.err) })
}
ui <- rxode2::rxode2(base)
# Add WT power covariate on CL — pure R, no text editing
ui2 <- rxode2::model(ui, cl <- tvcl * (WT/70)^0.75 * exp(eta.cl))
# Translate both to ferx, run, compare
```

**2. The ferx_ir is already a programmable R object**

The `ferx_ir` produced by `rxui_to_ir()` is itself a structured R list
(`new_ferx_ir()`). It contains `$thetas`, `$omegas`, `$indiv_params`, `$odes`,
etc. as plain R lists. You can inspect and modify it between translation and
emission:

```r
ir  <- rxui_to_ir(ui)
# Inspect what was translated
ir$indiv_params
# Fix a parameter after the fact
ir$thetas[[1]]$fixed <- TRUE
# Emit and run
ferx_fit(emit_ferx(ir), data)
```

**3. The modify-translate-run loop is entirely R-native**

The pipeline `nlmixr2 fn → rxUi → ferx_ir → .ferx text → ferx_fit()` has two
programmable R objects in it. This means:

```
write nlmixr2 fn
     → rxode2() → rxUi            (inspectable R object)
     → rxui_to_ir() → ferx_ir    (modifiable R object)
     → emit_ferx() → .ferx text  (only materialises here)
     → ferx_fit() → ferx_fit obj (inspectable R object)
     → update inits, fix params, add covariate
     → repeat
```

Every step is R. No file editing, no text patching. This is a qualitatively
different workflow from "edit the .ctl, submit to NONMEM, parse .ext output."

**4. What this enables that NONMEM cannot**

| Capability | NONMEM | nlmixr2 + ferx |
|---|---|---|
| Modify parameters programmatically | Edit text | Modify `iniDf` data frame |
| Add/remove covariates | Edit $PK text | Modify `$lstExpr` AST |
| Compare two model variants | Text diff | Compare `ferx_ir` lists field-by-field |
| Generate model from R data | Write text template | Build `ferx_ir` or rxUi directly |
| Automate stepwise covariate search | External R scripts patching text | Pure R loop over `ferx_ir` variants |
| Update initials from previous fit | Parse `.ext` file, edit `.ctl` | `ir$thetas <- update_from_fit(fit)` |

**5. The ferx_ir as a model builder API**

`new_ferx_ir()` already works as a direct API. Users who do not want to write
NONMEM or nlmixr2 syntax at all can build a model entirely in R:

```r
ir <- new_ferx_ir(
  thetas      = list(list(name="TVCL", init=0.134, lower=0.001, upper=10),
                     list(name="TVV",  init=8.1,   lower=0.1,   upper=500)),
  omegas      = list(list(type="diagonal", names="ETA_CL", values=0.07)),
  sigmas      = list(list(name="PROP_ERR", value=0.01, scale="sd")),
  indiv_params= list(list(lhs="CL", rhs="TVCL * exp(ETA_CL)")),
  structural  = list(type="pk_macro", pk_call="one_cpt_oral",
                     pk_args=list(cl="CL", v="V", ka="KA")),
  error_model = list(list(dv="DV", type="proportional", params="PROP_ERR")),
  fit_options = list(method="focei", maxiter=300L)
)
ferx_fit(emit_ferx(ir), data)
```

This is a ferx model with no `.ferx` syntax, no NONMEM, no nlmixr2. Pure R. The
`ferx_ir` is already the native R model API — it just isn't documented as such.

### What is currently missing to make this story complete

The `ferx_ir` modification API exists implicitly (you can manipulate the list
fields) but there are no helper functions for common operations:

| Missing function | What it does |
|---|---|
| `ir_add_covariate(ir, param, cov, model, ref)` | Adds a covariate relationship to an `[individual_parameters]` expression |
| `ir_fix_theta(ir, name)` | Marks a theta as FIX |
| `ir_update_inits(ir, fit)` | Replaces theta/omega/sigma initials with converged estimates from a previous `ferx_fit` result |
| `ir_set_method(ir, method)` | Switches estimation method in `fit_options` |
| `ir_compare(ir1, ir2)` | Prints a structural diff of two `ferx_ir` objects |
| `ir_from_fit(fit)` | Reconstructs an `ferx_ir` with initials set to a fit's estimates (for sequential fitting) |

These are each 10-30 lines of R. Together they form a proper ferx model builder API.

### The nlmixr2 engine registration — still Phase 5

Registering ferx as an `est="ferx"` engine inside nlmixr2 proper (so
`nlmixr2::nlmixr2(model, data, est="ferx")` works) remains a longer-term goal
that requires coordination with the nlmixr2 team. But the R-object workflow
described above does not require that registration — it works today through
`rxui_to_ir()` directly. The engine registration is the last-mile UX polish, not
the capability unlock.

---

## Feasibility summary by phase

### Phase 1 — Minimal wrapper (feasible now, ~1 week)

**What:** A single `nm_run(ctl_file, data_file, ...)` function in ferxtranslate.

**What it does:**
1. Calls `nm_to_ferx(ctl_file)` internally
2. If `result$unsupported` is non-empty: stops with a clear, itemised error.
   Does NOT silently proceed with a broken translation.
3. If warnings exist: prints them prominently before fitting.
4. Writes `.ferx` to a tempfile (and optionally to disk via `show_ferx=`).
5. Calls `ferx_fit(ferx_file, data_file, ...)` and returns the fit.
6. Attaches `ferx_text`, `warnings`, and `unsupported` to the returned object.

**What it does NOT do:**
- Does not preprocess NONMEM data (IGNORE filters, $INPUT renaming).
  The data_file must already be a clean NONMEM-format named-header CSV.
- Does not reformat output to look like NONMEM output.
- Does not handle models that nonmem2rx cannot parse.

**Scope check:** This is a thin wrapper. It should be treated as a convenience
function, not a compatibility claim. Documentation must be explicit: "works for
models that `nm_to_ferx()` already supports." No new translation logic in Phase 1.

### Phase 2 — Data preprocessing layer (~2-4 weeks)

**What:** A `nm_preprocess_data(ctl_file, data_file)` function that:
- Parses `$DATA` for the IGNORE= condition
- Parses `$INPUT` for column renaming / DROP
- Applies IGNORE filters row by row
- Returns a clean data frame ready for ferx_fit

This makes `nm_run()` work for real-world NONMEM datasets, not just clean
reference ones.

**Hard sub-problems:**
- `IGNORE=(condition)` with compound expressions: `IGNORE=(ID.EQ.0.AND.TIME.GT.0)`
  requires a mini expression evaluator
- Headerless data files: positional column assignment from $INPUT order
- The `$DATA` path is relative to the run directory, not the R working directory

### Phase 3 — NONMEM-style output formatter (~2-3 weeks)

**What:** A `print.nm_run_result` method that presents:
- OFV, number of observations, number of subjects
- THETA table with initial / final / RSE% / 95% CI
- OMEGA lower triangular with RSE%
- SIGMA with RSE%
- Condition number, shrinkage per ETA
- Convergence message

This is display-only work. The underlying data is already in the ferx_fit result.

### Phase 4 — PKPD / multi-endpoint translation (~3-6 weeks)

**What:** Extend `rxui_to_ir()` to map multiple-DVID NONMEM models to ferx's
per-CMT `[error_model]` syntax.

This requires:
- Detecting multiple error expressions keyed by CMT or DVID in nonmem2rx's output
- Mapping NONMEM's DVID convention to ferx CMT indices
- Emitting `CMT=N: DV ~ ...` error model lines
- Testing against the `pkpd_ir.mod` model already in the test suite

### Phase 5 — ferx_ir modifier API (~1-2 weeks, high leverage)

**What:** A small set of helper functions for programmatic `ferx_ir` manipulation.
Because ferx_ir is already a plain R list, these are thin wrappers — the value
is documentation and discoverability, not complexity.

Functions to ship in this phase:
- `ir_update_inits(ir, fit)` — replace theta/omega/sigma initials from a converged
  `ferx_fit` result; enables warm-restarts and sequential model building
- `ir_fix_theta(ir, name, unfix=FALSE)` — mark/unmark a theta as FIX
- `ir_add_covariate(ir, param, covariate, model="power", reference=NULL)` — appends
  a covariate relationship to the relevant `$indiv_params` expression; supports
  `"power"`, `"proportional"`, `"linear"`, `"categorical"` model types
- `ir_set_method(ir, method)` — switches `fit_options$method`
- `ir_compare(ir1, ir2)` — prints a field-by-field structural diff

Document `new_ferx_ir()` explicitly as a model builder, not just an internal
constructor. Add examples showing end-to-end model definition in pure R.

**Why now:** This phase unlocks the programmatic workflow described in Layer 8
without any external dependencies or coordination. It is the highest-leverage work
relative to effort in the roadmap.

### Phase 6 — nlmixr2 engine integration (months, external collaboration)

As described in Layer 8 above. The R-object workflow (Phases 1-5) works without
this. Phase 6 is last-mile UX for users who want `nlmixr2::nlmixr2(model, data, est="ferx")`.

---

## What will never work, and should not be promised

| Feature | Why |
|---|---|
| Exact NONMEM OFV reproduction | Different gradient method; AD vs FD by design |
| `$PRED`-only models | nonmem2rx does not parse; no path to ferx |
| MIXTURE / latent class | Not in ferx-core |
| NONPARAMETRIC | Not in ferx-core |
| `$PRIOR` / Bayesian | Not in ferx-core |
| FO estimation (`METHOD=0`) | Not in ferx-core; silent upgrade to FOCE is a lie |
| MTIME event scheduling | Not in nonmem2rx, not in ferx |
| Producing NONMEM table files (`.tab`, `.sdtab`) | Different output object |
| Bit-for-bit identical results to NONMEM | Not the goal; not possible |

---

## Architecture decision: where does this live?

Two options:

**Option A: In ferxtranslate** — `nm_run()` and `nlmixr2_run()` added to the
existing package. Clean: keeps all translation and execution in one place.
The ferxtranslate package becomes both the translator and the run wrapper.

**Option B: Separate package** — A `ferxrun` or `ferxcompat` package that depends
on both ferxtranslate and ferx. Cleaner separation of concerns. ferxtranslate
stays as a pure translation library. But adds a package dependency layer that
fragments the user experience.

**Recommendation:** Option A for Phase 1 and 2. If the run-wrapper grows a
substantial codebase of its own (output formatters, data preprocessors, nlmixr2
integration), revisit splitting at that point.

---

## The honest risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Users mistake "runs in ferx" for "identical to NONMEM" | High | High | Prominent documentation + numeric difference warning in every run |
| Silent translation errors produce wrong parameter estimates | Medium | High | Hard-stop on any `$unsupported`; verbose warning on WARN-level gaps |
| Real-world data files break on IGNORE filters | High | Medium | Phase 2 preprocessing layer; Phase 1 documents the data requirement |
| nonmem2rx fails or mis-parses complex models | Medium | High | Wrapper catches nonmem2rx errors cleanly; does not pass malformed UI downstream |
| Users skip reading translation warnings | High | Medium | Warnings cannot be suppressed by default; require explicit `suppress_warnings=TRUE` |
| OFV difference questions dominate support load | High | Medium | Reference document explaining AD vs FD, link in every print.nm_run_result |
| Phase 4 PKPD work is larger than estimated | Medium | Low | It is scoped as its own phase; not blocking Phase 1-3 |

---

## Recommended starting point

Start with Phase 1. Write `nm_run()` and `nlmixr2_run()` with the following
constraints enforced by design:

1. **Hard-stop on unsupported**: `if (length(result$unsupported) > 0) cli::cli_abort(...)`.
   No silent proceeds on known translation failures.

2. **Warnings always visible**: do not pipe WARN-level translation warnings into
   a suppressed log. Print to console before calling ferx_fit.

3. **ferx_text always accessible**: expose the intermediate `.ferx` as a field on
   the returned object. Print it when `verbose=TRUE`.

4. **Numeric difference disclaimer**: print a one-line note on every successful
   fit: "ferx results may differ numerically from NONMEM due to exact-AD vs
   finite-difference gradients; see vignette('nonmem-differences')."

5. **Data requirement documented**: document that data_file must be a named-header
   NONMEM-format CSV with IGNORE rows already removed. Do not silently accept
   files that will silently produce wrong fits.

Do not attempt Phases 2-6 in the same PR as Phase 1. Each phase is independently
testable and releasable.

---

## Phase ordering rationale

Phase 5 (ferx_ir modifier API) could be done before Phase 2-4 and still deliver
substantial value. The modifier API enables programmatic workflows for nlmixr2
users immediately, without requiring the NONMEM data preprocessing or output
formatting work. If the target audience is "nlmixr2 users who want to run in
ferx," the priority order should be: 1 → 5 → 3 → 2 → 4 → 6.

If the target audience is "NONMEM users who want zero friction entry to ferx,"
the order is: 1 → 2 → 3 → 4 → 5 → 6.
