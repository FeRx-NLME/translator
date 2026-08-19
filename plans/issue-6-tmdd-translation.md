# Issue #6 -- TMDD model translates to an invalid .ferx

Evaluation and implementation plan.

Status: **evaluated, not implemented**. Every claim below was reproduced or
disproved against `ferxtranslate` @ `2a7df47` (main), `ferx` 0.2.0,
`nonmem2rx` 0.1.9, R 4.5.1. Baseline test suite before any change:
`FAIL 0 | WARN 5 | SKIP 1 | PASS 271` (123.7 s, concordance included).

---

## 1. Verdict

All seven reported defects reproduce exactly as described. The issue
under-reports the problem in four ways:

1. **The blast radius is wider than TMDD.** The root cause of defect 4 (`$DES`
   conditional emits an undefined `cf`) is that *every* `if` statement in the
   model body is silently discarded, wherever it appears. `IF (SEX.EQ.1) TVCL =
   THETA(1)*THETA(3)` in `$PK` -- the most common categorical-covariate idiom in
   NONMEM -- is dropped with zero warnings. Any covariate model written that way
   currently translates to a wrong model that validates clean.

2. **Seven further defects exist that the issue does not mention**, six of them
   producing a model that fits and returns *wrong numbers* rather than failing
   loudly.

3. **The worst defect found is not in the issue at all, and it is live in our
   own test suite.** In every ferx block where thetas are in scope
   (`[individual_parameters]`, `[scaling]`), a name that matches a declared
   theta always resolves to the theta -- even after being redefined on a
   preceding line. Because `nonmem2rx` names thetas after the `$PK` parameters
   whenever `$THETA` carries labels, we routinely emit `theta CL` alongside
   `CL = CL * exp(ETA_CL)`, and the second line is **dead**. Numerically proven
   in section 3.14: `pk_1cmt_oral.mod` -- a tier-4 concordance model -- is
   currently fitted with no inter-individual variability on clearance at all,
   and the test passes because it only asserts structural thetas.

4. **One bundled test model is already wrong for a defect the issue does
   report.** `pkpd_ir.mod` contains `A_0(4)=BL`. It is dropped, the effect
   compartment starts at 0 instead of 100, and the output passes
   `ferx_model_validate()`. Defect 3 is shipping in our own corpus today.

The good news is decisive: **no ferx-core change is required.** A faithful,
lossless translation of the reporter's model is expressible in ferx 0.2.0 as it
stands, including the two-endpoint `FLAG` dispatch. This is entirely a
`ferxtranslate` bug.

The issue's suggested fix for defect 5 (emit per-`CMT` blocks and have the user
remap their dataset) should **not** be implemented -- ferx dispatches on any
column directly. See section 5.3.

---

## 2. Reproduction

Reporter's model saved as `qss_tmdd.mod`, three-row `data.csv`, then
`nm_to_ferx("qss_tmdd.mod", output = "out.ferx")`. Output (abridged):

```
[individual_parameters]
  KEL = KEL * exp(ETA1)
  VC = VC * exp(ETA2)
  RBASE = RBASE * exp(ETA3)
  KSYN = RBASE * KDEG
  FB = cf/(KSS + cf)          # cf undefined; FB is state-dependent
  W1 = 0                      # $ERROR scaffolding
  W2 = 0

[structural_model]
  ode(obs_cmt=c.RTOT, states=[CENT, TISS, c.RTOT])

[odes]
  d/dt(CENT) = -(KEL + KPT) * cf * VC - c.RTOT * KINT * FB * VC + KTP * TISS
  ...

[error_model]
  DV ~ combined(EPS1, EPS2)
```

One warning (the `obs_cmt` guess); `unsupported` empty. All seven confirmed.

The decisive artefact is `ui$lstExpr` from `nonmem2rx`: **nothing is missing
from the input.** Every dropped feature is present and recoverable.

```
[13] rxini.rxddta3. <- rbase        # A_0(3) -- present
[14] c.RTOT(0) <- rxini.rxddta3.    #        -- present, dropped by translator
[15] scale1 <- vc                   # S1=VC  -- present
[20] if (dsc < 0) dsc <- 0          # dropped
[22] if (bb >= 0) { cf <- ... } else { cf <- ... }   # dropped
[27] f <- CENT/scale1               # the drug readout -- present
[31] if (FLAG == 2) ipred <- rtot   # endpoint dispatch -- dropped
[34] if (FLAG == 1) w1 <- 1         # dropped
[36] y <- ipred * (1 + w1*eps1 + w2*eps2)
```

Every defect is a translator omission, not an upstream parsing loss.

---

## 3. Defect inventory

Severity uses the only distinction that matters: does the user find out?

- **SILENT-WRONG** -- output validates and fits, answers are wrong.
- **LOUD** -- output fails `ferx_model_validate()`; the user is blocked, not
  misled.

| # | defect | class | source |
|---|---|---|---|
| 14 | theta name shadows identically-named individual parameter | **SILENT-WRONG** | new |
| 1 | state id `c.RTOT` contains a dot | LOUD | issue |
| 2 | thetas used only in `$DES` get no pass-through | LOUD | issue |
| 3 | `A_0(n)` initial condition dropped | **SILENT-WRONG** | issue |
| 4 | `IF/THEN/ELSE` dropped, emits undefined `cf` | LOUD here, SILENT-WRONG in `$PK` | issue |
| 5 | two endpoints collapse to one, error mistyped | **SILENT-WRONG** | issue |
| 6 | `$ERROR` indicator vars leak to `[individual_parameters]` | cosmetic | issue |
| 7 | output never validated before returning | meta | issue |
| 8 | `IF` in `$PK` dropped -- categorical covariates vanish | **SILENT-WRONG** | new |
| 9 | `F1` dropped; `ALAG1` emitted as an orphan parameter | **SILENT-WRONG** | new |
| 10 | `combined()` args in traversal order, not (prop, add) | **SILENT-WRONG** | new |
| 11 | proportional error written additively classified `additive` | **SILENT-WRONG** | new |
| 12 | `S1=VC` scaling dropped whenever `obs_cmt` is guessed wrong | **SILENT-WRONG** | new |
| 13 | `.inline_aux_vars()` depth-30 cutoff silently un-inlines | LOUD | new |
| 15 | `obs_scale` + Form C `y` double-scale if both emitted | **SILENT-WRONG** | new (design constraint) |

### 3.14 Theta shadowing (new, worst)

`parse_atom` resolves identifiers in strict precedence order --
**theta, then eta, then defined variables (individual parameters, states), then
covariate** -- `ferx-core/src/parser/model_parser.rs:16026-16050`:

```rust
// Check if it's a theta
if let Some(idx) = ctx.theta_names.iter().position(|n| n == name) {
    return Ok((Expression::Theta(idx), pos + 1));
```

So in any block where thetas are in scope, a redefinition is silently dead.
Minimal reproduction -- one compartment, dose 100, `PRED` at t = 1:

```
[parameters]                    | [parameters]
  theta CL(1.0, ...)            |   theta CL(1.0, ...)
  theta V(10.0, ...)            |   theta V(10.0, ...)
[individual_parameters]         | [individual_parameters]
  CL = CL * 5                   |   iCL = CL * 5
  K20 = CL / V                  |   K20 = iCL / V
```
```
PRED = 90.48374   (= 100*exp(-0.1),  K20 = 1/10  -- the `* 5` was ignored)
PRED = 60.65307   (= 100*exp(-0.5),  K20 = 5/10  -- correct)
```

`ferx_model_validate()` reports **no diagnostic at all** for the left-hand
model. Phase 1 will not catch this one; only the rename will.

**Scope of the bite.** The shadow only matters where a shadowed name is
*referenced* in a theta-visible block. It does **not** affect `[odes]` (thetas
out of scope) and it does **not** affect `pk ...(cl=CL, ...)` arguments, which
are resolved by direct name lookup against the individual-parameter list rather
than through `parse_atom`.

Two bundled models emit a theta whose name matches an individual-parameter LHS:

| model | shadowed names | actually wrong? |
|---|---|---|
| `pk_1cmt_oral.mod` | `CL`, `KA` | **yes** -- `K20 = CL/V` sits inside `[individual_parameters]` |
| `pk_1cmt_oral_ampsim.ctl` | `CL`, `KA`, `V` | no -- latent; every reference is a pk macro arg |

Both are in the tier-4 concordance suite. Verified for ampsim by the same
dead-line test below: its `CL` line *is* live (`PRED` moves from `0.0776` to
`0.0361`), so the pattern is present but harmless there today. It becomes
harmful the moment anyone adds a derived parameter to that model.

Proof for `pk_1cmt_oral.mod`, using its own concordance dataset. Take the
emitted `.ferx` and replace `CL = CL * exp(ETA_CL)` with `CL = CL * 5`:

```
pk_1cmt_oral.mod.ferx  PRED[1:3] = 0.019409423 0.030702439 0.040499662
pk1_dead.ferx          PRED[1:3] = 0.019409423 0.030702439 0.040499662
```

Multiplying clearance by five changes nothing. The line is dead; `K20 = CL/V`
reads the theta. `ETA_CL` therefore has no effect on the model, and the tier-4
test passes only because it asserts structural thetas, not omegas.

Note the scope is exactly the blocks where thetas are visible. Inside `[odes]`
thetas are **not** in scope, so the same name there correctly resolves to the
individual parameter -- which is why the ODEs are right and the derived
parameter is wrong, in the same file.

**Fix:** never emit an individual parameter whose name matches a declared theta.
Rename the theta, not the parameter -- `theta TVCL` + `CL = TVCL * exp(ETA_CL)`
is idiomatic ferx and matches the existing style in `1cpt_oral.ctl`. Fall back
to `THETA_<name>` then a numeric suffix if `TV<name>` is taken.

### 3.8 `IF` in `$PK` (new, severe)

```
$PK
  TVCL = THETA(1)
  IF (SEX.EQ.1) TVCL = THETA(1)*THETA(3)
  CL = TVCL*EXP(ETA(1))
```
emits `TVCL = THETA1` / `CL = TVCL * exp(ETA1)`, with `warnings` and
`unsupported` both empty. The sex effect is gone; `THETA3` is estimated against
nothing. ferx supports the construct natively -- verified valid, with and
without a trailing `else`:

```
  if (SEX == 1) { TVCLI = TVCL * COVF }
  if (RACE == 1) { KF = 1 } else if (RACE == 2) { KF = 2 } else { KF = 3 }
```

### 3.9 `F1` / `ALAG1` (new, severe)

`F1 = THETA(4)` vanishes entirely -- `f(depot) <- ...` has a non-symbol LHS and
hits the silent `next` at `rxui_to_ir.R:432`. `ALAG1` survives as a dead
`RXALAG_DEPOT_` parameter because `.infer_pk_macro()` looks for an lhs named
`alag`/`lagtime`/`tlag` and `nonmem2rx` names it `rxalag.depot.`.
`pk one_cpt_oral(..., f=F, alag=ALAG)` is accepted by ferx. Translator-side
naming bug only.

This is the case `ferx_model_validate()` catches for free:

```
W_UNUSED_PARAM: theta 'THETA4' is declared in [parameters] but not referenced
W_PARSE: [individual_parameters] `RXALAG_DEPOT_` is computed but never used
```

Two warnings the user never saw, because we never ran the validator.

### 3.10 / 3.11 Error-model misclassification (new, severe)

`.classify_error_assignment()` counts epsilons instead of analysing structure.

| `$ERROR` | correct | emitted |
|---|---|---|
| `Y = F + EPS(1) + F*EPS(2)` | `combined(EPS2, EPS1)` | `combined(EPS1, EPS2)` -- **swapped** |
| `Y = F + F*EPS(1)` | `proportional(EPS1)` | `additive(EPS1)` |

ferx `combined(a, b)` is `(Proportional, Additive)` --
`ferx-core/src/types.rs:1606`:
```rust
ErrorModel::Combined => vec![SigmaType::Proportional, SigmaType::Additive],
```
When the additive term is written first, the two sigmas are transposed: the
additive SD is applied proportionally and vice versa. The fit converges and the
answer is wrong.

### 3.12 `[scaling]` lost with a wrong `obs_cmt` (new)

`scaling_hint` is keyed on the *index of `obs_cmt` in the state list*
(`rxui_to_ir.R:77-79`). The reporter's model has `S1=VC`, but `obs_cmt` was
guessed as `c.RTOT` (index 3), so `scaling_hint[["3"]]` is `NULL` and no
`[scaling]` block is emitted. Concentrations are then predicted as amounts.
Defects 5 and 12 compound.

Note the same code prefers a theta match over an individual-parameter match
(`rxui_to_ir.R:84-88`), which walks straight into defect 14.

---

## 4. Root causes

Fifteen symptoms, six causes.

**RC-A -- `if` statements are not parsed at all.** `.parse_model_exprs()`
(`rxui_to_ir.R:392`) handles only `~` and `.is_assignment()`. An `if` call
matches neither and falls off the end of the loop with no warning and no
`unsupported` entry. Causes defects 4, 6, 8, half of 5.

**RC-B -- identifiers are not sanitised to the ferx grammar.** Only names in
`name_map` (built from `iniDf`) pass through `.norm()`; state names and any
unregistered symbol are emitted verbatim. ferx identifiers are
`[A-Za-z_][A-Za-z0-9_]*` (`model_parser.rs:15776-15784`). `c.RTOT` is one
instance of a general `nonmem2rx` behaviour -- it prefixes `c.` onto a
compartment name that collides with a variable name (here the `RTOT` in
`$ERROR`). Causes defect 1.

**RC-C -- non-symbol assignment targets are silently discarded.**
`rxui_to_ir.R:432`: `if (!is.symbol(lhs_expr)) next`. One line drops
`STATE(0) <- ...`, `f(STATE) <- ...`, `alag(STATE) <- ...`, `dur()`, `rate()`.
Causes defects 3 and 9.

**RC-D -- the emitter never checks its own output.** `ferx` is not even in
`Suggests`. Causes defect 7, and is why 1, 2, 4, 9 reached a user instead of a
test.

**RC-E -- `$ERROR` is classified by epsilon-counting, not by structure.** Both
`.classify_error_assignment()` and `.parse_error_rhs()` fall back to
"proportional, plus a soft WARN" when they do not recognise a shape. That
fallback is the disease: a guess that looks like a translation. Causes 5, 10, 11.

**RC-F -- emitted names are not checked against ferx's scoping rules.** The
translator picks names by normalising the source and never asks whether the
result means what it says in the block it lands in. Causes defect 14, and makes
15 possible.

Enabling all of these: **the IR cannot represent an ordered statement list.**
`indiv_params` is a flat `list(lhs, rhs-string)` and `odes` a flat
`list(state, rhs-string)`. Neither can hold an `if`, an `init()`, or an
intermediate at a required position. This must change before RC-A and RC-C can
be fixed properly.

---

## 5. Ground truth: what ferx 0.2.0 actually supports

Verified empirically against installed `ferx` 0.2.0, cross-read against
`ferx-core` @ `54a7d44`.

### 5.1 Everything needed is already there

| need | ferx construct | evidence |
|---|---|---|
| initial condition (ODE) | `init(<state>) = <expr>` in `[odes]`, any position | `model_parser.rs:8629`, validated |
| initial condition (analytic) | `[initial_conditions]` block | `model_parser.rs:6839` |
| `$DES` intermediates | plain assignments in `[odes]` | `model_parser.rs:8901` |
| `$DES` / `$PK` conditionals | `if {} else if {} else {}` in both blocks, `else` optional | `model_parser.rs:16299`, validated |
| inline conditional value | `if (cond) a else b`, `else` mandatory | `model_parser.rs:15945` |
| bioavailability / lag | `pk ...(f=F, alag=ALAG)` (`lagtime=` is an alias) | `types.rs:692-710`, validated |
| non-`CMT` endpoint dispatch | covariate-selected `[error_model]` if/else (#658) | `model_parser.rs:10470`, validated |
| non-`CMT` prediction dispatch | `[scaling]` Form C `y = if (FLAG==2) ... else ...` | validated |
| `CMT`-based dispatch | `y[CMT=N]` + `CMT=N: DV ~ ...` | `model_parser.rs:6471`, `10606` |
| no observed compartment | `ode(states=[...])` with no `obs_cmt`, **only** when `y` is given | `model_parser.rs:2443` |

### 5.2 Target output for the reporter's model -- validated

Faithful, no data changes. `ferx_model_validate()` returns `ok = TRUE` with zero
diagnostics, model-only and with the three-row dataset.

Note the `i` prefix on every individual parameter: it is **required**, not
cosmetic. Without it the `y` readout's `VC` resolves to the theta and the IIV
on volume silently disappears from the concentration endpoint (defect 14,
confirmed numerically -- `y = CENT/VC` gave `PRED = 24.57` against `CENT/iVC`'s
`12.29`, a clean factor of two).

```
[individual_parameters]
  iKEL   = KEL * exp(ETA1)
  iVC    = VC * exp(ETA2)
  iKTP   = KTP                      # theta pass-through (defect 2)
  iKPT   = KPT
  iKSS   = KSS
  iKINT  = KINT
  iKDEG  = KDEG
  iRBASE = RBASE * exp(ETA3)
  iKSYN  = iRBASE * iKDEG

[structural_model]
  ode(states=[CENT, TISS, C_RTOT])     # no obs_cmt: the y readout defines it

[odes]
  init(C_RTOT) = iRBASE                                 # defect 3
  CT  = CENT/iVC                                        # defect 4: intermediates,
  RT  = C_RTOT                                          #   not inlining
  BB  = CT - RT - iKSS
  DSC = BB * BB + 4 * iKSS * CT
  if (DSC < 0) { DSC = 0 }
  DD  = sqrt(DSC)
  if (BB >= 0) { CF = 0.5 * (BB + DD) } else { CF = 2 * iKSS * CT/(DD - BB + 1e-30) }
  FB  = CF/(iKSS + CF)
  d/dt(CENT)   = -(iKEL + iKPT) * CF * iVC - RT * iKINT * FB * iVC + iKTP * TISS
  d/dt(TISS)   = -iKTP * TISS + iKPT * CF * iVC
  d/dt(C_RTOT) = iKSYN - iKDEG * RT - (iKINT - iKDEG) * RT * FB

[scaling]
  y = if (FLAG == 2) C_RTOT else CENT/iVC               # defects 5 + 12

[error_model]
  if (FLAG == 1) { DV ~ proportional(EPS1) }            # defect 5
  else { DV ~ proportional(EPS2) }
```

(The recommended fix in 3.14 renames the *theta* rather than the parameter, so
the shipped form would be `theta TVVC` + `VC = TVVC * exp(ETA2)`. Either
direction removes the shadowing; the `i` prefix is used above only to keep the
diff against the reporter's names readable.)

### 5.3 Correction to the issue's suggested fix for defect 5

The issue proposes per-`CMT` blocks plus a user-side data remap from `FLAG` to
`CMT`. **Do not implement this.** ferx's covariate-selected error model
(ferx-core #658) plus a Form C `y` readout dispatches on any data column
directly, dataset unchanged. A remap would be a lossy, user-hostile fix for a
problem ferx does not have.

Both forms are still needed: dispatch on `CMT` -> per-`CMT` form; dispatch on
any other column -> selector form. Decision rule in phase 6.

### 5.4 Traps that constrain the design

All verified; several numerically.

- **Theta shadowing** (defect 14, section 3.14). Applies to
  `[individual_parameters]`, `[scaling]` (`y` and `obs_scale` alike), and
  `[initial_conditions]`. Not to `[odes]`, where thetas are out of scope, and
  not to `pk ...(cl=CL)` arguments, which resolve by direct name lookup. No
  diagnostic in any case.
- **`obs_scale` stacks with a Form C `y` -- silently double-scaling.**
  `apply_scaling` guards only on `ScalingSpec::None` (`src/pk/mod.rs:177`);
  there is no Form-C guard, despite a comment at `:161-164` claiming Form C
  does not reach it. Measured on identical data: `y = CENT/VC` -> 28.108;
  `y = CENT` + `obs_scale = VC` -> 28.108; **both together -> 9.369**
  (divided twice). All validate clean. **Emit one or the other, never both.**
- **`y = <expr>` sees states and individual parameters only.** ODE
  intermediates are *not* in scope (`model_parser.rs:6898-6903`) and silently
  become covariate references. `y = ... CT` where `CT` is an ODE intermediate
  validates `ok = TRUE` model-only and fails only with data
  (`E_MISSING_COVARIATE`). Inline into `y`; never reference an intermediate.
  Same restriction on `init()`.
- **Model-only validation is necessary but not sufficient.** An undefined name
  anywhere a covariate is legal is silently accepted. Phase 1 must pass `data=`
  whenever the `$DATA` file resolves, and say so when it does not.
- **Covariate name matching is case-SENSITIVE** (`datareader.rs:690-699` uses
  case-insensitive matching only for the standard columns). We uppercase every
  name via `.norm()`. A dataset column `Flag` against an emitted `FLAG` is an
  `E_MISSING_COVARIATE` at fit time. Covariate references must preserve the
  case as written in `$INPUT`.
- **`[odes]` has no use-before-def check, and the failure is silent and total.**
  Moving `CT = CENT/VC` below the `d/dt` line that uses it keeps the model
  valid and turns `PRED` from `28.11 24.52 19.93 14.82 9.32` into a constant
  `33.33`. Source order must be preserved exactly. (`[individual_parameters]`,
  by contrast, *does* have a forward-reference check -- `collect_forward_refs`,
  `model_parser.rs:1255` -- so ordering errors there are loud.)
- **State / individual-parameter / ODE-intermediate names must be distinct,
  case-insensitively.** An intermediate named `CL` beside an individual
  parameter `CL` gives `E_PARSE: [odes]: name 'CL' collides ...`. A `$DES` block
  that reassigns a `$PK` name is legal NONMEM and illegal ferx -- it needs a
  rename.
- **`init()` may reference individual parameters and states, not thetas.** The
  theta pass-through (defect 2) is a prerequisite for the `A_0` fix (defect 3)
  whenever the baseline is a fixed effect.
- **Unknown function names are accepted and silently evaluate to the identity.**
  The whitelist (`exp`, `log`/`ln`, `sqrt`, `abs`, `floor`, `ceil`, `round`,
  `inv_logit`/`expit`, `logit`) is enforced only at eval. `min()`/`max()` do not
  exist in expressions, so NONMEM models using them are untranslatable and must
  raise an `ERROR` rather than emit a call ferx will quietly ignore.
- **`pk` macro arguments:** `cl`, `v`/`v1`, `q`/`q2`, `v2`, `ka`, `f`,
  `q3`, `v3`, `lagtime`/`alag`, plus `n`, `mtt`, `mat`, `cv2`
  (`types.rs:692-710`). There is **no `dur=`, `rate=`, `d1=`, `r1=`** --
  infusion comes from the data `RATE` column, so `D1`/`R1` in `$PK` is genuinely
  untranslatable. Unknown argument names are a hard parse error, which is the
  good case.
- **`[odes]` alongside a `pk ...(...)` structural block is silently ignored.**
  Never emit both.

---

## 6. Implementation plan

Eight phases, one PR each, in dependency order. Phase 0 is first because it is
the smallest change and the most severe defect. Phase 1 is second because it
converts most of the remainder from silent to loud, in CI, before any of them
are fixed.

Every phase carries `roxygen2::roxygenize()`, the non-ASCII check,
`R CMD check --as-cran`, and a full `devtools::test()` including tier 4.
Baseline to preserve: 271 passing.

### Phase 0 -- De-shadow theta names (defect 14, RC-F) -- DONE

Implemented on branch `worktree-fix+theta-shadowing`. Two things the design
above got wrong, both caught by inspecting output rather than assuming:

- **The rename must be applied to `name_map` before the expressions are parsed,
  never to the emitted strings afterwards.** A textual pass cannot distinguish
  the two meanings: in `cl <- t.CL * exp(eta2)` the RHS means the theta, but in
  the following `k20 <- cl/v` it means the individual parameter, and both
  deparse to the same token once normalised. Renaming after the fact turned
  `K20 = CL/V` into `K20 = TVCL/V`, preserving the bug in a new spelling, and
  also clobbered `pk one_cpt_oral(cl=CL)` into `cl=TVCL`.
- **The collision pre-scan must skip theta-alias self-assignments.**
  `nonmem2rx` emits `tvcl <- t.TVCL`; counting that as a future individual
  parameter renamed a perfectly good `TVCL` to `TVTVCL` and churned six models
  that were never affected. With the alias filter, exactly the two shadowed
  models change and the other eight are byte-identical.

Result: `FAIL 0 | PASS 300` (was 271; tier 4 included), `R CMD check --as-cran`
0 ERROR / 0 WARNING / 2 pre-existing NOTEs. Only `pk_1cmt_oral.mod` and
`pk_1cmt_oral_ampsim.ctl` changed, both snapshot diffs reviewed line by line.
`pk_1cmt_oral.mod` now reads `CL = TVCL * exp(ETA_CL)` / `K20 = CL/V`, and the
`CL` line is live -- the dead-line probe that previously changed nothing now
moves `PRED[1:3]` from `0.0776 0.1228 0.1620` to `0.0361 0.0382 0.0366`.

**Direct confirmation from the fit.** Fitting the old and new
`pk_1cmt_oral.mod` output against the same data:

```
old (shadowed)  theta KA=0.1018 CL=2.054 V=1.021   omega 0.01026, 0.02
new (fixed)     theta TVKA=0.1018 TVCL=2.053 V=1.021   omega 0.01027, 6.14e-06
```

The old omega for `ETA_CL` comes back as exactly `0.02` -- its initial value,
untouched. That is the signature of a parameter with no gradient: the optimiser
never moved it because it had no effect on the likelihood. Structural thetas
barely move, which is why the tier-4 test passed through all of this.

**Two follow-ups this exposed, neither fixed here.**

1. `data-raw/generate_concordance_data.R` substituted theta initials with
   `sub("theta KA[(]...")`, which silently no-ops under the new names --
   it would have simulated from the wrong parameters without saying so. Replaced
   with a `set_theta()` helper that accepts either spelling and `stop()`s on no
   match.
2. `inst/testdata/ode_1cpt_oral_concordance.csv` was itself simulated from the
   shadowed model, so it carries no IIV on clearance, and the corrected model
   correctly recovers `omega_CL` near zero. **The tier-4 test therefore does not
   yet exercise the repaired path.** Regenerating it is blocked on a separate,
   pre-existing breakage: `ferx_simulate()` in the generator returns zero rows
   against `ferx` 0.2.0, failing on the very first dataset
   (`1cpt_oral.ctl`, which this change leaves byte-identical). Until that is
   fixed the guard is the tier-1 test asserting `K20 = CL/V` with de-shadowed
   thetas, plus the corpus-wide invariant in `test-integration.R`.

Original design notes follow.

- After thetas and individual parameters are both known, rename any theta whose
  name collides (case-insensitively) with an individual-parameter LHS:
  `<name>` -> `TV<name>`, then `THETA_<name>`, then a numeric suffix. Rewrite
  every reference. Emit one `INFO` warning per rename.
- Add a tier-1 invariant test asserting no emitted theta name equals any
  emitted individual-parameter name, and run it over the whole bundled corpus.
- **Expect snapshot and concordance churn.** The emitted text of
  `pk_1cmt_oral.mod` and `pk_1cmt_oral_ampsim.ctl` both change. Only
  `pk_1cmt_oral.mod` changes *numerically* -- IIV on `CL` becomes live there for
  the first time -- so its tier-4 reference omegas must be re-baselined in the
  same commit and its structural-theta assertions re-checked, not assumed.
  ampsim should be text-only; if its fit moves, something else is wrong and the
  PR should stop.

### Phase 1 -- Validate the output (defect 7, RC-D) -- DONE

Implemented. `to_ferx(..., validate = TRUE, strict = TRUE)`; the TMDD reproducer
now aborts with the engine's own `E_PARSE` instead of returning a broken file.

Three things worth recording:

- **Validation notes must not reach the emitted file.** Folding the INFO lines
  ("validated without data", "ferx is not installed") into `ir$warnings` made
  every clean model emit `# Warnings: 1` for a note saying nothing is wrong, and
  churned every snapshot. Only engine findings *about the model* go into
  `$warnings` and the file; the notes stay in `$validation`. With that split all
  ten bundled models stay byte-identical.
- **The gap report already existed** -- in `test-concordance.R`, not
  `test-integration.R`. It reported `$unsupported` and then called `succeed()`,
  so it could only ever report what the translator already knew it could not do.
  It now validates each model and fails on engine errors. Proven to fire: adding
  the TMDD reproducer to the corpus turns it red with the exact `E_PARSE`.
- **The CI engine pin was stale.** It pinned `ferx-r@54f25d4` = 0.1.5 while all
  of this was baselined on 0.2.0. Bumped to `ferx-r@731adc9` = 0.2.0, and
  `CLAUDE.md` now records that the bundled `inst/testdata/` datasets are tied to
  the pin as well as the reference omegas.

The ferx-absent branch is behind `.has_ferx()` purely so it can be mocked and
asserted -- it is how the fast PR job runs, so "an optional dependency never
breaks a translation" has to be a test, not an intention.

`FAIL 0 | PASS 332`. `R CMD check --as-cran`: 0 ERROR, 0 WARNING, 2 pre-existing
NOTEs.

Original design notes follow.

### Phase 1 design notes

- Add `ferx` to `Suggests` (already used by `test-concordance.R` via
  `skip_if_not_installed()` but undeclared). Add `amp.sim` too -- `CLAUDE.md`
  documents it as a `Suggests` dependency and it is absent.
- `to_ferx(..., validate = TRUE, strict = TRUE)`:
  - `ferx` not installed -> one `INFO` warning explaining that validation was
    skipped. Never fail for a missing optional dependency.
  - resolve the `$DATA` path relative to the control stream. If it exists, call
    `ferx_model_validate(tmp, data = <path>)`; otherwise validate model-only and
    emit an `INFO` warning that covariate and endpoint coverage were **not**
    checked. This distinction is load-bearing -- see 5.4.
  - map every diagnostic into `ir$warnings` with its `E_`/`W_` code verbatim,
    and every `error`-severity one into `result$unsupported`.
  - `strict = TRUE` (default) aborts on any `error`; `strict = FALSE` downgrades
    to `cli_warn` and returns the result.
- Extend the existing **translation gap report**. It lives in
  `test-concordance.R` (not `test-integration.R`), translates every model in
  `inst/testmodels/nonmem/`, and prints a `model -> gap` table of
  `$unsupported`. Two limits: it never runs the engine, and it always passes
  (`succeed()`), so it reports only what the translator already knows it cannot
  do -- never what it got wrong. Add a validation pass over each emitted model
  and make it *fail* on any `error`-severity diagnostic while still merely
  reporting warnings. That is the regression net for phases 2-7.
- Be explicit in the PR that this does not catch defect 14 (no diagnostic
  exists) -- phase 0 is the only guard there.

### Phase 2 -- Identifier sanitisation (defect 1, RC-B) -- DONE

- One `.ferx_ident()` helper mapping any name to `[A-Za-z_][A-Za-z0-9_]*`,
  applied at a single choke point covering state names, `obs_cmt`,
  `states=[...]`, individual-parameter LHS, and every symbol in every emitted
  expression. Not ad hoc per call site.
- **Collision-aware:** `c.RTOT -> C_RTOT` is safe only if nothing else
  normalises to `C_RTOT`. Keep a rename map, compare case-insensitively (ferx's
  rule), disambiguate with a numeric suffix, and emit one `INFO` per rename.
- Covariate references are the exception: preserve source case (5.4).
- Note the asymmetry -- a dotted name is *accepted* inside `states=[...]` but
  rejected at every reference site. "It parsed" is not evidence.

### Phase 2 design notes

Two things the plan above got wrong, found while building it.

- **Covariate case was already correct, and the plan's framing risked breaking
  it.** `nonmem2rx` preserves the `$INPUT` case for data items while lowercasing
  assigned variables, and a covariate is by definition absent from `name_map`, so
  `.normalise_expr()` leaves the symbol alone and the source spelling survives.
  A "single choke point" applied to every emitted symbol -- which is how phase 2
  was written -- would have uppercased them and turned every covariate into an
  `E_MISSING_COVARIATE` at fit time. Nothing tested it. The correct action was a
  regression test plus an `ERROR` for an illegal covariate name, which cannot be
  renamed at all.
- **Do not reserve theta names when choosing a state name.**
  `.deshadow_theta_names()` is the single owner of theta naming; reserving thetas
  in the state sanitiser too gives one collision two owners, and they rename
  against each other (a state `central` beside a theta labelled `CENTRAL` became
  `central_1` *and* the theta became `TVCENTRAL`). The de-shadow loop already
  reserves state names, so the theta side is where the clash is resolved.

Also worth recording for later phases: renaming only what must be renamed
matters. Sending states through `.norm()` rather than `.ferx_ident()` would
uppercase `depot`/`central` in every ODE model in the corpus -- names users read
and index by -- for no gain.

### Phase 3 -- Theta pass-through, generalised (defect 2) -- DONE

- Replace the `pk_candidates` name whitelist (`rxui_to_ir.R:103`), which covers
  only `linCmt` models and only hardcoded PK names, with a general rule: collect
  every symbol referenced by `[odes]`, `init()`, `[scaling]`, and the pk macro
  args; any that resolves to a declared theta with no individual-parameter
  definition gets a pass-through.
- With phase 0 in place the pass-through reads `CL = TVCL`, not `CL = CL`.
- Required before phase 4: `init()` rejects a bare theta.

#### Phase 3 design notes -- what the plan got wrong

**`[scaling]` does not need pass-throughs.** The bullet above lists it alongside
`[odes]`; section 5.4 says the opposite, and 5.4 is right. Measured against both
ferx 0.2.0 and 0.3.0, a theta is readable from `[individual_parameters]`, from
`[scaling] y` and from `obs_scale`; it is NOT readable from a `d/dt` right-hand
side, an ODE-block intermediate, an `init()` expression, or a pk macro argument.
Only the second group gets a carrier. Treating `[scaling]` as out of scope would
rename thetas and add parameters for references that were already correct.

**The whitelist was not replaced.** `.PK_CANDIDATES` still drives the linCmt
pass-through and should stay: those models have no `[odes]` block for the general
rule to read, and the two mechanisms answer different questions -- "which pk macro
arg is missing" versus "which emitted ODE symbol resolves to nothing".

**Discovery cannot key on the theta's source name.** In any de-shadowed model the
source name IS the individual parameter's name, so `d/dt(ABS) = -KA * ABS` beside
`theta TVKA` and `KA = TVKA * exp(ETA_KA)` reads the parameter. Matching `KA`
invented a second carrier for a correct reference and moved a bundled snapshot.
The rule that works is: the theta's *current emitted* name appears in the emitted
ODE text, OR its source name appears and nothing else declares that symbol. It
must be recomputed each round of the de-shadow fixpoint, not accumulated -- a
theta matches on its source name only until the rename lands.

**A carrier may reuse an existing parameter only if that parameter is a pure
alias of the theta.** `frac <- central/cl` above `cl <- cl*exp(eta.cl)` reads the
theta; pointing it at `CL` substitutes the IIV-applied value. Both forms parse and
both fit, so nothing catches it. Such a reference gets its own carrier.

**The scope has to be applied to the emitted text, not during the walk.** An
ODE-block intermediate that touches a state is inlined, and the inlined text was
normalised for a different context -- `ki <- KTP*CENT` arrived as `TVKTP * CENT`.
Resolving only the `d/dt` line also made the output depend on statement order.
This is the same defect as the phase-2 state map, one layer further in.

**nonmem2rx does not bind a theta referenced by its `$THETA` label.** `FLUX =
KTP*A(1)` for `(0,0.2) ; KTP` leaves a free `KTP` in `lstExpr` and records the
theta in `rxmissingvars1 <- t.KTP`. We were emitting that placeholder as
`RXMISSINGVARS1 = TVKTP` while the reference stayed dangling. Placeholders are now
dropped and the reference bound. `inst/testmodels/nonmem/ode_theta_ref.ctl` is the
fixture. An undeclared-but-legal identifier passes the phase-2 legality check --
that check tests the grammar, not whether anything declares the name -- but ferx
does catch it: `[odes]: RHS references undefined name(s): KTP`, with or without a
dataset, measured against 0.3.0. So the fixture fails the concordance corpus sweep
and aborts `to_ferx(strict = TRUE)` if the carrier regresses.

**Carrier naming.** Source name where free (`KTP = TVKTP`) -- ferx's own examples
do this and it keeps the emitted `[odes]` diffable against `$DES`. Where taken,
derived from the theta's *emitted* name: `TVCL_ODE = TVCL`. Not `CL_ODE`, which
reads as the individual value; not `CL_1`, because `_1` already means "state
disambiguated" and is `.free_theta_name()`'s last resort, and `CL_1`/`V_1` are
plausible model variables. The numbered last resort is `TVCL_ODE_1` and warns,
because the number is positional. Measured: ferx accepts all of these, including a
leading underscore, so this is a readability decision, not a grammar constraint.

Phase 4 caveat: `init()` needs carriers too and lives inside `[odes]`, so `_ODE`
stays accurate there. A pk macro argument does not -- if phase 4 needs a fallback
name on that path, either pick a block-neutral suffix then or accept two.

#### The declaredness check, and the premise that was wrong

Added in this phase: every name the emitted `[odes]` block references must resolve,
or it is an `ERROR` with an `$unsupported` entry.

Two wrong premises died here, in order, and both are worth recording because each
one produced a check that looked finished.

**First: checking against the covariate set is circular.** The plan said "`[odes]`
is the unambiguous block, so check it against states, individual parameters and
covariates". But **the covariate set is defined as the symbols nothing else
binds**. `.covariate_names()` and rxode2's `ui$allCovs` both classify a name the
translator failed to bind as a legitimate covariate -- measured, both call the
unbound `CF` in `qss_tmdd.mod` and the unbound `KTP` in `ode_theta_ref.ctl`
covariates. That version passed the entire test suite and could not fire on either
defect it was written for.

**Second: the fix for that was also wrong.** The conclusion drawn was "the
authority must come from outside the translator", and NONMEM `$INPUT` was parsed
(`.extract_nm_input()`) to supply it, with a `WARN`-level degradation for sources
that have no column list. Then the engine rejected the fixture:

> `[odes]: RHS references undefined name(s): WT. An ODE RHS may only reference
> declared states, individual parameters, ODE-block intermediates, or the reserved
> TIME/TAFD/TAD/MACHEPS variables. If one of these is a covariate, pre-compute the
> covariate-dependent term in [individual_parameters] and reference that variable
> here instead.`

**A covariate is not in scope in `[odes]` either.** The set is closed, so there is
no ambiguous case, no column list is needed, every leftover is an ERROR regardless
of source format, and `.extract_nm_input()` was deleted unused. The earlier scope
probe measured thetas, etas, sigmas and TIME and never measured a covariate --
it tested the cases that were expected to matter.

Two reports, because the remedies differ: a theta/eta/sigma needs a carrier,
anything else needs the term pre-computed one block earlier. The first is what
makes an eta-in-ODE reportable for the first time.

Corpus false-positive rate is zero. Measured on ferx 0.3.0 (ferx-r tag `v0.3.0`),
the version shipped to users; no bundled model references a covariate from
`[odes]`, so nothing that translated before changes either way.

Worth acting on separately: the `engine` CI job still pins `ferx-r@731adc9` =
0.2.0, so it validates the translator against a build no user runs. The full suite
passes against `v0.3.0` unchanged (`FAIL 0 | PASS 518`, engine tier included),
which suggests the bump needs no re-baselining of the concordance references or the
`inst/testdata/` datasets -- CLAUDE.md warns that both are tied to the pin, so that
should be confirmed rather than assumed, and CLAUDE.md updated in the same commit.

For phase 5 and beyond: extending this to `[individual_parameters]` is a much
weaker check (thetas ARE in scope, so the legitimate set is far larger) and needs
the same column list plus per-block scope tables. It should follow phase 5, not
precede it, because phase 5 stops inlining and emits ODE intermediates -- which
changes the `[odes]` declared set (see `odes_intermediates`, already read here and
NULL until then) and makes statement order significant. Also relevant to phase 5:
ferx emits `computed but never used` for a dead individual parameter, which is the
diagnostic its "drop dead intermediates" requirement needs.

Still open, found while doing this and out of scope here: `obs_cmt` is not
inferred from `$MODEL COMP=(X, DEFOBS)` (the guess happens to be right in the
fixture), and an inline scaling division in `$ERROR` (`Y = A(2)/V*(1+EPS(1))`) is
dropped silently -- only `S1`/`S2` assignments are picked up.

### Phase 4 -- Non-symbol assignment targets (defects 3, 9, RC-C)

**Branch phase 4 off `main` after phase 3 lands, NOT off the phase 3 branch.**
Phase 3 is expected to be squash-merged, which is safe only while nothing stacks
on it. Branching phase 4 off it recreates exactly the trap phase 2 hit: the parent
PR's commits never become ancestors of `main`, `merge-base` stays behind, and the
child PR's diff re-shows the parent's work as new (measured on that pair: 12 files
/ 3353 lines instead of 6 / 1158). Phase 4 also modifies functions phase 3
introduces (`.emitted_ode_symbols()`, `.scope_odes_to_params()`, the `[odes]`
declaredness check), so the conflict would be semantic as well as textual.

Replace the silent `next` at `rxui_to_ir.R:432` with explicit dispatch:

| `lstExpr` form | action |
|---|---|
| `STATE(0) <- expr` | `init(STATE) = expr` in `[odes]`, or `[initial_conditions]` for a pk macro |
| `f(STATE) <- var` | `f=` pk macro arg |
| `alag(STATE) <- var` | `alag=` pk macro arg |
| `dur(...)`, `rate(...)` | `ERROR` + `unsupported`; ferx has no such arg (5.4) |
| anything else | `ERROR` + `unsupported`, naming the construct |

The default arm must be loud. A silently unhandled LHS form is how defect 3
shipped in `pkpd_ir.mod`.

Also fix `.infer_pk_macro()` so `nonmem2rx`'s `rxf.depot.` / `rxalag.depot.`
resolve -- match on the `RXF_` / `RXALAG_` prefix rather than a fixed name list.

Regression target: `pkpd_ir.mod` must emit `init(EFFECT) = BL`. Verified valid.

#### Phase 4 measurement pass (done before implementation)

Measured against ferx 0.3.0 and nonmem2rx, because the phase 3 text was wrong
about `[scaling]` and the phase 3 implementation was then wrong twice about what
`[odes]` accepts. Nothing below is inferred from the plan or from ferx-core source
alone.

**The two init paths have DIFFERENT scope.** The table above treats them as
interchangeable ("in `[odes]`, or `[initial_conditions]` for a pk macro"). They are
not, and the difference decides whether phase 3's carrier machinery is needed:

| reference | `init()` in `[odes]` | `[initial_conditions]` (pk macro) |
|---|---|---|
| constant | yes | yes |
| individual parameter | yes | yes |
| expression of individual parameters | yes | yes |
| another state | yes | n/a |
| `TIME` | yes | not measured |
| **theta** | **no** | **yes** |
| **eta** | **no** | **yes** |
| **covariate** | **no** | **yes** (re-checked with the column present in data) |

So `init()` in `[odes]` has exactly the closed set of a `d/dt` right-hand side, and
the phase 3 carrier + declaredness machinery extends to it with **no rule change**
-- `.emitted_ode_symbols()`, `.scope_odes_to_params()` and the `[odes]` check just
need init expressions folded in, which the code comment there already anticipates.
`[initial_conditions]` needs none of it.

The covariate row was taken twice on purpose: without a dataset an unknown name is
read as a covariate and reported valid, so the first reading proved nothing. With
the column present in the data, `[initial_conditions] init(central) = WT` is
genuinely accepted and `[odes] init(CENT) = WT` is genuinely rejected.

**`[initial_conditions]` compartment names are not the model's state names.** It
takes `central`, `depot` (oral models only) or a 1-based number; `init(CENT)` is
rejected with "unknown compartment". The translator has to map, not pass through.

**pk macro parameter names**, from the engine's own error text:
`cl, v/v1, q/q2, v2, ka, f, q3, v3, lagtime/alag`. So `alag=` and `lagtime=` are
aliases -- the plan's `alag=` is valid, and the bundled ferx examples use
`lagtime=`. `dur=` and `rate=` are genuinely absent, as the plan says. `f=` also
scales IV bolus/infusion doses (ferx-core #327), so it is not oral-only.

**A pk macro argument must be a single declared individual parameter.** `f=TVB` (a
theta) and `f=BL*2` (an expression) are both rejected. The plan does not mention
this: it makes phase 3's carriers a prerequisite on the pk path too, not only for
`[odes]`.

**nonmem2rx always emits these as a PAIR** -- a helper assignment with a dotted
name, plus the non-symbol-LHS statement that reads it:

```
A_0(4)=BL     ->  rxini.rxddta4. <- bl        +  EFFECT(0) <- rxini.rxddta4.
F1 = THETA(3) ->  rxf.central.   <- theta3    +  f(central) <- rxf.central.
ALAG1 = ...   ->  rxalag.central. <- theta4   +  alag(central) <- rxalag.central.
D1 = ...      ->  rxdur.central.  <- theta5   +  dur(central)  <- rxdur.central.
R1 = ...      ->  rxrate.central. <- theta6   +  rate(central) <- rxrate.central.
```

The regression target is reachable, but not directly: `EFFECT(0)` reads
`rxini.rxddta4.`, not `BL`, so emitting `init(EFFECT) = BL` means resolving that
alias chain. The plan says "Verified valid" without mentioning it. The dotted
helper names are a third instance of the `rxmissingvars` shape already handled in
phase 3 -- they normalise to `RXINI_RXDDTA4_`, `RXF_CENTRAL_` and so on, and must
be consumed rather than emitted.

**What ships today, both silently.** Confirmed by translating:

- `pkpd_ir.mod`: no `init(...)` at all, so the effect compartment starts at 0
  instead of `BL`. `$unsupported` is EMPTY. (The consequence for the fit was not
  measured -- what is measured is that the initial condition is absent from the
  emitted model and nothing reports it.)
- A model with `F1`/`ALAG1`/`D1`/`R1`: the pk macro emits `one_cpt_iv(cl=CL, v=V)`
  with no `f=` and no `lagtime=` **although ferx supports both**;
  `RXALAG_CENTRAL_`, `RXDUR_CENTRAL_` and `RXRATE_CENTRAL_` are emitted as
  meaningless individual parameters; `RXF_CENTRAL_` disappears entirely.
  `$unsupported` is EMPTY.

Note that phase 3's `[odes]` declaredness check does **not** catch either of these.
It reports dangling references; these are silent *omissions*, which leave nothing
dangling. Only the loud default arm catches them, which is the point of this phase.

### Phase 5 -- Ordered statements and `if` support (defects 4, 6, 8, 13, RC-A)

The largest phase. Two coupled changes.

**5a. IR carries ordered statement lists.** `indiv_params` and `odes` each
become an ordered list of tagged statements: `assign(lhs, rhs)`,
`if(cond, then, else)`, `init(state, expr)`, `ddt(state, rhs)`. Keep accepting
the current flat `list(lhs, rhs)` shape as an `assign` so hand-built IRs in
`test-emit.R` and `test-ir.R` keep working. `emit_ferx()` gains a statement
renderer for both blocks.

**5b. Stop inlining `$DES` intermediates; emit them as ODE-block
intermediates.** Delete `.inline_aux_vars()`. Reasons in order:
- it cannot represent a variable defined inside an `if` -- which is defect 4;
- its depth-30 cutoff silently returns the un-inlined expression (defect 13),
  leaving an undefined name in the output;
- repeated substitution is exponential in nesting depth (`pkpd_ir.mod` already
  duplicates `CENTRAL/V2`);
- ferx supports intermediates directly and the output becomes diffable against
  the source `$DES`.

Constraints, all verified in 5.4:
- **Preserve source order exactly.** There is no use-before-def check and the
  failure is silent and total (`PRED` collapses to a constant). `pkpd_ir.mod`
  already needs this -- `C2`/`EFF` sit between `DADT(3)` and `DADT(4)`.
- **Rename on collision** with a state or individual-parameter name, with a
  warning.
- **Partition `$DES` from `$ERROR`.** Both reference states, so a variable is an
  ODE intermediate iff it is transitively reachable from a `d/dt` RHS or an
  `init()`; otherwise it is readout/error scaffolding and belongs to phase 6 --
  not `[odes]`, not `[individual_parameters]`. This is what fixes defect 6
  (`W1`/`W2`) and the misplaced `FB`.
- **Drop dead intermediates.** The same reachability rule must delete unused
  ones. `pk_1cmt_oral.mod` has an unused `CP = A(2)/V` in `$DES` that inlining
  currently discards for free; emitting it would produce
  `W_PARSE: computed but never used` and trip the phase-1 gap report.
- Emit `if` statements verbatim into whichever block their target belongs to.

Regression target: defect 8's model must emit `if (SEX == 1) { ... }`.

#### Phase 5 measurement pass (done before implementation)

Measured against ferx 0.3.0 with the CLI at `maxiter = 0`, because the phase
text asserts two things the engine has to actually do and the plan has been
wrong three times already about what `[odes]` accepts.

**`if`/`else` inside `[odes]` parses AND branches.** This is the linchpin: 5.4
records that unknown function names are accepted and silently evaluate to the
identity, so "it validates" would not have been evidence. The discriminating
test is a conditional whose two arms give different answers, run with the
condition forced each way against a control with the branch hardcoded:

```
[odes]
  CT = CENTRAL / V
  if (<cond>) { SCL = 1.0 } else { SCL = 0.25 }
  d/dt(CENTRAL) = -KE * CENTRAL * SCL

TIME      cond true (SCL=1)  cond false (SCL=0.25)  no `if`, hardcoded 1.0
 0.5           1.902459             1.975156              1.902459
 2.0           1.637462             1.902459              1.637462
24.0           0.181450             1.097625              0.181450
```

The true arm equals the hardcoded control to every printed digit and the false
arm does not, so both arms are really evaluated and the condition really
selects. Braces on both arms; `else` on the same line as the closing brace.

**ODE-block intermediates in source order validate alongside the `if`.** The
same model declares `CT`, `BB`, `FB` above the `d/dt` line that reads them and
reports `ok -- no errors (0 warning(s))`. Note this is NOT independent evidence
that the order is right: 5.4 records that `[odes]` has no use-before-def check
and that the wrong order stays valid while collapsing `PRED` to a constant. The
order constraint has to be honoured by construction, not confirmed by the
validator.

**Consequence for 5b.** Both halves of the phase-5 target output are expressible
today, so nothing here needs a ferx-core change.

#### Phase 5b design notes -- what the plan above gets wrong

Written after reading `.parse_model_exprs()` with 5a landed, before implementing
5b. Five corrections; the first two change the rule itself, not its wording.

**1. The intermediate rule is a CONJUNCTION, not the backward reachability
above.** The text says an ODE intermediate is a variable "transitively reachable
from a `d/dt` RHS or an `init()`". That alone is wrong: a `$PK` assignment is
routinely reachable from a `d/dt`, and it must stay an individual parameter --
individual parameters ARE in `[odes]` scope, so nothing is gained by moving it,
and moving it loses the per-subject evaluation. The real test is

    ODE intermediate  ==  references a state (transitively)
                          AND is reachable backward from a d/dt RHS or init()

The first half already exists and is computed: it is `aux_vars`, the pass-2
fixpoint. Only the second half is new. The distinction is exactly what separates
`FB` (reads `CF`, which reads `CT = A(1)/VC`, which reads a state -> must move)
from `KSYN = RBASE*KDEG` (reads no state -> stays an individual parameter).

**2. Defect 6 is NOT covered by that rule, in either direction.** `W1 = 0`
references nothing, so it is not in `aux_vars`; and nothing in `[odes]`
references it, so it is not backward-reachable either. It reaches
`[individual_parameters]` through pass 3's DEFAULT -- "everything not in
`aux_vars` becomes an individual parameter" -- which no reachability rule
touches. Evicting it needs its own rule, and the narrow one is safest:

    drop an assignment that is referenced by NOTHING emitted
    AND is referenced by the error/readout expression

That is precisely `W1`/`W2` and does not disturb anything else. The broader
"drop every unused individual parameter" is tempting and should be resisted in
this phase: it changes output for models unrelated to this issue, and ferx
already reports those as `W_UNUSED_PARAM`.

**3. Source order between an intermediate and a `d/dt` line is currently
UNRECOVERABLE.** `all_assigns` and `odes` are separate lists built in the same
pass, and neither records position. Since order in `[odes]` is a correctness
property with a silent failure (5.4), 5b has to add a `pos` counter to both in
pass 1 and interleave on it. This is the cheapest of the five and the easiest to
forget, because the emitted file looks right either way until a model puts a
`$DES` intermediate after a `DADT` line -- `pkpd_ir.mod` already does.

**4. There is no `if` branch in pass 1 to upgrade -- an `if` falls off the end
of the loop.** The statement matches none of the branches (`cmt`,
`rxmissingvars`, linCmt tilde, tilde, assignment), so it is discarded with no
diagnostic. `.flatten_stmts()` DOES walk into `if` bodies but is only called
from the covariate census, not from the parse. So this is a new branch, not a
changed one, and defect 8 ("`IF` in `$PK` dropped, zero warnings") is the same
missing branch seen from `$PK`.

**5. An `if` cannot be half-emitted, so reachability granularity is the whole
statement.** If any name assigned in EITHER arm is needed, the entire `if` is
emitted. Treating the arms independently would emit a conditional that assigns a
name in one branch and not the other, which is a different model. Nested
assignments therefore have to register their LHS in `name_map` for later
references to resolve, while the emission decision stays at the `if` level.

**Consequence for sequencing.** None of this changes the phase 5 -> phase 6
dependency; it sharpens it. The partition in correction 2 is what hands phase 6
the `W1`/`W2` statements as endpoint-selection input, so phase 6 consumes a
THIRD bucket that phase 5b must produce -- not merely the two the plan names.

### Phase 6 -- Error model and readout (defects 5, 10, 11, 12, 15, RC-E)

**6a. Structural classification.** Decompose `Y = ...` into terms and classify
each as additive (`+ EPS`) or proportional (`* EPS` against the prediction).
Emit `combined(prop_sigma, add_sigma)` in ferx's fixed order, not traversal
order. Fixes 10 and 11.

**6b. Endpoint dispatch.** Fold the sequential assignments feeding `Y` into a
single expression, converting each conditional reassignment into a nested inline
`if`: `ipred <- ctot; if (FLAG==2) ipred <- rtot` becomes
`if (FLAG == 2) RTOT else CTOT`. That expression is the `[scaling]` `y` readout.
Then:
- dispatch column is `CMT` -> `y[CMT=n]` + `CMT=n: DV ~ ...`
- any other column -> Form C `y = ...` + covariate-selected `[error_model]`
  if/else. No data change required.
- single endpoint -> unchanged from today.

**6c. `obs_cmt`.** When a `y` readout is emitted, drop `obs_cmt` entirely --
`ode(states=[...])` alone is valid *only* in that case, and the guess
disappears. When no readout is derivable, parse `DEFOBS` from the raw `$MODEL`
block instead of `tail(state_names, 1)`. Only if both fail should a guess
happen, at `ERROR` level.

**6d. Never emit `obs_scale` and `y` together** (defect 15). Once a `y`
expression exists, scaling is part of it (`CENT/VC`); the separate `obs_scale`
path survives only for the no-readout case. Add a tier-1 invariant test -- the
double-scaling is silent and a factor-of-`V` error in every prediction.

**6e. Scope discipline for `y`.** Every symbol must resolve to a state, an
individual parameter, a theta, an eta, or a real data column. Inline down to
those; never reference an ODE intermediate (5.4).

**6f. No fallback guess.** `.parse_error_rhs()`'s
`"complex $ERROR -- classified as proportional, verify"` becomes an `ERROR`
plus an `unsupported` entry, which phase 1 then blocks on. A guess that looks
like a translation is the failure mode this whole issue is about.

### Phase 7 -- `inits = c("control", "final")`

`nonmem2rx(..., updateFinal = FALSE)` is the knob; its roxygen reads "Update the
parsed model with the model estimates from the '.lst' output file." Add
`inits = c("control", "final")` to `to_ferx()`, default `"control"`. Report the
choice as a `WARN` when `"final"` is selected and a `.lst`/`.ext` was found, as
the issue asks.

**Verify before implementing.** I could not construct a valid `.ext` by hand and
had no real NONMEM run directory available, so the mapping from `updateFinal` to
the reported `i change initial estimate of theta1 to 0.029065` message is
inferred from the signature and roxygen, not observed. Confirm against a real
run directory first.

---

## 7. Test plan

Per `CLAUDE.md`'s four tiers, written with each phase, not at the end.

**Tier 1 (unit).** Theta/parameter de-shadowing invariant; `.ferx_ident()`
collision cases; statement-list emission for `if` / `init` / intermediates;
error-term decomposition including the swapped-`combined` and `F + F*EPS(1)`
cases; readout folding of conditional reassignment; the
"never `obs_scale` and `y` together" invariant.

**Tier 2 (integration).** Add `inst/testmodels/nonmem/qss_tmdd.mod` (the
reporter's model -- self-contained, round numbers, no proprietary data) plus
three narrower models isolating single root causes: `IF` in `$PK`, `F1`/`ALAG1`,
and an `A_0` baseline. Snapshot each.

**Tier 3 (reference).** Reference `.ferx` for the TMDD model. Phase 0 changes
`pk_1cmt_oral.mod` and `pk_1cmt_oral_ampsim.ctl`; phase 5b changes
`pkpd_ir.mod`. `_snaps/integration.md` needs review -- `CLAUDE.md`'s snapshot
rule applies; review the diff, do not bulk-accept.

**Tier 4 (concordance).** The only tier that catches defects 3, 8, 10, 11, 14,
15 -- they all produce output that parses.
- Re-baseline `pk_1cmt_oral.mod` and the amp.sim benchmark after phase 0, and
  state in the PR that the previous references were fitted with dead IIV.
- Add `pkpd_ir.mod` against simulated data, asserting the effect-compartment
  baseline is recovered. Fails today with the dropped `A_0(4)=BL`.
- Add a `combined` error model with the additive term written first, asserting
  the two sigmas come back the right way round.
- Add an assertion that a model with an ETA on a `$PK` parameter used only via a
  derived parameter (the `K20 = CL/V` shape) recovers a non-zero omega. This is
  the direct regression test for defect 14.

**New standing check.** The phase-1 gap report should make a green CI mean "the
engine accepted every bundled model", the role `CLAUDE.md` assigns to the
`engine` job for concordance.

---

## 8. Risks and open questions

1. **Phase 0 changes numbers, not just text.** Fits that previously ran without
   IIV on a shadowed parameter will move. Re-baseline deliberately, and check
   whether any structural-theta assertion was passing *because* the IIV was
   dead.
2. **`ferx` in `Suggests` and CI.** The gap report needs the engine. Run it in
   the existing `engine` job, not the fast PR job, and skip cleanly elsewhere.
3. **Snapshot churn** across phases 0 and 5b. Deliberate and correct, but review
   the diffs; do not `snapshot_accept()` in bulk.
4. **`$ERROR` idiom coverage is a catalogue, not an algorithm.** Phase 6 handles
   the common idioms; everything else must reach the user as an `ERROR`. Resist
   widening the fallback -- that is how we got here.
5. **`nonmem2rx` version coupling.** Phases 2, 4, 5 depend on `lstExpr` shapes
   (`rxini.rxddta3.`, `c.` prefixing, `rxf.depot.`). Assert the `nonmem2rx`
   version and add a tier-1 test that fails loudly if those shapes change,
   rather than silently dropping features again.
6. **Local `ferx-r` checkout is stale** -- `../ferx-r` is on
   `feat/sim-horizon-526` at 0.1.6 while installed `ferx` is 0.2.0, and its
   `Cargo.lock` pins a `ferx-core` 483 commits behind local `main`. Read
   `../ferx-core` for format truth, not `../ferx-r`.
7. **Two upstream reports worth filing**, both silent-wrong in ferx itself and
   both hit by this work:
   - `obs_scale` stacking with a Form C `y` (5.4) -- no guard, no warning, and
     the code comment claims the opposite. A hard error would be right.
   - a theta shadowing an individual parameter (3.14) produces no diagnostic.
     `W_PARSE: computed but never used` would make it self-correcting for every
     downstream tool, not just this one.
   Cosmetic third: `ferx-r/R/model.R:1103` lists `"initial_values"` in
   `optional_sections`; the real block is `initial_conditions`, so phase 4's
   output will print as `[unknown section]`.

---

## 9. Explicitly out of scope

- Any ferx-core or ferx-r format change. None is required (section 5). The two
  upstream reports in risk 7 are diagnostics, not format changes, and they are
  not prerequisites.
- Rewriting the user's dataset. `FLAG` dispatch needs no remap.
- `mrgsolve` support (v0.2).
- General symbolic algebra over `$ERROR`. Catalogue plus loud failure, per
  risk 4.

---

## 10. Reply to the reporter

- All seven confirmed. Seven more found, six silent-wrong, including one
  (theta shadowing) that is worse than anything reported and is live in our own
  concordance suite.
- No ferx-core change needed -- their model is fully expressible today.
- Their suggested fix for defect 5 is superseded: ferx dispatches on `FLAG`
  directly, so no data remapping is required.
- Take up their offer to test full TMDD, Michaelis-Menten, Wagner, and
  irreversible-binding structures. Those are exactly the phase 5/6 stress cases,
  and the current bundled corpus contains no model with a conditional in it at
  all.
