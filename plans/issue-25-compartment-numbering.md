# Issue #25 -- compartment numbering: ferx counts states, NONMEM counts `$MODEL`

Plan for [#25](https://github.com/FeRx-NLME/translator/issues/25). Written
against `main` at `80d4400`, measured with `ferx` 0.3.0 (`ferx-r@7889719`, the
SHA the `engine` CI job pins and the version users install).

## 1. Verdict

The issue is real, reproduces exactly as filed, and **understates itself**. It
describes an `obs_scale` bound to the wrong variable. The same root cause also
sends the **dose** to the wrong compartment, which the issue does not mention
and which is the worse half by a wide margin.

The issue also names the wrong trigger as the only one. It says a reordered
`$DES` is needed. A `$MODEL` that declares a compartment with no `DADT` -- an
ordinary NONMEM construct, in-order `$DES`, nothing unusual anywhere -- produces
the same divergence, and there the measured result is **every prediction exactly
zero**.

The fix is smaller than the defect: one emitted line, no ferx-core dependency.

## 2. What was measured

Environment, stated because a stale one produced a wrong answer twice in this
investigation:

- `ferxtranslate` at `main` `80d4400`, clean tree, 0 ahead / 0 behind
  `origin/main`, loaded with `pkgload::load_all()`. The **installed** package
  was pre-#24 and gave pre-#24 answers; measurements taken against it were
  discarded.
- `ferx` 0.3.0 (`ferx-r@7889719`, tag `v0.3.0`) -- the `engine` job's pin.
- `ferx-core` checked at `origin/main` `e7517a01`. The local checkout was 9
  commits stale. Those 9 commits are mixture models and zero-omega; the diff
  touches `types.rs` and `model_parser.rs` with **no** `cmt` / `states` /
  `compartment` line outside comments, and does not touch `predictions.rs`,
  `datareader.rs`, or `dosing.rs` at all. The semantics below are current.

### 2.1 The issue's own table reproduces

Three models, identical but for `$PK` scaling and `$DES` statement order.
`$ERROR` names `A(2)` outright in all three, so `obs_cmt` resolves at the top of
the cascade and is **correct every time** -- confirming this is not an `obs_cmt`
defect:

| `$PK` | `$DES` order | `obs_cmt` | emitted `obs_scale` | what the user is told |
|---|---|---|---|---|
| `S2 = V` | `DADT(1)`, `DADT(2)` | `CENTRAL` ok | `V` ok | `INFO \| S2 = V detected` |
| `S2 = V` | `DADT(2)`, `DADT(1)` | `CENTRAL` ok | *none emitted* | `WARN` naming the contradiction |
| `S1 = VD`, `S2 = V` | `DADT(2)`, `DADT(1)` | `CENTRAL` ok | **`VD`** WRONG | `INFO \| S1 = VD detected` |

Row 2's `WARN` reads `... not for the observed compartment 1 ('CENTRAL')` while
`[structural_model]` says `obs_cmt=CENTRAL`, which is NONMEM compartment 2. The
message contradicts the file it is describing.

### 2.2 ferx numbers compartments by `states=[...]` POSITION

Read in `ferx-core`: a dose carries `cmt: usize`, `cmt_idx()` is
`cmt.saturating_sub(1)`, and `src/ode/predictions.rs` applies it as
`u[cmt_idx] += f_bio * dose.amt` straight into the ODE state vector. So the
data's `CMT` column is a 1-based index into the state list.

Measured rather than left at that. Two `.ferx` files, identical equations,
identical parameters, identical data, differing only in state order, with the
dose row's `CMT` swapped between runs:

```
                                              t=1     t=4    t=12    t=24
inord / dose CMT=1 (DEPOT positionally)   1.34855 1.63868 1.02463 0.49868
inord / dose CMT=2 (CENTRAL positionally) 1.88351 1.57319 0.97339 0.47376
rev   / dose CMT=1 (CENTRAL positionally) 1.88351 1.57319 0.97339 0.47376
rev   / dose CMT=2 (DEPOT positionally)   1.34855 1.63868 1.02463 0.49868
```

Identical to every printed digit across the diagonal. The `CMT` column indexes
`states=[...]` by position, full stop.

### 2.3 `states=[...]` alone controls it -- `[odes]` statement order does not

This is the finding the whole fix rests on, so it was measured directly:

```
                                        t=1     t=4    t=12    t=24
inord odes + states=[DEPOT,CENTRAL] 1.34855 1.63868 1.02463 0.49868
rev   odes + states=[CENTRAL,DEPOT] 1.88351 1.57319 0.97339 0.47376
rev   odes + states=[DEPOT,CENTRAL] 1.34855 1.63868 1.02463 0.49868
```

Row 3 -- `[odes]` statements in one order, `states=[...]` in the other --
reproduces row 1 exactly. `d/dt(NAME)` binds by name; `states=[...]` sets the
numbering. `ferx_model_validate()` on that file: `VALID`, **zero diagnostics**.

That matters because `[odes]` statement order is a correctness property this
package already protects. From `.emit_odes_section()`:

> The rest of the block is emitted in list order and NOTHING here reorders it.
> `[odes]` has no use-before-def check: an intermediate placed below the d/dt
> line that reads it stays valid and reads a stale slot, and PRED collapses to a
> constant with no diagnostic. Source order is the correctness property.

So the naive fix -- reorder the `[odes]` block to match `$MODEL` -- would trade
this defect for that one. It is not available, and `states=[...]` makes it
unnecessary.

### 2.4 The realistic trigger is a `$MODEL` gap, not a reordered `$DES`

An ordinary control stream. In-order `$DES`. One declared compartment with no
`DADT`:

```
$MODEL
  COMP=(DEPOT, DEFDOSE)
  COMP=(DUMMY)
  COMP=(CENTRAL, DEFOBS)
$DES
  DADT(1) = -KA*A(1)
  DADT(3) =  KA*A(1) - (CL/V)*A(3)
```

nonmem2rx materialises the gap as `d/dt(DUMMY) = 0` and places it **first**, so
`nm_to_ferx()` emits `states=[DUMMY, DEPOT, CENTRAL]` against
`COMP=[DEPOT, DUMMY, CENTRAL]`. NONMEM's `CMT=1` dose therefore lands in DUMMY,
whose derivative is zero. Simulated from the translator's own emitted file:

```
                                     t=1      t=4     t=12    t=24
emitted as-is                    0.00000 0.000000 0.000000 0.00000
only states=[] put in COMP order 1.34855 1.638677 1.024627 0.49868
```

The dose vanishes. The file validates clean and says nothing. Changing that one
line -- and only that line -- restores it completely.

This is the case that decides the priority. A reordered `$DES` is legal but
rare; a `$MODEL` compartment without a `DADT` is routine.

### 2.5 Corpus census -- nothing bundled can show this

Every ODE model in `inst/testmodels/nonmem/`, `$MODEL` COMP order against
emitted `states=[...]`:

| model | COMP | states | diverges? |
|---|---|---|---|
| `defobs_expression_wins.ctl` | DEPOT,CENT | DEPOT,CENT | no |
| `defobs_not_last.ctl` | CENT,PERIPH | CENT,PERIPH | no |
| `dotted_state.ctl` | CENT,RTOT | CENT,c_RTOT | no (rename, same position) |
| `ode_theta_ref.ctl` | CENT,PERIPH | CENT,PERIPH | no |
| `ode_warfarin.ctl` | DEPOT,CENTRAL | DEPOT,CENTRAL | no |
| `pk_1cmt_oral.mod` | ABS,CENTRAL | ABS,CENTRAL | no |
| `pkpd_cmt.mod` | CENT,PD | CENT,PD | no |
| `pkpd_ir.mod` | DOSE,CENTRAL,PERIPH,EFFECT | same | no |
| `qss_tmdd.mod` | CENT,TISS,RTOT | CENT,TISS,c_RTOT | no (rename) |
| `s_scaling_not_last.ctl` | DEPOT,CENTRAL,PERIPH | same | no |

Zero divergences. Two apparent mismatches are `.same_cmt_name()` renames at the
same position. Consequences, both of which the plan must act on:

- **No snapshot should change.** If `_snaps/integration.md` moves, the fix did
  something it was not asked to.
- **No existing test can show the fix works.** New fixtures are mandatory, and
  by the discriminating-fixture rule each must be shown to fail first.

### 2.6 What is NOT affected -- checked, not assumed

- `[initial_conditions]`: `A_0(2) = 5` on a reversed `$DES` emitted
  `init(CENTRAL) = 5`. Resolved by name. Correct.
- `[scaling] y[CMT=n]` and `[error_model] CMT=n:`: the `n` comes from the
  source's own `CMT.EQ.n` literals, not from state positions. `pkpd_cmt.mod`
  with a reordered `$DES` emitted the same `y[CMT=1] = CENT/VC`,
  `y[CMT=2] = PD` as shipped. The numbers are already NONMEM-correct -- but see
  4.2, they are correct against a numbering the rest of the file gets wrong.
- nlmixr2 / rxode2 / Monolix sources: `obs_hint` is `NULL` unless
  `format == "nonmem"` and the file is readable (`R/translate.R:82`). No COMP
  list, no reordering, d/dt order remains -- which is right, because those
  languages have no separate compartment numbering to disagree with.

## 3. Root cause

One sentence: **`$MODEL` COMP order is NONMEM's compartment numbering, d/dt
statement order is ferx's, the translator emits the second while the source's
data uses the first, and nothing reconciles them.**

Every symptom is a consequence:

| site | what it does | why it is wrong |
|---|---|---|
| `rxui_to_ir.R:677` `structural$states <- state_names` | emits `states=[...]` in d/dt order | sets ferx's numbering to the wrong permutation -- **the dose** |
| `rxui_to_ir.R:531` `obs_cmt_num <- match(explicit, state_raw)` | a d/dt position | read below as a NONMEM compartment number -- **`obs_scale`** |
| `rxui_to_ir.R:690` `which(state_names_uc == toupper(obs_cmt))[1L]` | another d/dt position | same, on the fallback path |
| `rxui_to_ir.R:733` the no-scaling `WARN` | prints `obs_idx` | names a position as a compartment number |
| `.assemble_endpoints()` `rest <- setdiff(seq_along(state_names), cmts)` | positions | emitted as `CMT=n` keys alongside source-derived ones |

`#24` fixed exactly this confusion in the DEFOBS and `S<n>` tiers, by
cross-checking the ordinal against the name before trusting it. Those two tiers
**decline** on disagreement. This plan does the same thing one layer up, except
that here there is a real permutation available, so it can reconcile rather than
decline.

## 4. Design

### 4.1 One canonical compartment order, computed once

Build the NONMEM-ordered state vector immediately after `state_names` is
available, and use it everywhere a compartment NUMBER is meant:

```
cmt_order <- .nm_cmt_order(state_raw, state_names, obs_hint$comps)
```

- For each `comps[i]`, find the unique `state_raw[j]` with
  `.same_cmt_name(state_raw[j], comps[i])`. Measured total and well-defined on
  both divergent fixtures (`perm=[2,1,3]`, `perm=[2,1]`).
- Returns `NULL` -- meaning "keep d/dt order" -- when the permutation is not a
  bijection: no `comps` at all (non-NONMEM source), a COMP that matches no
  state, a state that matches no COMP, a name matching more than one, or
  differing lengths. Declining is the phase 6c policy and it is the right one
  here: a partial permutation is a guess, and guessing is the defect.
- When it is a bijection AND differs from d/dt order, emit an `INFO` saying the
  state list was put in `$MODEL` compartment order so the source's `CMT` values
  still select the compartments the source meant. A silent correctness change is
  exactly what the warning system exists to prevent.
- When it is NOT a bijection and the two orders would have differed, emit a
  `WARN` plus an `unsupported` entry: the compartment numbering could not be
  reconciled, so a NONMEM dataset's `CMT` column may not mean what the file
  says. That is a user action item, which is what `unsupported` is for.

Then:

- `structural$states <- cmt_order %||% state_names`
- `obs_cmt_num <- match(obs_cmt, cmt_order)` when `cmt_order` exists -- a real
  NONMEM compartment number, which is what `scaling_hint`'s keys are. This
  replaces both `:531` and the `:690` fallback with one resolution against one
  list.
- `.assemble_endpoints()` receives the same list, so `rest` produces ordinals
  rather than positions.

`[odes]` statement order is **left alone**, deliberately (see 2.3). The emitted
`states=[...]` will therefore sometimes list states in a different order from
the `d/dt` lines below it. Measured VALID with zero engine diagnostics. The
readability cost is real and is accepted: `[odes]` order carries a correctness
property that `states=[...]` does not, and trading a silent wrong dose for a
silent stale intermediate is not a trade.

### 4.2 The per-CMT dispatch needs NOTHING -- verified, not assumed

`y[CMT=n]` / `CMT=n:` numbers are source literals and already NONMEM-correct
(2.6). But they were correct against a file whose *numbering* was wrong, so on a
divergent model they were being read against the wrong compartments. Fixing the
numbering fixes them -- for free, and that is the point of doing it at the root.

An earlier draft of this plan claimed `rest` (the fall-through compartment set
in `.assemble_endpoints()`) needed the canonical list threaded into the
expression pass. Checked line by line: `.assemble_endpoints()` and
`.build_endpoints()` use their `state_names` argument only for
`length()`/bounds (`cmts > length(state_names)`, `seq_along`) and for
name-set membership (`setdiff(aux_vars, toupper(state_names))`). Never order.
`rest` is pure integers whose meaning becomes correct the moment the emitted
`states=[...]` numbering is correct. No plumbing. The expression pass is not
touched.

### 4.3 The #24 decline branches must resolve, not decline, when cmt_order exists

Phase 6c gave the DEFOBS tier and the S<n> tier a name cross-check that
DECLINES on d/dt-vs-COMP disagreement, because at the time there was no
reconciliation available -- declining was the only honest move. This fix builds
the reconciliation. Keeping the declines as-is would leave the translator
refusing to infer `obs_cmt` from evidence it now knows how to read, and worse:
their WARN text reasons from a disagreement that the emitted file no longer
has once `states=[...]` follows COMP order. A user reading "the two orderings
disagree, so DEFOBS was not used" next to a file whose numbering IS the COMP
order is being told something false about the output.

So: when `cmt_order` is non-NULL, both tiers resolve through it
(`obs_cmt <- cmt_order[[obs_hint$index]]`, same for `s_idx`) and their
disagreement branches become unreachable by construction. When `cmt_order` is
NULL, the existing cross-checks and declines stay, verbatim -- they are then
still the only honest move. The #24 tests asserting the declines are updated
to pin the NULL-cmt_order path (via a fixture whose COMP list cannot be
reconciled), and the disagreement fixtures flip to asserting correct
resolution. That is deliberate, visible churn: the behaviour is strictly
better, and the sabotage re-run in section 5 is what proves the re-pointed
tests still discriminate.

### 4.4 Implementation traps, named so they are not re-discovered

- The permutation is computed on RAW d/dt names (`.state_raw_names()`, e.g.
  `c.RTOT`) matched against COMP names via `.same_cmt_name()`, but APPLIED to
  the SANITISED list (`state_names`, e.g. `c_RTOT`). Mixing the two lists up
  compiles fine and reorders garbage; the dotted_state/qss_tmdd unit test in
  section 5 exists to catch exactly this.
- The `:690` fallback is case-insensitive today
  (`which(state_names_uc == toupper(...))`). Its replacement must stay
  case-insensitive -- resolve through `.same_cmt_name()`, not bare `match()`,
  or a lowercased nonmem2rx state name silently stops matching.
- `obs_cmt_num` values set by the DEFOBS and S<n> tiers are ALREADY NONMEM
  numbers; only the explicit tier (`:531`) and the `:690` fallback produce
  positions. Do not "fix" the two that are right.

### 4.5 What this deliberately does not do

- Does not touch `[odes]` ordering.
- Does not renumber anything for non-NONMEM sources.
- Does not attempt the trailing-`$MODEL`-gap case in 6.1 -- separate defect,
  separate issue.
- Requires no ferx-core change. `states=[...]` in a non-`[odes]` order is
  already accepted, measured clean.

## 5. Test plan

Every fixture below must be shown to FAIL against unfixed code before it counts.
The corpus census (2.5) proves no existing test can distinguish fixed from
unfixed, so this is not a formality here -- it is the only evidence the tests
are worth anything.

**Tier 1 (unit)** -- `.nm_cmt_order()` directly:
- bijection with a permutation -> the permuted list
- bijection in the identity order -> unchanged, and no `INFO`
- COMP name matching no state -> `NULL`
- state matching no COMP -> `NULL`
- length mismatch -> `NULL`
- renamed state (`c.RTOT` vs `RTOT`) still matches -- guards the
  `.same_cmt_name()` dependency
- `comps = NULL` -> `NULL`, no warning

**Tier 2 (integration)** -- new bundled fixtures. Two are required, because one
cannot show both halves:
- `cmt_order_gap.ctl` -- a middle `$MODEL` compartment with no `DADT`,
  in-order `$DES`. The realistic trigger (2.4). Asserts `states=[...]` in COMP
  order and the `INFO`.
- `cmt_order_reordered_des.ctl` -- reordered `$DES`, **two** `S<n>` entries.
  Two is not optional: with a single `S2` the defect drops the scaling loudly
  instead of substituting the wrong variable, so a one-scaling fixture cannot
  tell "wrong variable" from "no variable". Asserts `obs_scale = V`, not `VD`.
- a fixture whose COMP list cannot be reconciled -> `WARN` + `unsupported`, and
  d/dt order retained.
- Both new fixtures carry a DO-NOT-REORDER comment, as `ode_theta_ref.ctl` and
  `s_scaling_not_last.ctl` do.

**Tier 3 (snapshot)** -- assert the *absence* of churn. Per 2.5 no existing
snapshot should move. If one does, stop and find out why before accepting.

**Tier 4 (concordance)** -- the tier that would have caught this. The 2.4 case
predicts exactly zero, so a fit on it is unmissable. Add a concordance test that
fits the gap fixture against pre-simulated data and recovers the true thetas;
without the fix, IPRED is identically 0 and it cannot pass by accident. This is
the strongest discriminating fixture available and is worth the ~2 minutes.

**Regression risk to check explicitly.** Section 4.3 makes the phase 6c
disagreement branches unreachable whenever `cmt_order` exists, by design. The
#24 disagreement fixtures flip to asserting correct resolution, and a new
irreconcilable-COMP fixture pins the decline path. Re-run the #24 sabotage
mutations against the re-pointed tests: if any mutation stops failing, the
test has gone vacuous and must be re-aimed. Three tests have shipped vacuous
in this repo already, all three through review.

## 6. Follow-ups found while measuring -- file separately, do not fold in

### 6.1 A TRAILING `$MODEL` compartment with no `DADT` destroys every state name

Same shape as 2.4, gap moved to the end:

```
$MODEL
  COMP=(DEPOT, DEFDOSE)
  COMP=(CENTRAL, DEFOBS)
  COMP=(DUMMY)
```

nonmem2rx drops DUMMY entirely and loses the remaining names, so the emitted
file reads:

```
  ode(obs_cmt=rxddta2, states=[rxddta1, rxddta2])
  d/dt(rxddta1) = -KA * rxddta1
  d/dt(rxddta2) = KA * rxddta1 - (CL/V) * rxddta2
```

Internal nonmem2rx placeholders in user-facing output; three NONMEM
compartments become two ferx ones, so a `CMT=3` row has nowhere to go. It also
fires the #24 DEFOBS cross-check into a false accusation -- "the observation
expression reads compartment 'rxddta2' but $MODEL declares 'CENTRAL'" -- which
blames a source contradiction that does not exist. Different root cause
(upstream name resolution, not numbering), so it gets its own issue.

### 6.2 NONMEM doses by DEFDOSE when the data has no CMT column; ferx doses compartment 1

Found while evaluating this plan, read in current `ferx-core` `origin/main`:
`src/io/datareader.rs` resolves a dose row's compartment as
`cmt_col.and_then(...).unwrap_or(1)` -- a missing CMT cell, or no CMT column at
all, doses **compartment 1**. NONMEM's rule for the same dataset is the
`$MODEL` DEFDOSE compartment. ferx has no model-side dose-compartment concept
to translate DEFDOSE into (`dose_cmt` exists only in the `[simulation]` block,
which is design generation, not fitting), and the translator never reads
DEFDOSE at all.

So a control stream whose `$INPUT` lists no CMT and whose DEFDOSE is not
compartment 1 doses the wrong compartment after translation, silently. Not
made worse or better by this fix -- it is a missing-input default, not a
numbering permutation -- and the translator can DETECT it (`$INPUT` is in the
ctl, DEFDOSE is one `grepl` away from the DEFO code that exists since #24).
Per CLAUDE.md's rule for constructs ferx cannot express: ERROR warning plus
`unsupported` entry. Separate issue, separate PR.

### 6.3 Coupling with #16

[#16](https://github.com/FeRx-NLME/translator/issues/16) will emit `F1`,
`ALAG1`, `D1`, `R1`, whose trailing digit is a compartment number resolved the
same positional way. Landing #16 before this makes the numbering defect wider.
Order matters: **this first**.

## 7. Review pass -- what evaluating this plan changed

Adversarial pass over the plan itself before publication, per house rule.
Three claims did not survive:

1. Section 4.2 originally required threading the canonical list into the
   expression pass for `rest`. False -- verified order-independent; the
   section now records the verification instead of the requirement.
2. The plan originally left the #24 decline branches untouched and only
   scheduled a vacuity re-check on their tests. That preserved declines whose
   stated reason becomes false once the emitted numbering follows COMP order.
   Section 4.3 now resolves through `cmt_order` and re-points the tests.
3. The plan said "replace :531 and :690 with one resolution" without naming
   that `match()` loses `:690`'s case-insensitivity and that the permutation
   must be computed on raw names but applied to sanitised ones. Section 4.4
   names both traps.

One new adjacent defect surfaced (6.2: DEFDOSE vs ferx's missing-CMT default
of compartment 1) -- filed separately, not folded in.

## 7b. During implementation -- where the plan above was wrong

Recorded because the plan text stays in the repo and a wrong line in it costs
the next reader more than it saved.

**Section 5 claimed the Tier 4 test "cannot pass by accident". It could, and on
first writing it did.** The reasoning was that the broken case predicts
identically 0, so a fit on it is unmissable. Measured, that is only half true:
with IPRED == 0 there is no gradient at all, so the optimiser returns every
theta EXACTLY where it started -- and `cmt_order_gap.ctl`'s `$THETA` initials
are the simulation truth. Broken, the fit returned `TVCL=3, TVV=50, TVKA=1.2`
against a reference of `(3.0, 50.0, 1.2)`, deviation ~0, green. A sabotage run
found it: mutating the renumbering away produced ONE failure, the `states=[...]`
pin, and all three theta assertions passed.

The fix is the one the ODE omega guard already uses and that CLAUDE.md states
as a rule: start away from the truth. Each theta now starts at half its true
value, and the test additionally asserts each MOVED. Measured after: broken
freezes at `(1.5, 25, 0.6)` and the same mutation now produces seven failures.

This is the fourth vacuous test this repo has produced, and the first that a
sabotage run caught before merge rather than after. "A model's initials are
often the simulation truth" was written down for omegas; it is not an omega
property, it is a property of every fixture generated from its own model.

**The `length(comps) != length(state_raw)` guard in `.nm_cmt_order()` was
undiscriminated as first tested.** Every case the plan listed for it -- a
dropped compartment, a COMP matching no state -- is also caught by the
per-COMP check or the bijection check, so removing the guard changed no test
result. Rather than delete it as dead, the case only it catches was found:
a REPEATED COMP name gives `perm = c(1, 2, 1)`, every COMP finds exactly one
state so the per-COMP check passes, and `setequal()` ignores duplicates so the
bijection check passes too -- returning a three-element state list for a
two-state model. The unit test now pins that, and mutation confirms the guard
is the only thing catching it.

**Section 4.2's verification held.** `.assemble_endpoints()` needed nothing;
the expression pass was not touched.

**No snapshot moved**, as section 2.5 predicted. Suite went 1038 -> 1069 with
zero churn in `_snaps/`.

## 8. Risks

- **`.same_cmt_name()` becomes load-bearing for the dose.** It already gates
  `obs_cmt`; this makes a false negative reorder nothing (safe, warns) and a
  false positive build a wrong permutation (not safe). Its one transformation
  is strip `^c[.]`, `.ferx_ident()`, uppercase. Worth a unit test aimed at a
  false positive -- two COMPs that normalise to the same string -- which the
  bijection check should catch as non-unique.
- **Measured through `ferx_simulate()`, not `ferx_fit()`.** The dose application
  read in `predictions.rs` is the shared ODE path, so the fit is expected to
  follow, but it was not separately measured. The Tier 4 test closes this
  properly rather than by argument.
- **`nonmem2rx`'s placement of an invented state is behaviour, not contract.**
  DUMMY-first (2.4) and name-loss (6.1) are both observed, neither documented.
  A nonmem2rx upgrade could change either. The fix does not depend on the
  placement -- it reorders by name whatever it is handed -- but the fixtures do,
  so a nonmem2rx bump should re-run them rather than assume.
- **Scope temptation.** `[odes]` reordering, 6.1, 6.2, and #16 are all
  adjacent and all out. One PR per logical change.
- **4.3 widens the diff into #24's tiers.** The alternative -- numbering fix
  now, tier upgrade later -- was considered and rejected: it ships a release
  where the file's numbering and the tiers' warnings contradict each other,
  which is a worse intermediate state than a larger reviewed diff. Still one
  logical change: compartment-number reconciliation.
