# Issues #26 and #27 -- two `$MODEL`/`$INPUT` constructs the translator mishandles

Plan for [#26](https://github.com/FeRx-NLME/translator/issues/26) and
[#27](https://github.com/FeRx-NLME/translator/issues/27). Written against `main`
`5668318`, measured with `ferx` 0.3.0 (`ferx-r@9c97c13`, the current engine pin).

**Two PRs, not one.** They share a subject area and nothing else: #26 is upstream
name loss we have to work around, #27 is a semantic default ferx cannot express.
Different root causes, different fixes, different fixtures.

**Order: #27 first.** It is self-contained and touches no code #26 touches. #26
modifies `.nm_cmt_order()`, which shipped four days ago with its own discriminating
tests, and that is the riskier edit of the two.

---

## Part 1 -- #27: a CMT-less dataset doses the wrong compartment

### What was measured

An ordinary control stream with **no `CMT` in `$INPUT`** and DEFDOSE on
compartment 2. `$MODEL` order is already COMP order, so this is **not** #25:

```
$INPUT ID TIME AMT DV MDV          ; <- no CMT
$MODEL
  COMP=(PERIPH)
  COMP=(CENTRAL, DEFDOSE, DEFOBS)
```

Emitted `ode(obs_cmt=CENTRAL, states=[PERIPH, CENTRAL])` -- correct in every
respect the translator currently checks. Simulated against the source's own
CMT-less dataset:

```
                                         t=1      t=4     t=12     t=24
emitted (ferx doses cmt 1 = PERIPH) 0.265796 0.641802 0.694758 0.511258
NONMEM meaning (dose into CENTRAL)  1.626028 1.020358 0.605598 0.424725
```

A factor of ~6 at t=1. Warnings emitted about it: **none**.

### Mechanism

`ferx-core/src/io/datareader.rs` resolves a dose row's compartment as
`cmt_col.and_then(...).unwrap_or(1)` -- absent column or absent cell means
compartment 1. There are three such sites (1618, 1748, 1868). NM-TRAN's rule for
the same dataset is the `$MODEL` DEFDOSE compartment.

ferx has no model-side dose-compartment binding on the fitting path: `dose_cmt`
exists only on `SimulationSpec`, parsed from the `[simulation]` block and read
only by the simulation path. So this is **not** a matter of emitting the right
keyword -- there is no keyword.

### Design

Detect and report. Per CLAUDE.md's rule for a construct we detect but ferx cannot
express: `ERROR`-level warning, `# WARNING:` comment in the output, and a
`result$unsupported` entry naming the compartment NONMEM would have dosed.

Fires only when BOTH hold:

1. `$INPUT` lists no `CMT` data item, **and**
2. the DEFDOSE compartment is not number 1.

Either alone is silence: with a `CMT` column the data decides and DEFDOSE is
irrelevant; with DEFDOSE on compartment 1 the two rules agree.

Infrastructure is nearly all present. `.nm_block_lines()` in `R/utils.R` already
owns block extraction -- note its comment about one owner deciding where a record
ends, which is the reason `$INPUT` should be read through it rather than by a
fresh `grep`. `.extract_nm_defobs()` already parses the full COMP list WITH
attributes and matches `^DEFO`; DEFDOSE is `^DEFD`, unambiguous because the two
attributes diverge at the fourth character. Extending the existing function to
return `defdose` alongside `defobs`/`comps` is cheaper than a second parser and
keeps one owner of the COMP list -- which is the property #25 depends on.

### Deliberately NOT done

**Do not reorder `states=[...]` to float DEFDOSE into position 1.** That would
fight #25, which pins the emitted order to `$MODEL` COMP order precisely so the
source's own `CMT` values keep meaning what the source meant. The two rules
collide only when there is no CMT column, and there the honest move is to say so.

If ferx-core later grows a model-side dose-compartment binding, this ERROR
becomes an emission. Worth a linked ferx-core issue if the team wants it -- the
peer working there has already raised the expressiveness gap with their user.

### Tests

- **Tier 1**: the `$INPUT`-has-CMT predicate, and DEFDOSE extraction -- including
  the abbreviation cases (`DEFD`, `DEFDOS`, `DEFDOSE`) that the `^DEFO` bug
  taught us to cover, and a `COMP=(X)` with no attributes.
- **Tier 2**: a bundled fixture with BOTH halves. Neither alone discriminates --
  with a CMT column the code path is dead; with DEFDOSE on compartment 1 the
  wrong code gives the right answer. Assert the ERROR, the `# WARNING:` comment
  and the `unsupported` entry.
- **Negative fixtures**: CMT present + DEFDOSE elsewhere -> silent; no CMT +
  DEFDOSE on 1 -> silent. Both must be shown to fail if the guard is widened.
- Every bundled model in the corpus currently agrees by accident (depot is both
  compartment 1 and DEFDOSE), so **nothing existing can show this** -- new
  fixtures are mandatory and each must be shown to fail first.

### During implementation -- where this section was incomplete

Landed on `fix/issue-27-defdose-no-cmt`, measured against `main` `f55d1fe` with
`ferx` installed locally. The design above survived intact; five things it did
not say came up.

**Silence on a missing DEFDOSE is CORRECT, not merely cautious.** The plan
framed the guard as "fires only when both hold" and left `defdose = NA` implicit.
It matters, because two bundled models reach it (`ode_warfarin.ctl`,
`s_scaling_not_last.ctl`): `$MODEL` present, no `DEFDOSE` attribute anywhere.
NONMEM's own default there is compartment 1, which is exactly what ferx does, so
the two agree and there is nothing to report. `.extract_nm_defobs()` returns
`NA_integer_` rather than filling in the 1, on the same principle its `defobs`
already follows -- filling it in would turn an agreement into a claim the control
stream never made.

**`CMT=DROP` needed a decision the plan did not anticipate.** NM-TRAN's `DROP` /
`SKIP` marks a column it will not read, so for the dose-compartment question a
DROPped `CMT` is identical to an absent one and `.nm_input_has_cmt()` answers
FALSE. That leaves one narrower divergence deliberately unreported and named in
the code: ferx binds columns by CSV header and knows nothing of `$INPUT`, so a
DROPped `CMT` column is still read by ferx even when `DEFDOSE` is compartment 1,
where this guard stays quiet. Chasing that would have widened #27 into a
different defect.

**`$INP`, not `$INPUT`.** NM-TRAN accepts any unambiguous record prefix, and
`$IN` collides with `$INFN`. `$INP` is the shortest that does not. The corpus
already contains a continuation-style `$INPUT` whose data items are on the NEXT
line (`pk_1cmt_oral.mod`), so reading only the `$INPUT` line would have answered
FALSE for a real model -- pinned by a unit test and by mutation m6.

**The `# WARNING:` comment needed a new mechanism.** `emit_ferx()` had a
per-statement `note` for `[odes]` and `init()` lines but nothing for
`[structural_model]`, so the header block was the only place a warning could
land. `ir$structural$note` now renders as a comment line directly above the
`ode(...)` line. The test for it is anchored by POSITION, not by slicing the
block out with a regex: `sub()` returns its whole input when it fails to match,
so a slice that silently degrades to the entire file still finds the note in the
header and passes.

**Scoped to ODE models.** `DEFDOSE` is a `$MODEL` attribute and a `$MODEL` block
means a general ADVAN, so the evidence only exists there. What ferx's pk macro
path does with a dose compartment was not measured, and firing on it would have
been a guess.

Reported the way CLAUDE.md prescribes, with the upstream gap now filed as
ferx-core#1009 (open, no comments): if that lands, this report becomes an
emission.

#### Mutation results

Ten mutations, each run against the full `rxui_to_ir` + `integration` suites.

| mutation | verdict |
|---|---|
| `isFALSE(has_cmt_col)` -> `!isTRUE(...)`, so NA fires | CAUGHT |
| drop the CMT half of the guard | CAUGHT |
| drop the `dd != 1L` half of the guard | CAUGHT |
| `^DEFD` -> `^DEFDOSE`, rejecting the legal abbreviations | CAUGHT |
| drop the `DROP`/`SKIP` check | CAUGHT |
| read only the first `$INPUT` line, not the record | CAUGHT |
| do not set `structural$note` | CAUGHT |
| `parts == "CMT"` -> `grepl("CMT", parts)` | CAUGHT |
| delete the `unsupported` append | CAUGHT |
| emit the note BELOW the `ode()` line | CAUGHT |

Two harness faults surfaced first and are worth recording, because both read as
a green result:

- `reporter = "summary"` prints failures as `1. Failure (...)`, so a harness
  matching lines that START with `Failure` found none and reported the first
  mutation VACUOUS when the suite had in fact caught it. Parse the
  `[ FAIL n | ... | PASS n ]` line instead -- every reporter emits it.
- The first attempt at the `unsupported` mutation renamed the entry's PREFIX and
  left `DEFDOSE` in the string, which is exactly what the assertions match on.
  It read VACUOUS while testing nothing. A mutation has to remove the thing the
  test asserts, not edit text next to it.

That is three harness faults across two PRs now (the `str.replace(..., 1)` that
hit the wrong function on #32 was the first). A VACUOUS verdict should be
treated as a claim about the harness until the mutation is confirmed applied at
the intended site.

---

## Part 2 -- #26: a trailing `$MODEL` gap destroys every state name

### The trigger, pinned by measurement

Not "a `$MODEL` compartment with no `DADT`" -- that is #25 and is already fixed.
The trigger is narrower and exact: **the highest-numbered `$MODEL` compartment has
no `DADT`.**

| fixture | `$MODEL` | `DADT` for | resulting states | names |
|---|---|---|---|---|
| `gap_mid.ctl` | DEPOT, CENTRAL, HOLE, TERM | 1, 2, 4 | `[HOLE, DEPOT, CENTRAL, TERM]` | **kept** -- HOLE materialised as `d/dt=0` and hoisted first (the #25 case) |
| `gap_tail.ctl` | DEPOT, CENTRAL, TAIL | 1, 2 | `[rxddta1, rxddta2]` | **lost**, TAIL dropped entirely |

So any gap below the top is materialised with names intact; only a gap AT the top
falls off the end and takes the name map with it.

### The names are recoverable, and the key survives a reordered `$DES`

`rxddta<i>` corresponds to `DADT(i)` -- the suffix IS the NONMEM compartment
number. Verified by dynamics on a 3-equation model:

```
d/dt(rxddta1) = -KA * rxddta1                                   <- DADT(1), depot
d/dt(rxddta2) = KA*rxddta1 - (CL/V)*rxddta2 - (Q/V)*rxddta2 ...  <- DADT(2), central
d/dt(rxddta3) = (Q/V)*rxddta2 - (Q/V)*rxddta3                    <- DADT(3), peripheral
```

and then re-verified against the one thing that could have made it coincidence:
writing `$DES` in the order 3, 2, 1 yields `[rxddta3, rxddta2, rxddta1]`, with
`d/dt(rxddta3)` still the peripheral equation. **The suffix follows the `DADT`
index, not statement order.** That is what makes the recovery rule safe.

So: `rxddta<i>` -> `comps[i]`.

### Design

Rename before anything else reads the state list, so every downstream consumer --
`.nm_cmt_order()`, the `obs_cmt` cascade, `[scaling]`, the per-CMT dispatch --
sees real names and needs no knowledge of this at all.

1. Detect: any raw state name matching `^rxddta[0-9]+$`, with a COMP list
   available.
2. Map each to `comps[i]` by its suffix. Decline the whole rename unless EVERY
   placeholder maps to an in-range, non-empty COMP name -- a partial rename is
   worse than none, since it produces a file that looks correct and is not.
3. Report at `INFO` that the names were recovered from `$MODEL`, naming them.
   This is a correctness-relevant change to identifiers the user will read.
4. Report at `ERROR` plus `unsupported` that compartment(s) `k..n` are declared in
   `$MODEL` but absent from the emitted model, so a `CMT=k` row cannot be
   expressed. Do not paper over it: ferx genuinely has fewer compartments than
   NONMEM here.

Point 4 is bounded good news. Since ferx-core#899, `check_dose_compartments`
**rejects** `cmt > n_states` on ODE models rather than silently dropping the row,
so a dataset naming the missing compartment fails loudly at fit time. The
`unsupported` entry exists so the user learns it before then, not instead of.

### The trap in this fix, named so it is not rediscovered

Recovering the names makes `.nm_cmt_order()` (from #25) suddenly applicable to
these models -- and it will **decline**, because `length(comps) != length(state_raw)`
(3 declared, 2 states). Relaxing that guard is the obvious next step and is where
the danger is: **the length check is the only thing catching a repeated COMP
name.** A duplicate gives `perm = c(1, 2, 1)`; every COMP finds exactly one state
so the per-COMP check passes, and `setequal()` ignores duplicates so the bijection
check passes too, returning a three-element state list for a two-state model.
That was found by mutation testing on #25 and has a unit test pinning it.

If `.nm_cmt_order()` is relaxed to allow a longer `comps`, it must gain a
replacement guard -- `!anyDuplicated(perm)` plus full coverage of
`seq_along(state_raw)` -- and the existing duplicate-COMP test MUST still fail
against the relaxed version. Re-run the #25 mutations; if any stops failing, the
guard has gone vacuous.

An alternative worth weighing at implementation time: since `rxddta<i>` gives the
NONMEM number directly, these models do not need name matching at all -- order by
the recovered suffix. That keeps `.nm_cmt_order()` untouched, at the cost of a
second ordering path. Prefer whichever leaves ONE owner of the ordering.

### Tests

- **Tier 1**: the placeholder-to-COMP mapping -- in range, out of range, empty
  COMP name, mixed placeholder and real names, no COMP list.
- **Tier 2**: a bundled fixture whose gap is **trailing**. A middle gap exercises
  #25 instead and passes whether this is fixed or not, which is exactly the
  discriminating-fixture rule. The fixture needs a DO-NOT-REORDER comment saying
  so, as `ode_theta_ref.ctl` and `s_scaling_not_last.ctl` do.
- **Tier 2**: assert the recovered names, the `INFO`, and the `ERROR` +
  `unsupported` for the dropped compartment.
- **Tier 4**: worth one. With names recovered the model should fit; the
  discriminating property is that it fits *to the right compartment*, so the
  fixture must dose `CMT=1` and observe a compartment whose identity would change
  if the recovery mapped wrongly. Start thetas away from truth -- a fixture
  generated from its own model makes "recovered the truth" and "never moved"
  indistinguishable, which has now bitten this repo four times.

---

## Risks

- **#26 edits code that shipped four days ago** with mutation-verified tests.
  Every #25 mutation must be re-run afterwards, not just the new tests.
- **`rxddta` is an undocumented nonmem2rx internal.** The prefix, the suffix
  convention, and the trailing-gap trigger are all observed, none contracted. A
  nonmem2rx bump can change any of them. The fixtures pin the behaviour; a bump
  should re-run them rather than assume. Recovery declining safely (point 2) is
  what keeps a changed convention from producing a wrong file instead of no
  rename.
- **#27's reachability is unmeasured in the wild.** Every bundled model agrees by
  accident. Oral PK usually has depot as both compartment 1 and DEFDOSE; IV-into-
  central models that number a peripheral first, and PKPD models, break it. If
  the ERROR proves noisy in practice, the response is better targeting, not
  downgrading it to a WARN -- the failure it describes is silent in the numbers.
