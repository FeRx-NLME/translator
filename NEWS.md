# ferxtranslate 0.1.0.9000

## Breaking changes

* `[scaling] obs_scale` and `ode(obs_cmt=...)` are no longer emitted alongside a
  `y` readout (issue #6, defects 15 and 12). ferx applies `obs_scale` on TOP of
  the readout rather than instead of it, which validates clean while moving
  every prediction -- by up to 11.37 units on the reporter's model -- and
  ignores `obs_cmt` entirely once `y` is present. `validate_ferx_ir()` now
  rejects an IR carrying both. `ferx_ir$scaling` gains `y` and `per_cmt`, and
  an `error_model` entry may carry `cmt` or `cond`. Only models whose `$ERROR`
  dispatches between endpoints are affected; every other model's output is
  unchanged.

* Warning text changed for untranslatable `$ERROR` expressions. The prediction
  named in "neither 1 nor the prediction ..." is now the expression with its
  epsilons zeroed, as the words say -- it previously printed the epsilon terms
  back, so the message named a "prediction" containing a sigma. Code matching
  on that message text needs updating; `result$warnings` structure is unchanged.

* The `sigma` declarations in `[parameters]` are now emitted in the order the
  `[error_model]` consumes them, rather than in `$SIGMA` order. ferx binds a
  single-endpoint error model's sigmas POSITIONALLY from the declaration order
  and discards the names written in `DV ~ combined(A, B)`, so a `$ERROR` that
  writes its additive term first -- `Y = F + EPS(1) + F*EPS(2)` -- previously
  produced a file in which the additive SD was applied proportionally and vice
  versa, with no diagnostic. Only models declaring more than one sigma are
  affected; the emitted output of every bundled model is unchanged.

* An error expression that cannot be translated no longer produces a guessed
  `proportional` model. `emit_ferx()` omits the `[error_model]` block entirely
  and writes a commented-out suggestion in its place, so the engine rejects the
  file (`E_MISSING_BLOCK`) and `to_ferx(strict = TRUE)` aborts. Previously this
  was a `WARN` on a file that fit and returned numbers. `ferx_ir` gains an
  `error_suggestion` field carrying those comment lines.


* Thetas whose name matches an individual parameter are now renamed
  (`theta CL` becomes `theta TVCL`, with `CL = TVCL * exp(ETA_CL)` unchanged).
  In ferx a theta silently shadows an identically named individual parameter in
  every block where thetas are in scope, so the individual definition was
  written and never read. Code that indexes `fit$theta` by name may need
  updating; the emitted output of `pk_1cmt_oral.mod` and
  `pk_1cmt_oral_ampsim.ctl` changes, and `pk_1cmt_oral.mod` now fits with
  inter-individual variability on clearance for the first time
  (`inst/testdata/ode_1cpt_oral_concordance.csv` was regenerated accordingly).

* `to_ferx()` and its wrappers now validate the emitted `.ferx` with
  `ferx::ferx_model_validate()` and, by default, abort when the engine rejects
  it. Previously a translation that could not be parsed was returned without
  comment. Pass `strict = FALSE` to warn instead, or `validate = FALSE` to skip
  the check.

## New features

* A NONMEM `$ERROR` block that dispatches between two endpoints is now
  translated instead of being reduced to one of them (issue #6, defect 5). The
  readout chain -- `CTOT = A(1)/VC`, `IPRED = CTOT`, `IF (FLAG.EQ.2) IPRED =
  RTOT` -- is rebuilt into a `[scaling]` readout, and the indicator structure of
  `Y = IPRED*(1 + W1*EPS(1) + W2*EPS(2))` into one error model per endpoint. A
  source switching on `CMT` emits `y[CMT=N]` and `CMT=N: DV ~ ...`; a source
  switching on any other column emits a Form C `y = if (...) ... else ...` and a
  covariate-selected `[error_model]`. The two forms are not interchangeable:
  ferx does not expose `CMT` to the covariate scope, so Form C cannot dispatch
  on it. The conditions must be equality tests on ONE column; anything else is
  reported rather than guessed at, and the pre-existing single-endpoint path
  runs unchanged, so no existing model's output moves.

* An endpoint a dispatch cannot express no longer sinks the whole translation.
  The `[scaling]` readout is derived per endpoint and independently of the error
  model, so it is emitted in full; the `[error_model]` becomes a commented block
  in the same dispatch shape with the one branch that needs a human marked
  `DV ~ ???` beside the source expression. ferx still rejects the file for the
  missing block, so nothing is silently accepted, but the mechanical half of the
  work is done. The report also names the endpoint and the epsilon actually at
  fault: previously such a model fell back to the single-endpoint path, which
  re-diagnosed the un-substituted `Y`, blamed whichever epsilon the indicators
  touched first, and suggested a single-endpoint `combined()` model of the wrong
  shape entirely.

* NONMEM `$ERROR` expressions are classified by STRUCTURE rather than by counting
  epsilons (issue #6, defects 10 and 11). The coefficient of each epsilon decides
  its role, so `Y = F + F*EPS(1)` is now proportional rather than additive, and a
  two-sigma model is emitted as `combined(proportional, additive)` whichever
  order the source wrote them in. Expressions with no single ferx equivalent --
  indicator-weighted multi-endpoint errors, exponential error, a scaled sigma --
  are reported with the offending coefficient named, instead of being reduced to
  a confident guess.


* `ferx_ir`'s `indiv_params` and `odes` are now ordered STATEMENT lists rather
  than flat assignment lists, and `emit_ferx()` renders them. A statement
  carries a `kind` -- `assign`, `ddt`, `init`, or `if` (with `cond`, `then` and
  an optional `else_` holding nested statement lists). This is what makes a
  source-model conditional representable at all; the translator does not yet
  produce `if` statements, so no translated model changes. An entry with no
  `kind` is read as the kind its block used before, so existing hand-built IRs
  keep working. Order is preserved exactly and nothing reorders it: ferx has no
  use-before-def check in `[odes]`, and an intermediate below the `d/dt` line
  that reads it stays valid while collapsing the prediction to a constant.

* The emitted `.ferx` records ODE state renames as `# renamed: state <from> ->
  <to>` provenance comments, and `ferx_ir` gains a `$state_renames` field. The
  `.ferx` is the artefact that gets shared, and a reader holding only that file
  could not otherwise map a sanitised state such as `c_RTOT` back to the
  `$MODEL` compartment or `A(n)` index it came from. These are provenance, not
  diagnostics: they do not count towards `# Warnings:`, and renamed thetas are
  deliberately excluded because `TVCL` is self-evidently derived from `CL`.

* `to_ferx(validate = TRUE, strict = TRUE)`. When the NONMEM `$DATA` file can be
  resolved relative to the control stream it is passed to the validator, which
  is what allows covariate references to be checked at all -- without a dataset
  every unknown identifier is read as a covariate and reported as valid.

* `ferx_translate_result` gains a `$validation` field holding `ok`,
  `diagnostics`, the `data_file` used, and the notes on how validation ran.

## Bug fixes

* **`$DES` and `$PK` conditionals are translated instead of discarded, and
  ODE intermediates are emitted rather than inlined.** An `IF`/`THEN`/`ELSE`
  anywhere in the model body used to be dropped with no diagnostic. In `$DES`
  that left the name it defines undefined in the output (the file did not
  parse); in `$PK` it silently deleted covariate effects, so any covariate model
  written the standard way was mistranslated with `warnings` and `unsupported`
  both empty. Conditionals now emit into `[odes]` or `[individual_parameters]`
  according to whether they depend on a compartment amount, with both arms and
  in source position.

  `$DES` intermediates are emitted as ODE-block intermediates in source order
  instead of being substituted into the `d/dt` lines. Substitution could not
  represent a variable defined inside a conditional, and its depth cap silently
  returned the un-inlined expression. The emitted `[odes]` block now reads like
  the source `$DES`. **The output of `pkpd_ir.mod` changes**: `C2` and `EFF` are
  emitted as intermediates rather than substituted, which is algebraically
  identical and computes `CENTRAL/V2` once instead of twice.

  A variable that reads a compartment amount now lands in `[odes]` rather than
  `[individual_parameters]` -- the latter is evaluated once per subject, so such
  a variable was being evaluated at the wrong frequency in a file that parsed
  and fitted. `$ERROR` indicator variables (`W1 = 0` / `W2 = 0`) are dropped
  rather than emitted as parameters the engine reports as never used.

  The TMDD model from issue #6 now translates to a `.ferx` the engine accepts
  with zero errors and zero warnings. Endpoint dispatch (issue #6 defect 5) is
  still not translated.

* **An ordinary parameter named like a ferx dose attribute is renamed.** ferx
  reads an `[individual_parameters]` name of the form `F{n}`, `D{n}`, `R{n}`,
  `ALAG{n}`, `LAGTIME{n}` (compartment `n` >= 1), or a bare `F`, `ALAG` or
  `LAGTIME`, as a dose attribute and applies it to the dose in addition to
  whatever the source model does with it. The emitted file validated clean,
  `$unsupported` was empty, and no engine diagnostic mentioned bioavailability,
  so every prediction was wrong by exactly the parameter's value with nothing to
  notice it -- measured at a factor of ten on the reported model. Such names are
  now renamed (`F1` becomes `F1_PAR`), every reference follows, and a `WARN`
  says what ferx would have done with the original. NONMEM sources are largely
  protected by NONMEM's own reservation of these names in `$PK`; nlmixr2 and
  Monolix reserve nothing, which is where this was reported from. No bundled
  model was affected, so no snapshot changed. Deliberate translation of a
  genuine `F1`/`ALAG1` remains unimplemented (#16).

* **A bundled model's fitted behaviour changes.** `pkpd_ir.mod` sets
  `A_0(4)=BL`, initialising the effect compartment to baseline. That statement
  was being discarded without a diagnostic, so the translated model started the
  compartment at 0 instead of `BL` and fitted different numbers from the source
  model. The emitted `.ferx` was valid either way and no test was looking for
  it, so nothing surfaced it. It now emits `init(EFFECT) = BL`. If you have
  fitted a translated indirect-response model with a non-zero baseline, the
  results change -- for the better, but they change.

* The observed compartment is decided by the DV expression first, `$MODEL`'s
  `DEFOBS` second, and a positional guess only as a last resort. `DEFOBS` is
  NONMEM's default for data records with no `CMT`; it says nothing about what
  `$ERROR` reads, so a model whose `$ERROR` names `A(2)` outright is observed on
  compartment 2 even when `DEFOBS` names another. When the DV expression carries
  no compartment of its own (`IPRE = F`), `DEFOBS` is the authority, because
  that is what NONMEM's `F` means. A disagreement between the two is reported
  rather than resolved silently.

* `$MODEL`'s `DEFOBS` attribute is now read from the control stream. `ui$central` is `NULL` for both nonmem2rx and
  rxode2, so `obs_cmt` was always the LAST declared compartment -- a positional
  guess that every bundled model happened to satisfy. The same index selects the
  `$PK` scaling variable, so both failed together: a model declaring
  `COMP=(CENT, DEFDOSE, DEFOBS)` before `COMP=(PERIPH)` with `S1 = V` attributed
  observations to the peripheral compartment AND emitted no `[scaling]` block at
  all, discarding a correctly parsed scaling with no diagnostic. The `$MODEL`
  ordinal is cross-checked against the differential-equation order and refused
  with a warning if the two disagree, rather than trusted. New bundled model
  `defobs_not_last.ctl` is the regression guard.

* ODE state initial conditions are translated. `A_0(n) = <expr>` in `$PK`
  becomes `init(<state>) = <expr>` inside `[odes]` (see the behaviour-change
  note above for the bundled model this corrects). ferx resolves an init
  expression against individual parameters, other states and literals but not
  thetas, so one referencing anything else is dropped with an explanation naming
  what was out of scope, rather than emitting a file the engine rejects.

* Statements whose assignment target is a call rather than a name -- `F1`,
  `ALAG1`, infusion `D1`/`R1` -- are no longer dropped in silence. Each is
  reported by name in `$warnings`. They are NOT added to `$unsupported`: ferx
  supports all of them (`f=`/`lagtime=` on a pk macro, the reserved `F` and
  `LAGTIME` names for ODE models), so not emitting them is a `ferxtranslate`
  limitation rather than a ferx feature gap, and `$unsupported` is the
  ferx-core prioritisation signal. Assignments that set the value ferx already
  uses (`F1 = 1`, `ALAG1 = 0`) are noted at `INFO`. An alias the source binds
  more than once -- the `F1 = 1` / `IF (FORM.EQ.2) F1 = THETA(4)` idiom -- is
  never treated as a constant.

* A name reachable only through such a dropped statement is no longer reported
  as an illegal covariate. `f(depot) <- BIO.AV` raised
  `ERROR | covariate reference(s) BIO.AV ...` telling the user to rename a data
  column, a remedy that fixed nothing because `BIO.AV` was never emitted, and
  added a phantom entry to `$unsupported`.

* Two source names that normalise onto the same spelling no longer merge into
  one random effect. `$OMEGA 0.09 ; CL.IIV` beside `$OMEGA 0.04 ; CL_IIV` both
  emitted `omega CL_IIV`, and every reference resolved to the first -- one
  inter-individual variability term silently dropped, the other double-counted,
  with the engine reporting no problem. Thetas and states already had a
  uniqueness check; the eta/omega/kappa/sigma channel did not.

* `validate_ferx_ir()` now rejects an IR that declares a name which is not a
  legal ferx identifier, an `obs_cmt` that is not one of the declared states, or
  an unnamed `state_renames` vector. Legality was previously enforced only where
  each name was minted, so a name reaching the file through an unchecked channel
  was emitted verbatim into a file the engine could not parse. Expression text
  is deliberately exempt: covariate references must keep the data column's exact
  spelling.

* A shadowed theta can no longer be "renamed" to its own name.
  `.free_theta_name()` could return the name it was asked to replace once both
  the `TV` and `THETA_` prefixes were taken, reporting a rename that had not
  happened while the shadowing survived.
* A theta referenced from an ODE is now carried in by an individual parameter.
  ferx resolves a `d/dt` right-hand side against states, individual parameters and
  covariates only -- a theta is not in scope there -- so a `$DES` line naming a
  theta produced an `[odes]` block ferx could not resolve (issue #6, defect 2).
  Measured against ferx 0.2.0 and 0.3.0, a theta IS readable from
  `[individual_parameters]`, `[scaling] y` and `obs_scale`, and is NOT readable
  from a `d/dt` right-hand side, an ODE-block intermediate, an `init()` expression
  or a pk macro argument; only the second group gets a carrier, so a theta read
  from `[individual_parameters]` is left exactly as it was.

  Where the source already supplies the assignment nothing is added: `KTP =
  THETA(3)` arrives as the alias `KTP <- KTP`, which the alias filter drops as a
  self-assignment but which theta de-shadowing turns into the surviving `KTP <-
  TVKTP` -- the same carrier, spelled by the model. `inst/testmodels/nonmem/ode_theta_ref.ctl`
  covers the shape that has nothing to convert.

  A carrier is named after the source (`KTP = TVKTP`) where that name is free,
  which is the form ferx's own examples use and keeps the emitted `[odes]` diffable
  against the source `$DES`. Where it is taken the name is derived from the theta's
  *emitted* name instead: `TVCL_ODE = TVCL`. Not `CL_ODE`, which reads as "the CL
  used in the ODE" -- that is the individual value, and distinguishing it from the
  theta is the entire point. Not `CL_1` either: `_1` already means "state
  disambiguated" (`central_1`) and is `.free_theta_name()`'s last resort, so a
  third meaning on that spelling leaves a reader unable to tell a renamed
  compartment from a renamed theta from a carrier, and `CL_1`/`V_1` are plausible
  model variables in their own right. A numbered form is still possible if both
  spellings are taken, and now warns: the number is positional, so it moves when
  another carrier is added before it.

  A carrier never takes a name ferx reserves for the solver. The de-shadow pass
  renames a theta off `TIME`/`T`/`TAFD`/`TAD`/`MACHEPS`, and the carrier used to
  take the theta's source name back, putting the collision one block lower with
  the ODE term reading the integrator's clock.

  A duplicate `$THETA` label produces one carrier, not two. The discovery
  predicate is per-theta and two thetas can share a source name, so collapsing
  its result to names let both through and defined a second, dead parameter
  beside the one the reference resolves to.

  Carriers are appended before `[scaling]` is resolved, and the order is
  load-bearing. `[scaling]` resolves `S2 = VC` by looking the name up among the
  theta names and then among the individual parameters, and a carrier moves the
  name that answers to `VC` from the first list to the second. Resolved first,
  the lookup found the theta renamed to `TVVC` and no parameter named `VC`, so
  `[scaling]` was dropped with no diagnostic -- an emitted model that predicts
  amounts against concentration data and that the engine validates clean.

* A carrier reuses an existing individual parameter only when that parameter is a
  *pure alias* of the theta. Matching on the name alone gets the arithmetic wrong:
  `frac <- central/cl` in a model that later writes `cl <- cl*exp(eta.cl)` reads
  the theta, so pointing that reference at the individual parameter `CL` silently
  substitutes the IIV-applied value. Both forms parse and both fit. Such a
  reference now gets its own carrier (`TVCL_ODE = TVCL`).

* The `[odes]` scope is applied to the emitted right-hand side rather than during
  the walk, because that is where every path converges. An ODE-block intermediate
  that touches a state is inlined instead of emitted, and its text was normalised
  for a different context: `ki <- KTP*CENT` reached `[odes]` as `TVKTP * CENT`, a
  bare theta, with the carrier defined and unreferenced beside it. Resolving the
  `d/dt` line alone also made the result depend on statement order -- the alias
  written above the `d/dt` bound in time, the identical line below it did not.

* A theta is recognised as ODE-referenced by its *current emitted* name, or by its
  source name when nothing else declares that symbol -- not by source name alone.
  In any de-shadowed model the source name is the individual parameter's name, so
  `d/dt(ABS) = -KA * ABS` beside `theta TVKA` and `KA = TVKA * exp(ETA_KA)` reads
  the parameter, and matching on `KA` invented a second carrier for a reference
  that was already correct. The test is recomputed each round rather than
  accumulated, for the same reason: a theta matches on its source name only until
  de-shadowing renames it.

* nonmem2rx's `rxmissingvars` placeholders are no longer emitted. When a NONMEM
  model names a theta by its `$THETA` label (`FLUX = KTP*A(1)` for `(0,0.2) ;
  KTP`), nonmem2rx does not bind the symbol: it leaves a free `KTP` and records
  the theta in `rxmissingvars1 <- t.KTP`. That produced the meaningless individual
  parameter `RXMISSINGVARS1 = TVKTP` while the reference itself stayed dangling.
  The placeholder is dropped and the reference bound by the carrier.

  Unlike the covariate case noted above, `[odes]` does not need a dataset for the
  engine to object: it reports `[odes]: RHS references undefined name(s)` either
  way, because thetas and etas are out of scope there and a bare name can only be
  a state, a parameter or a covariate.

* Every name the emitted `[odes]` block references is now checked against what the
  model declares, and an unresolved one is an `ERROR` with an `$unsupported` entry.
  ferx does report this itself, but only where the engine runs -- and the
  identifier-legality check beside this one tests the *grammar*, not whether
  anything declares the name, so `KTP` passed it. Two defects reached the corpus
  through that gap (issue #6 defects 2 and 4).

  The check is possible because the set of names a ferx ODE right-hand side may
  reference is *closed*: declared states, individual parameters, ODE-block
  intermediates, and `TIME`/`TAFD`/`TAD`/`MACHEPS`. Not thetas, etas or sigmas --
  and not covariates either. ferx's own message gives the remedy, which this one
  repeats: pre-compute the covariate-dependent term in `[individual_parameters]`
  and reference it from the ODE by that name.

  No covariate list is consulted, and that is deliberate rather than an omission.
  `.covariate_names()` *defines* a covariate as a symbol nothing else binds, and
  rxode2's `ui$allCovs` does much the same, so both classify a name the translator
  failed to bind as a legitimate covariate -- measured, both call the unbound `CF`
  in the QSS TMDD model and the unbound `KTP` in `ode_theta_ref.ctl` covariates.
  An earlier version of this check consulted them, passed the entire test suite,
  and could not fire on either defect it was written for.

* A state and a model parameter whose names differ only in case are now reported
  as the source-name collision they are. ferx compares them case-insensitively:
  measured against ferx 0.3.0, `d/dt(central) = -central * (CENTRAL/V)` beside
  `states=[central]` and `theta CENTRAL` validates clean and reports
  `W_UNUSED_PARAM` for the theta -- the engine read the ODE's `CENTRAL` as the
  compartment, the theta went dead, and the term became the amount squared over
  `V`. The collision only becomes visible once something drags the reference into
  `[odes]`: `kk <- CENTRAL/v` is an honest read of the theta in
  `[individual_parameters]` scope, but its right-hand side "references a state"
  under the same case-folded comparison, so the inliner moves the text into the
  one block where that spelling means something else.

  No function-name whitelist is consulted either, for a simpler reason: a call
  head is never collected. The symbol walk recurses over the *arguments* of a
  call, so `exp(-K*CENT)` yields `K` and `CENT` and nothing else, and a list of
  function names could only ever whitelist ordinary identifiers that happen to
  share a spelling. One did: with `EXP`/`LOG`/`ABS`/`MIN`/`MAX`/`SIGN` and the
  rest declared, an undeclared covariate named `MAX` passed this check while `WT`
  was reported, and ferx rejected the file the check had just called clean.

  Two reports, because the remedies differ: a theta, eta or sigma needs a carrier;
  anything else resolves to nothing at all. The first makes an eta referenced from
  an ODE reportable for the first time -- ferx rejects it as `E_PARSE` and it was
  previously emitted without comment.

  Scoped to `[odes]` on purpose. In `[individual_parameters]` thetas *are* in
  scope, so the set of legitimate names is far larger and a leftover carries much
  less information.

  Measured against ferx 0.3.0 (ferx-r tag `v0.3.0`), which is the version shipped
  to users. 0.2.0 was not re-measured; no bundled model references a covariate from
  `[odes]`, so nothing that translated before is affected either way. The `engine`
  CI job now pins that same build, so these measurements and the job agree.

* A theta renamed only because it is ODE-referenced is no longer reported as
  shadowing an individual parameter. At the point of the rename nothing of that
  name is one -- the trigger predicts a carrier that is created afterwards -- so
  the message described a collision the model never had.

* A theta whose name is also an ODE state name is not treated as an ODE reference:
  in `[odes]` the state wins. Without the exclusion the emitted ODE still came out
  right, but the theta was renamed for a reference that was never to it.


* Emitted identifiers are now sanitised to the ferx grammar
  (`[A-Za-z_][A-Za-z0-9_]*`). ODE state names never went through normalisation,
  so a compartment that `nonmem2rx` had renamed to `c.RTOT` -- which it does
  whenever a `$MODEL` compartment collides with a variable name -- made the whole
  file unparseable (issue #6, defect 1). The rename is applied to the
  declaration and to every reference at once (`obs_cmt=`, `states=[...]`, the
  `d/dt` target, and the state inlined into other ODE right-hand sides), because
  renaming the declaration alone leaves the references pointing at nothing.
  Only names that must change do change: an already-legal state keeps its name
  and its case.

* Covariate references are deliberately left alone, and an illegal one is now
  reported as an `ERROR` rather than renamed. ferx matches covariates to data
  columns by exact name, case included, so sanitising one would move the failure
  from `E_PARSE` to `E_MISSING_COVARIATE` at fit time instead of fixing it. The
  dataset column has to be renamed. A regression test pins the existing
  case-preserving behaviour, which nothing previously covered.

* Every ODE state -- renamed or not -- is pinned on the maps the model block is
  parsed with, so neither an `iniDf` key that happens to share a state's name nor
  an assignment walked earlier can rebind one. The general name map holds every
  parameter alias and grows as the walk proceeds, which made the emitted state
  name something nobody decided, on the declaration *and* on every reference:
  `d/dt(CENT) = -K * VC` when a theta keyed `CENT` was labelled `VC` -- the
  compartment amount replaced by a fixed theta, with no diagnostic and a file the
  engine accepts -- and an ODE whose right-hand side changed depending on whether
  `central <- 0` stood above or below the `d/dt` line. `obs_cmt` is resolved the
  same way; it previously read a different map from the `d/dt` target, so the two
  could disagree and emit an `obs_cmt` naming no declared state.

* A state whose source name is also a model parameter's is now resolved by
  scope -- the state inside `[odes]`, where thetas and etas are out of scope in
  ferx too, and the parameter everywhere else -- and reported as an `ERROR`,
  because the source cannot say which reading it meant. It was previously
  refused outright, which left the state sharing a name with the parameter: the
  assignment referencing it was absorbed as an ODE intermediate, dropped from
  `[individual_parameters]` (leaving the block empty), and self-inlined to the
  depth cap, emitting the same `exp(ETA_X)` factor fifteen times.

* A symbol referenced *before* its defining assignment is now emitted under the
  same name as one referenced after it. The alias was installed only once the
  walk had passed the assignment, so `d/dt(central) = -cl*central*f.rac` written
  above `f.rac <- 0.5` emitted the illegal `f.rac` verbatim beside a declared
  `F_RAC` -- two spellings of one variable, an unparseable file, and no warning,
  since the legality check treats every assignment target as legal by
  construction regardless of position. Dotted local names are idiomatic nlmixr2.

* A theta named after a ferx solver builtin (`TIME`, `T`, `TAFD`, `TAD`,
  `MACHEPS`) is renamed, and an individual parameter with such a name is
  reported. ferx-core checks its reserved list against states, individual
  parameters and ODE intermediates alike; only states were covered. The theta
  half was silent rather than loud: ferx resolves the bare name to the value the
  solver injects, so `KA = TIME * exp(ETA1)` made `KA` read the integrator clock
  and left the estimated theta unreferenced, reported only as a
  `W_UNUSED_PARAM`.

* A NONMEM sigma whose source name needs sanitising still reaches
  `[error_model]`. `nonmem2rx` keeps sigma out of `iniDf`, so nothing bound its
  source spelling to its emitted name and the two agreed only while both were
  plain uppercase; once the declaration was sanitised, the error assignment
  matched no sigma and the block was omitted entirely -- a model with no residual
  error, and an `ERROR` misreporting the sigma as an illegal covariate.

* The state sanitiser reserves every assignment target, covariate and reserved
  ODE name, not just the random-effect names. A sanitised state landing on an
  assignment target was absorbed into the auxiliary-variable set, dropped from
  `[individual_parameters]`, and its references re-resolved to the state --
  turning `d/dt(A.B) = -A_B * A.B` into `d/dt(A_B) = -A_B * A_B`, the amount
  squared, with the rate constant and its IIV gone. ferx validates that output
  clean, so nothing downstream would have caught it.

* The identifier-legality check runs on every symbol emitted verbatim, not on
  the covariate set. Covariates are classified by *normalised* name, so an
  illegal raw symbol whose normalised form matched a known name was filtered out
  before the check saw it, while still being emitted verbatim.

* Smaller fixes from the same review: `ui$central` is length-guarded (a
  zero-length or multi-element value aborted the translation); `.ferx_ident()`
  maps `NA` to a legal identifier rather than passing it through; `.ddt_state()`
  takes the first element, so a non-symbol `d/dt()` argument cannot invent a
  second state; and the two hand-inlined copies of the old normalisation rule
  (in `[scaling]` matching and in `.extract_sigmas()`) now call `.norm()`.

* `inst/testmodels/nonmem/dotted_state.ctl` joins the corpus, and the sweep in
  `test-integration.R` now tokenises the *emitted text* rather than listing IR
  declaration fields. An illegal name parses where it is declared and fails only
  where it is referenced, so a declaration-only check inspected the half that
  works -- it passed cleanly on a file the engine rejects.

* The concordance data generator (`data-raw/generate_concordance_data.R`)
  produced empty datasets against `ferx` 0.2.0: its templates wrote `DV = "."`,
  which the engine reads as "no observation" and skips. It now uses a numeric
  placeholder, restores `.` on dose rows when writing, and asserts that one
  value came back per observation row. Reported upstream as
  FeRx-NLME/ferx-r#283 and FeRx-NLME/ferx-core#957.

* The translation gap report now validates every bundled model and fails on
  engine errors, instead of only listing features the translator already knew it
  could not express.

## Internal

* The CI engine pin moves to `ferx-r@7889719` (0.3.0, tag `v0.3.0`) -- the version
  users install. It had stayed on `ferx-r@731adc9` (0.2.0) after 0.3.0 shipped,
  which meant the `engine` job proved only that a build nobody runs accepted the
  emitted `.ferx`. That is the one job standing between a silent engine regression
  and a merge, so pointing it at a stale engine defeats its purpose.

  No re-baselining was needed: the full suite passes against 0.3.0 with the
  concordance references in `test-concordance.R` and the bundled datasets in
  `inst/testdata/` unchanged. `CLAUDE.md` says both are tied to the pin and must be
  re-baselined alongside it; that stays true as guidance -- this bump happening to
  need none is not evidence the next one won't.

  The pin is a SHA rather than the tag on purpose: a tag can be moved.

* `ferx` and `amp.sim` are declared in `Suggests`; both were already used by the
  test suite without being declared.

* The CI engine pin moves from `ferx-r@54f25d4` (0.1.5) to `ferx-r@731adc9`
  (0.2.0), the build the concordance suite is now verified against. No reference
  *value* changed: the theta assertions were renamed to follow the `TV` prefix,
  two omega assertions were added, and `ode_1cpt_oral_concordance.csv` was
  regenerated because its source model changed. The other three datasets
  regenerate byte-identical against this pin. The 2-cpt omega references still
  carry their original 0.1.x provenance and were re-verified, not re-fitted.

* The `engine` CI job now runs the full test suite rather than
  `filter = "concordance"`. Four of the validation tests require `ferx`, so a
  filtered run left the output-validation feature executed by neither job.
