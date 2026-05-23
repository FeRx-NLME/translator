<!-- Title format: type(scope): short description  [closes #N] -->
<!--
  type  : feat | fix | refactor | docs | test | chore
  scope : nonmem | nlmixr2 | monolix | mrgsolve | emit | ir | api
  e.g.  : feat(nonmem): map ADVAN4 to two_cpt_oral pk macro [closes #12]
-->

## Why
<!-- What problem does this solve? Which source format or translation gap does it address? -->

## What changed
<!-- Functions added/changed, new IR fields, new warning messages, emit format changes. -->

## Alternatives considered
<!-- What else was evaluated and why this approach won. Omit if obvious. -->

## Formats affected
- [ ] NONMEM
- [ ] nlmixr2
- [ ] Monolix
- [ ] mrgsolve (v0.2)
- [ ] Shared (ferx_ir / emit_ferx / result class)

## ferx format dependency
<!-- Does this PR generate new .ferx syntax? If so, confirm the ferx parser already
     accepts it, or link the ferx-core PR that adds support. -->
- [ ] No new `.ferx` syntax generated
- [ ] New syntax confirmed accepted by current ferx parser (tested locally)
- [ ] Requires ferx-core change — linked PR: FeRx-NLME/ferx-core#___

## Breaking changes
- [ ] `ferx_ir` structure changed (fields added / renamed / removed)
- [ ] `ferx_translate_result` fields changed
- [ ] Warning message text changed (affects users parsing `result$warnings`)
- [ ] `.ferx` output format changed for existing models (snapshot diffs expected)
- [ ] None

## Tests
- [ ] Tier 1 unit tests added / updated
- [ ] Tier 2 integration tests added / updated
- [ ] Tier 3 reference snapshots updated (`testthat::snapshot_accept()` run after review)
- [ ] `R CMD check --as-cran` passes — zero ERRORs, zero WARNINGs
- [ ] Non-ASCII check run on `R/` — clean

## Documentation
- [ ] `roxygen2::roxygenize()` run and `man/` committed
- [ ] `@examples` updated for any changed function
- [ ] `plans/v0.1-implementation.md` updated if scope or design changed
- [ ] `NEWS.md` entry added

## Translation quality
<!-- For changes that affect emit output: paste a before/after of the .ferx text
     for a representative model to show the change is correct. -->

<details>
<summary>Before / after .ferx output (if applicable)</summary>

```
# before

# after
```

</details>

## Reviewer hints
<!-- Where to focus. What is subtle. What can be skimmed. -->

## Open questions
<!-- Things you are uncertain about and want input on. -->
