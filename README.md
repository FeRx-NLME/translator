# ferxtranslate

<!-- badges: start -->
[![check](https://github.com/FeRx-NLME/translator/actions/workflows/check.yml/badge.svg)](https://github.com/FeRx-NLME/translator/actions/workflows/check.yml)
<!-- badges: end -->

Translate pharmacometric models written in NONMEM, nlmixr2, or Monolix into
[ferx](https://github.com/FeRx-NLME/ferx-r) `.ferx` format, so you can move
to ferx without rewriting your models by hand.

## How it works

```
NONMEM .ctl  ──nonmem2rx──┐
nlmixr2 fn   ──rxode2──── rxUI ──rxui_to_ir──► ferx_ir ──emit_ferx──► .ferx
Monolix .mlx ──monolix2rx─┘
```

A single adapter (`rxui_to_ir`) converts the rxode2 UI object to a ferx
intermediate representation, then `emit_ferx` writes the `.ferx` file.
All three source formats share the same translation path.

## Installation

Install with [pak](https://pak.r-lib.org). One call resolves CRAN and GitHub
packages and pulls the `nonmem2rx` and `rxode2` parsing backends automatically.

```r
# install.packages("pak")   # if you don't have it yet

# Development version from GitHub:
pak::pak("FeRx-NLME/translator")

# ...or once published on CRAN:
pak::pak("ferxtranslate")

# Optional add-ons:
pak::pak("monolix2rx")        # for Monolix .mlxtran sources
pak::pak("FeRx-NLME/ferx-r")  # the ferx engine, to fit the translated .ferx (needs Rust/cargo)
```

## Quick start

```r
library(ferxtranslate)

# -- From NONMEM --------------------------------------------------------------
result <- nm_to_ferx("run001.ctl", output = "run001.ferx")
result$warnings         # INFO/WARN/ERROR messages from the translation
result$unsupported      # features that could not be translated

# -- From nlmixr2 -------------------------------------------------------------
my_model <- function() {
  ini({ tvcl <- 0.134; eta.cl ~ 0.07; err.prop ~ 0.01 })
  model({ cl <- tvcl * exp(eta.cl); linCmt() ~ prop(err.prop) })
}
result <- nlmixr2_to_ferx(my_model, output = "my_model.ferx")

# -- From Monolix (requires monolix2rx) ---------------------------------------
result <- mlx_to_ferx("project.mlxtran", output = "project.ferx")

# -- Inspect without writing to disk ------------------------------------------
cat(to_ferx("run001.ctl", format = "nonmem")$ferx_text)
```

## Translation warnings

Every result carries a `$warnings` vector with prefixed messages:

```
INFO  | S2 = V detected -- emitting [scaling] obs_scale = V
WARN  | complex $ERROR -- classified as proportional, verify
WARN  | ETA 'KAPPA_CL' looks like inter-occasion variability but was emitted as IIV
ERROR | No structural model detected -- [structural_model] section omitted
```

`ERROR`-level items also appear in `$unsupported` so you can programmatically
check whether the output is complete:

```r
if (length(result$unsupported) > 0) {
  warning("Manual editing required before running ferx_fit()")
}
```

## What ferxtranslate can and cannot translate

See `vignette("translating-nonmem")` for the full catalogue. Short version:

| Source feature | Status |
|---|---|
| 1, 2, 3-compartment, oral and IV | Translated |
| Covariates (allometric power, linear) | Translated |
| Block omega | Translated |
| IOV (diagonal kappas) | Translated [^iov] |
| ODE models (`$DES` / `d/dt()`) | Translated |
| Proportional / additive / combined error | Translated |
| FIXED thetas | Translated (emits `FIX` in ferx) |
| Multiple DVIDs / joint PK-PD | ferx supports it; translator does not emit it yet [^dvid] |
| MIXTURE models | Not yet in ferx |
| IOV block omega | WARN -- diagonal kappas only emitted |

[^iov]: ferx and the emitter fully support IOV (`kappa` + `iov_column`), and it
    translates cleanly when the source exposes occasion-level random effects
    (e.g. nlmixr2 `iov`). For NONMEM, `nonmem2rx` commonly reads an ETA-coded
    IOV term (`KAPPA = ETA(n)`) as ordinary IIV, so it arrives as an extra
    `omega` and the occasion structure is lost. This is no longer silent: the
    translator emits a `WARN` when a `KAPPA*`/`IOV*`-named eta lands in the
    `[omega]` block, telling you to declare it as a `kappa` and set
    `iov_column`.

[^dvid]: ferx supports multiple observation types via a per-CMT error model
    (`CMT=2: DV ~ proportional(...)`, `CMT=3: DV ~ additive(...)`, with per-CMT
    `y[CMT=N] = ...` readouts). The translator currently maps every observation
    to a single `DV` and does not emit the per-CMT dispatch.

## Development

The test suite runs in four tiers -- unit, integration, reference snapshot, and
numerical concordance (which fits translated models against pre-simulated data
and checks that the recovered parameters match). CI runs two jobs on every pull
request:

- **fast** -- `R CMD check` plus the unit, integration, and snapshot tiers.
- **engine** -- builds a pinned [ferx](https://github.com/FeRx-NLME/ferx-r) and
  runs the concordance tier, so a green check means the engine actually accepted
  and fit the emitted `.ferx`.

```r
devtools::test()   # concordance runs locally when ferx is installed
```

## Learn more

```r
vignette("translating-nonmem")   # step-by-step NONMEM walkthrough
?to_ferx                         # main function reference
?new_ferx_ir                     # intermediate representation
```
