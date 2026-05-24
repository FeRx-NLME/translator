# ferxtranslate

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

```r
# Install ferxtranslate
install.packages("ferxtranslate")   # once on CRAN

# Install the parsing backends you need:
install.packages("nonmem2rx")       # for NONMEM
install.packages("rxode2")          # for nlmixr2
install.packages("monolix2rx")      # for Monolix (requires Monolix installation)
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
INFO  | THETA TVCL was FIXED in source -- treated as free in ferx
WARN  | Theta bounds not found for THETA2 -- set to (init, 0, Inf), please review
ERROR | MIXTURE model detected -- not supported in ferx, removed from output
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
| 1/2-cpt oral and IV bolus | Translated |
| Covariates (allometric power, linear) | Translated |
| Block omega | Translated |
| IOV (diagonal kappas) | Translated |
| ODE models (`$DES` / `d/dt()`) | Translated |
| Proportional / additive / combined error | Translated |
| FIXED thetas | Translated (emits `FIX` in ferx) |
| 3-compartment (oral / IV bolus) | Translated |
| Multiple DVIDs | ERROR -- not yet in ferx |
| MIXTURE models | ERROR -- not yet in ferx |
| IOV block omega | WARN -- diagonal only emitted |

## Learn more

```r
vignette("translating-nonmem")   # step-by-step NONMEM walkthrough
?to_ferx                         # main function reference
?ferx_ir                         # intermediate representation
```
