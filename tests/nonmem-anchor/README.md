# NONMEM anchors

Control streams whose NONMEM output pins a **NONMEM semantic** this package
depends on. They are not run by the test suite -- they need NONMEM, which CI does
not have. They are committed so the semantic is pinned rather than remembered,
and so anyone who doubts it can re-run it rather than re-derive it.

Run with NONMEM 7.6.0:

```bash
docker run --rm -v "$PWD/tests/nonmem-anchor":/w -w /w \
  --user "$(id -u):$(id -g)" nonmemdocker:V0.1 \
  bash -lc '/opt/NONMEM/nm760/run/nmfe76 anchorA.ctl anchorA.lst'
```

## `anchorA.ctl` / `anchorB.ctl` -- does `S<n>` apply to `A(n)`?

Identical but for one `$ERROR` line. All thetas `FIX`, `$OMEGA 0 FIX`,
`MAXEVAL=0`, `S2 = V` with `V = 10`, so `PRED` is the deterministic prediction:

```
A:  Y = F    *(1 + EPS(1))
B:  Y = A(2) *(1 + EPS(1))
```

Measured, NONMEM 7.6.0:

```
 TIME A_Y_eq_F B_Y_eq_A2 ratio_B_over_A
  0.5 3.631100  36.31100             10
  1.0 5.327700  53.27700             10
  2.0 5.906800  59.06800             10
  4.0 4.041100  40.41100             10
  8.0 1.291200  12.91200             10
 12.0 0.390250   3.90250             10
 24.0 0.010666   0.10666             10
```

`B/A = 10 = V` exactly, at every timepoint. So **`F` is the scaled prediction
(`A(n)/S<n>`) and a bare `A(n)` is the raw amount, which `S<n>` does not touch.**

That is what `.dv_is_scaled()` in `R/rxui_to_ir.R` decides, and issue #32 is the
defect it fixes: `obs_scale` used to be emitted for a bare `A(n)`, dividing by a
scale NONMEM never applied and putting every prediction out by a factor of
`S<n>`.

The corollary matters for fixtures: with `Y = A(n)` and no scaling, `CL` and `V`
are **not separately identifiable** -- `V` enters the ODE only as `CL/V`, so
amount data determines the ratio and nothing more. A concordance fit on such a
model sits at its starting values. Check identifiability before writing a Tier-4
test against an unscaled readout; `amount_readout.ctl` is checked by emitted text
and by simulation, deliberately not by fitting.
