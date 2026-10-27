$PROBLEM A $THETA labelled with the name of the covariate it scales
; Regression fixture for #52. THETA(3) is the WT exponent and its comment label
; is `; WT`, the name of the $INPUT data item it scales -- an ordinary habit.
; nonmem2rx keeps the two apart (`t.WT` and `WT`), but stripping the `t.` prefix
; named the theta WT, and ferx resolves theta before covariate. The covariate
; then vanished from the model: `CL = TVCL * (WT/70)^WT * exp(ETA1)`, both WTs the
; theta, validating clean with nothing in $unsupported.
;
; The label is the whole difference. The test swaps it for `; WT_CL` to build the
; control, so do not add anything else that depends on it.
;
; The $THETA initials are deliberately AWAY from the simulation truth (CL = 5,
; V = 50, exponent 0.75): with the covariate shadowed, the fit drives the
; exponent towards 0, and an initial at the truth would make "recovered it" and
; "never moved" the same observation.
$INPUT ID TIME DV AMT EVID MDV CMT WT
$DATA theta_label_covariate_concordance.csv IGNORE=@
$SUBROUTINES ADVAN1 TRANS2
$PK
  CL = THETA(1)*(WT/70)**THETA(3)*EXP(ETA(1))
  V  = THETA(2)
  S1 = V
$ERROR
  IPRED = F
  Y = IPRED*(1+EPS(1))
$THETA
  (0,3)    ; CL
  (0,30)   ; V
  (0.3)    ; WT
$OMEGA 0.09
$SIGMA 0.01
$EST METHOD=1 INTER
