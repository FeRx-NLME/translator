$PROBLEM A $PK intermediate that nonmem2rx renames (rxm_alias_update)
; Regression fixture for #51. An intermediate updated in place after its default.
; CLSEX starts with C, so on a solved-system model nonmem2rx renames it to
; `rxm.CLSEX`. The name is the whole trigger: do not rename it.
$INPUT ID TIME DV AMT EVID MDV CMT SEX
$DATA rxm_alias_concordance.csv IGNORE=@
$SUBROUTINES ADVAN1 TRANS2
$PK
  CLSEX = 1
  CLSEX = CLSEX*THETA(3)**SEX
  CL = THETA(1)*CLSEX*EXP(ETA(1))
  V  = THETA(2)
  S1 = V
$ERROR
  IPRED = F
  Y = IPRED*(1+EPS(1))
$THETA
  (0,3)    ; CL
  (0,30)   ; V
  (0,1.5)  ; SEXFAC
$OMEGA 0.09
$SIGMA 0.01
$EST METHOD=1 INTER
