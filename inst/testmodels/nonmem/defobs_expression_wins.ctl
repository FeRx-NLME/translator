$PROBLEM DEFOBS says DEPOT but $ERROR reads A(2)
;
; Regression fixture: $MODEL and $ERROR disagree about the observed
; compartment, and $ERROR wins. DEFOBS is NONMEM's default for data records
; with no CMT; it says nothing about what the prediction reads. Here $ERROR
; names A(2) outright, so CENT is observed and S2 = V is its scaling.
;
; DO NOT "FIX" THE DISAGREEMENT. It is the subject under test. Preferring
; DEFOBS here emitted obs_cmt=DEPOT with no [scaling] block and no warning --
; a model that regressed depot amount against concentration data, unscaled.
;
$INPUT ID TIME AMT DV MDV
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=6
$MODEL
  COMP=(DEPOT, DEFDOSE, DEFOBS)
  COMP=(CENT)
$PK
  KA = THETA(1)*EXP(ETA(1))
  CL = THETA(2)
  V  = THETA(3)
  S2 = V
$DES
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1) - CL/V*A(2)
$ERROR
  Y = A(2)/S2*(1+EPS(1))
$THETA (0,1) (0,5) (0,20)
$OMEGA 0.09
$SIGMA 0.04
$ESTIMATION METHOD=1
