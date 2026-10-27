$PROBLEM ODE readout divides the amount by a parameter (readout_ode_vol)
; The ONLY thing turning the compartment amount into a concentration is the
; division in $ERROR: there is no S1. The translation must carry that division
; into [scaling] y, or the file predicts the amount.
$INPUT ID TIME DV AMT EVID MDV CMT
$DATA readout_arith_concordance.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=6
$MODEL COMP=(CENT)
$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
$DES
  DADT(1) = -CL/V*A(1)
$ERROR
  IPRED = A(1)/V
  Y = IPRED*(1+EPS(1))
$THETA
  (0,3)    ; CL
  (0,30)   ; V
$OMEGA 0.09
$SIGMA 0.01
$EST METHOD=1 INTER
