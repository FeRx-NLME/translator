$PROBLEM ODE readout divides by S1 (readout_ode_s1_on_f)
; Control for readout_ode_vol_s1: here $ERROR divides by S1 itself, so S1 DOES
; apply and the existing `obs_scale = V` is already right. Nothing may change.
$INPUT ID TIME DV AMT EVID MDV CMT
$DATA readout_arith_concordance.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=6
$MODEL COMP=(CENT)
$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
  S1 = V
$DES
  DADT(1) = -CL/V*A(1)
$ERROR
  CC = A(1)/S1
  IPRED = CC
  Y = IPRED*(1+EPS(1))
$THETA
  (0,3)    ; CL
  (0,30)   ; V
$OMEGA 0.09
$SIGMA 0.01
$EST METHOD=1 INTER
