$PROBLEM ODE readout divides the amount by a parameter (readout_ode_vol_s1)
; S1 = V is declared, but $ERROR reads the AMOUNT (NONMEM applies S1 only to F),
; so S1 does nothing here and the division in $ERROR is what makes a
; concentration. [scaling] must say `y = CENT/V` and not also divide by V.
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
  IPRED = A(1)/V
  Y = IPRED*(1+EPS(1))
$THETA
  (0,3)    ; CL
  (0,30)   ; V
$OMEGA 0.09
$SIGMA 0.01
$EST METHOD=1 INTER
