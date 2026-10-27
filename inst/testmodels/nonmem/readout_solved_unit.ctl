$PROBLEM Solved system, readout converts units on F (readout_solved_unit)
; S1 = V is correct by itself (the pk macro already returns A/v); the factor under
; test is the unit conversion on F in $ERROR, which [scaling] obs_scale must carry.
$INPUT ID TIME DV AMT EVID MDV CMT
$DATA readout_arith_concordance.csv IGNORE=@
$SUBROUTINES ADVAN1 TRANS2
$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
  MW = 500
  S1 = V
$ERROR
  IPRED = F/MW*1000000
  Y = IPRED*(1+EPS(1))
$THETA
  (0,3)    ; CL
  (0,30)   ; V
$OMEGA 0.09
$SIGMA 0.01
$EST METHOD=1 INTER
