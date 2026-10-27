$PROBLEM Solved system, readout adds a constant to F (readout_solved_offset)
; $ERROR adds a constant to F. [scaling] has no way to say that, so the
; translation must report it rather than predict F.
$INPUT ID TIME DV AMT EVID MDV CMT
$DATA readout_arith_concordance.csv IGNORE=@
$SUBROUTINES ADVAN1 TRANS2
$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
  S1 = V
$ERROR
  IPRED = F + 5
  Y = IPRED*(1+EPS(1))
$THETA
  (0,3)    ; CL
  (0,30)   ; V
$OMEGA 0.09
$SIGMA 0.01
$EST METHOD=1 INTER
