$PROBLEM Solved system, S1 = V/1000 (scale_quot_solved)
; Observation scale is an EXPRESSION on a linCmt system. The pk macro already
; divides by V, so only what S1 carries beyond V may reach [scaling]. Anchored
; on the NONMEM definition F = A(1)/S1.
$INPUT ID TIME DV AMT EVID MDV CMT
$DATA scale_expr_concordance.csv IGNORE=@
$SUBROUTINES ADVAN1 TRANS2
$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
  S1 = V/1000
$ERROR
  Y = F*(1+EPS(1))
$THETA
  (0,3)    ; CL
  (0,30)   ; V
$OMEGA 0.09
$SIGMA 0.01
$EST METHOD=1 INTER
