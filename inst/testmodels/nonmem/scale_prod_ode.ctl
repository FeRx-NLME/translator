$PROBLEM ODE, S1 = V*MW (scale_prod_ode)
; Observation scale is an EXPRESSION, not a bare parameter. Anchored on the
; NONMEM definition F = A(1)/S1: the prediction is the amount over S1.
$INPUT ID TIME DV AMT EVID MDV CMT
$DATA scale_expr_concordance.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=6
$MODEL COMP=(CENT)
$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
  MW = 500
  S1 = V*MW
$DES
  DADT(1) = -CL/V*A(1)
$ERROR
  Y = A(1)/S1*(1+EPS(1))
$THETA
  (0,3)    ; CL
  (0,30)   ; V
$OMEGA 0.09
$SIGMA 0.01
$EST METHOD=1 INTER
