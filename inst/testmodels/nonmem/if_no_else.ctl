$PROBLEM  One-compartment IV bolus with a covariate effect on CL written as two IFs with no default

$INPUT  ID TIME AMT DV SEX

$DATA   if_no_else.csv IGNORE=@

$SUBROUTINE ADVAN1 TRANS2

$PK
  IF (SEX.EQ.0) TVCL = THETA(1)
  IF (SEX.EQ.1) TVCL = THETA(1) * THETA(3)
  CL = TVCL * EXP(ETA(1))
  V  = THETA(2)
  S1 = V

$ERROR
  IPRED = F
  Y     = IPRED * (1 + EPS(1))

$THETA
  (0, 3.0)   ; TVCL_F
  (0, 30.0)  ; TVV
  (0, 1.5)   ; SEX_ON_CL

$OMEGA
  0.04   ; ETA_CL

$SIGMA
  0.02   ; PROP_ERR

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5
