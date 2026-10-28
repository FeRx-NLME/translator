$PROBLEM  One-compartment IV bolus with IOV on CL, written with BLOCK(1) SAME

$INPUT  ID TIME AMT DV OCC

$DATA   iov_same.csv IGNORE=@

$SUBROUTINE ADVAN1 TRANS2

$PK
  CL = THETA(1) * EXP(ETA(1))
  IF (OCC.EQ.1) CL = THETA(1) * EXP(ETA(1) + ETA(2))
  IF (OCC.EQ.2) CL = THETA(1) * EXP(ETA(1) + ETA(3))
  IF (OCC.EQ.3) CL = THETA(1) * EXP(ETA(1) + ETA(4))
  V  = THETA(2)
  S1 = V

$ERROR
  IPRED = F
  Y     = IPRED * (1 + EPS(1))

$THETA
  (0, 3.0)   ; TVCL
  (0, 30.0)  ; TVV

$OMEGA
  0.04   ; ETA_CL (IIV)

$OMEGA BLOCK(1)
  0.04   ; IOV_CL (one variance, shared by the three occasions)

$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME

$SIGMA
  0.02   ; PROP_ERR

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5
