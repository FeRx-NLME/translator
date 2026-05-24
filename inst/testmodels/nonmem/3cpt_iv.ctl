$PROBLEM  Three-compartment IV bolus pharmacokinetics

$INPUT  ID TIME DV AMT EVID MDV

$DATA   3cpt_iv.csv IGNORE=@

$SUBROUTINE ADVAN11 TRANS4

$PK
  TVCL = THETA(1)
  TVV1 = THETA(2)
  TVQ2 = THETA(3)
  TVV2 = THETA(4)
  TVQ3 = THETA(5)
  TVV3 = THETA(6)

  CL = TVCL * EXP(ETA(1))
  V1 = TVV1 * EXP(ETA(2))
  Q2 = TVQ2 * EXP(ETA(3))
  V2 = TVV2 * EXP(ETA(4))
  Q3 = TVQ3
  V3 = TVV3
  S1 = V1

$ERROR
  IPRED = F
  Y     = IPRED * (1 + EPS(1))

$THETA
  (0.1,   3.0,  50.0)   ; TVCL
  (1.0,  10.0, 200.0)   ; TVV1
  (0.1,   5.0,  50.0)   ; TVQ2
  (1.0,  30.0, 500.0)   ; TVV2
  (0.1,   2.0,  20.0)   ; TVQ3
  (1.0,  20.0, 200.0)   ; TVV3

$OMEGA
  0.10   ; ETA_CL
  0.08   ; ETA_V1
  0.10   ; ETA_Q2
  0.08   ; ETA_V2

$SIGMA
  0.02   ; PROP_ERR

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5
$COVARIANCE PRINT=E
