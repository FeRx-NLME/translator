$PROBLEM  Two-compartment IV bolus pharmacokinetics

$INPUT  ID TIME DV AMT EVID MDV

$DATA   2cpt_iv.csv IGNORE=@

$SUBROUTINE ADVAN3 TRANS4

$PK
  TVCL = THETA(1)
  TVV1 = THETA(2)
  TVQ  = THETA(3)
  TVV2 = THETA(4)

  CL = TVCL * EXP(ETA(1))
  V1 = TVV1 * EXP(ETA(2))
  Q  = TVQ  * EXP(ETA(3))
  V2 = TVV2 * EXP(ETA(4))
  S1 = V1

$ERROR
  IPRED = F
  Y     = IPRED * (1 + EPS(1))

$THETA
  (0.1,   5.0,  100.0)  ; TVCL
  (1.0,  20.0,  500.0)  ; TVV1
  (0.1,   8.0,  100.0)  ; TVQ
  (1.0,  60.0,  500.0)  ; TVV2

$OMEGA
  0.10   ; ETA_CL
  0.10   ; ETA_V1
  0.08   ; ETA_Q
  0.05   ; ETA_V2

$SIGMA
  0.02   ; PROP_ERR

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5
$COVARIANCE PRINT=E
