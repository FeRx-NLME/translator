$PROBLEM  Two-compartment oral PK with weight and CRCL covariates

$INPUT  ID TIME DV AMT EVID MDV WT CRCL

$DATA   2cpt_oral_cov.csv IGNORE=@

$SUBROUTINE ADVAN4 TRANS4

$PK
  TVCL     = THETA(1)
  TVV1     = THETA(2)
  TVQ      = THETA(3)
  TVV2     = THETA(4)
  TVKA     = THETA(5)
  THETA_WT   = THETA(6)
  THETA_CRCL = THETA(7)

  CL = TVCL * (WT/70)**THETA_WT * (CRCL/100)**THETA_CRCL * EXP(ETA(1))
  V2 = TVV1 * (WT/70)**THETA_WT * EXP(ETA(2))
  Q  = TVQ  * EXP(ETA(3))
  V3 = TVV2 * EXP(ETA(4))
  KA = TVKA * EXP(ETA(5))
  S2 = V2

$ERROR
  IPRED = F
  Y     = IPRED * (1 + EPS(1))

$THETA
  (0.1,   5.0,  100.0)  ; TVCL
  (1.0,  50.0,  500.0)  ; TVV1
  (0.1,  10.0,  100.0)  ; TVQ
  (1.0, 100.0,  500.0)  ; TVV2
  (0.01,  1.2,   10.0)  ; TVKA
  (0.01,  0.75,   5.0)  ; THETA_WT
  (0.01,  0.50,   5.0)  ; THETA_CRCL

$OMEGA
  0.10   ; ETA_CL
  0.10   ; ETA_V1
  0.05   ; ETA_Q
  0.05   ; ETA_V2
  0.15   ; ETA_KA

$SIGMA
  0.02   ; PROP_ERR

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5
$COVARIANCE PRINT=E
