$PROBLEM  One-compartment oral PK with block omega on CL and V

$INPUT  ID TIME DV AMT EVID MDV

$DATA   1cpt_oral.csv IGNORE=@

$SUBROUTINE ADVAN2 TRANS2

$PK
  TVCL = THETA(1)
  TVV  = THETA(2)
  TVKA = THETA(3)
  CL   = TVCL * EXP(ETA(1))
  V    = TVV  * EXP(ETA(2))
  KA   = TVKA * EXP(ETA(3))
  S2   = V

$ERROR
  IPRED = F
  Y     = IPRED * (1 + EPS(1))

$THETA
  (0.001, 0.134, 10.0)  ; TVCL
  (0.1,   8.1,  500.0)  ; TVV
  (0.01,  1.0,   50.0)  ; TVKA

$OMEGA BLOCK(2)
  0.07         ; ETA_CL
  0.02  0.02   ; ETA_V (off-diagonal = covariance)

$OMEGA
  0.40   ; ETA_KA

$SIGMA
  0.01   ; PROP_ERR

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5
$COVARIANCE PRINT=E
