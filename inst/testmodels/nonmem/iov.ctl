$PROBLEM  One-compartment oral PK with inter-occasion variability

$INPUT  ID TIME DV AMT EVID MDV OCC

$DATA   1cpt_oral_iov.csv IGNORE=@

$SUBROUTINE ADVAN2 TRANS2

$PK
  TVCL = THETA(1)
  TVV  = THETA(2)
  TVKA = THETA(3)

  KAP_CL = ETA(4)

  CL = TVCL * EXP(ETA(1) + KAP_CL)
  V  = TVV  * EXP(ETA(2))
  KA = TVKA * EXP(ETA(3))
  S2 = V

$ERROR
  IPRED = F
  Y     = IPRED * (1 + EPS(1))

$THETA
  (0.001, 0.134, 10.0)  ; TVCL
  (0.1,   8.1,  500.0)  ; TVV
  (0.01,  1.0,   50.0)  ; TVKA

$OMEGA
  0.07   ; ETA_CL
  0.02   ; ETA_V
  0.40   ; ETA_KA

$OMEGA SAME   ; IOV on CL (occasion-specific, reuses structure)
  0.04        ; KAPPA_CL -- this line varies by implementation

$SIGMA
  0.01   ; PROP_ERR

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5
$COVARIANCE PRINT=E
