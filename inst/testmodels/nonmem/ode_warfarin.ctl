$PROBLEM  Warfarin one-compartment oral ODE model

$INPUT  ID TIME DV AMT EVID MDV

$DATA   warfarin.csv IGNORE=@

$SUBROUTINE ADVAN6 TOL=6

$MODEL
  COMP=(DEPOT)
  COMP=(CENTRAL)

$PK
  TVCL = THETA(1)
  TVV  = THETA(2)
  TVKA = THETA(3)

  CL = TVCL * EXP(ETA(1))
  V  = TVV  * EXP(ETA(2))
  KA = TVKA * EXP(ETA(3))
  S2 = V

$DES
  DADT(1) = -KA * A(1)
  DADT(2) =  KA * A(1) / V - (CL / V) * A(2)

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

$SIGMA
  0.01   ; PROP_ERR

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5
$COVARIANCE PRINT=E
