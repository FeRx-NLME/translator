$PROB  1 CMT oral PK with an M3 censored-data likelihood (BLQ records carry the limit in DV)
;; Initials sit away from the simulation truth (CL 2, V 20, KA 1, SD 0.2,
;; omega 0.09) so "recovered" and "never moved" are distinguishable.

$INPUT  ID TIME DV AMT EVID MDV CMT CENS

$DATA   dummy.csv IGNORE=@

$SUBROUTINE ADVAN2 TRANS2

$PK
  CL = THETA(1) * EXP(ETA(1))
  V  = THETA(2)
  KA = THETA(3)
  S2 = V

$ERROR
  IPRED  = F
  W      = IPRED * THETA(4)
  LLOQ   = 0.3
  F_FLAG = 0
  Y      = IPRED + W * EPS(1)
  IF (CENS.EQ.1) THEN
    F_FLAG = 1
    Y      = PHI((LLOQ - IPRED) / W)
  ENDIF

$THETA
  (0, 3.0)   ; CL
  (0, 30.0)  ; V
  (0, 2.0)   ; KA
  (0, 0.3)   ; proportional SD

$OMEGA
  0.2

$SIGMA
  1 FIX

$ESTIMATION METHOD=1 INTERACTION LAPLACIAN MAXEVAL=9999
