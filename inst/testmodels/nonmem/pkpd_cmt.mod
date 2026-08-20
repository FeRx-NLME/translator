;; Description: 1-cpt PK + indirect response PD, two endpoints dispatched on CMT
;; Exercises the per-CMT readout path: ferx does not expose CMT as a covariate,
;; so a source switching on it must become y[CMT=N] plus CMT=N: error bodies.
$PROBLEM PKPD dispatched on CMT
$INPUT ID TIME AMT EVID MDV CMT DV
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9
$MODEL
  COMP=(CENT, DEFDOSE, DEFOBS)
  COMP=(PD)
$PK
  KEL  = THETA(1)*EXP(ETA(1))
  VC   = THETA(2)*EXP(ETA(2))
  KOUT = THETA(3)
  BASE = THETA(4)*EXP(ETA(3))
  KIN  = BASE*KOUT
  A_0(2) = BASE
$DES
  CP = A(1)/VC
  DADT(1) = -KEL*A(1)
  DADT(2) =  KIN*(1.0 + CP) - KOUT*A(2)
$ERROR
  CONC = A(1)/VC
  RESP = A(2)
  IPRED = CONC
  IF (CMT.EQ.2) IPRED = RESP
  W1 = 0
  W2 = 0
  IF (CMT.EQ.1) W1 = 1
  IF (CMT.EQ.2) W2 = 1
  Y = IPRED*(1 + W1*EPS(1)) + W2*EPS(2)
$THETA
  (0.001, 0.1, 1.0)    ; 1 KEL
  (0.5, 5.0, 50.0)     ; 2 VC
  (0.01, 0.2, 5.0)     ; 3 KOUT
  (1.0, 10.0, 100.0)   ; 4 BASE
$OMEGA
  0.09
  0.04
  0.04
$SIGMA
  0.0225
  1.0
$ESTIMATION METHOD=1 INTER MAXEVAL=9999 PRINT=5 NOABORT
