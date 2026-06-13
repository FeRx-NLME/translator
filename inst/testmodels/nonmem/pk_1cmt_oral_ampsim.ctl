$PROB  1 CMT oral PK -- mirrors amp.sim PK.1CMT.ORAL IIV structure (ADVAN2)
;; ETA on KA and CL only; V is a fixed-effect theta with no IIV.
;; Used for amp.sim concordance tests; matches PK.1CMT.ORAL.mod structure
;; but uses ADVAN2 (linCmt path) to avoid the S2=V ODE scaling issue.

$INPUT  ID TIME DV AMT EVID MDV

$DATA   dummy.csv IGNORE=@

$SUBROUTINE ADVAN2 TRANS2

$PK
  KA = THETA(1) * EXP(ETA(1))
  CL = THETA(2) * EXP(ETA(2))
  V  = THETA(3)
  S2 = V

$ERROR
  IPRED = F
  Y     = IPRED * (1 + EPS(1))

$THETA
  (0, 0.1)   ; KA
  (0, 2.0)   ; CL
  (0, 1.0)   ; V

$OMEGA
  0.01   ; eta_KA
  0.02   ; eta_CL

$SIGMA
  0.1    ; prop_err

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999
