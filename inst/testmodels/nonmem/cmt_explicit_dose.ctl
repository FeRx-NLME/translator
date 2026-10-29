$PROBLEM 2-cpt IV, DEFDOSE on compartment 2, CMT column present, every dose row names its compartment
;
; Fixture for issue #35. Negative control: every dose row carries a compartment, so NONMEM and ferx agree.
;
; The pair differs in exactly one cell of the CSV beside it, the CMT of the
; first dose row (0 vs 2); the control stream is identical. Observation rows
; carry CMT 2 in both, so a check that counted any CMT 0 cell, rather than dose
; rows only, would still be right here and wrong on a real dataset.
; DO NOT EDIT THE CSV: either edit makes the pair stop discriminating.
;
$INPUT ID TIME AMT DV MDV CMT
$DATA cmt_explicit_dose.csv IGNORE=@

$SUBROUTINE ADVAN6 TOL=6

$MODEL
  COMP=(PERIPH)
  COMP=(CENTRAL, DEFDOSE, DEFOBS)

$PK
  TVCL = THETA(1)
  TVV1 = THETA(2)
  CL = TVCL*EXP(ETA(1))
  V1 = TVV1
  TVQ  = THETA(3)
  Q    = TVQ
  TVV2 = THETA(4)
  V2   = TVV2
  S2 = V1

$DES
  DADT(1) = (Q/V1)*A(2) - (Q/V2)*A(1)
  DADT(2) = -(CL/V1)*A(2) - (Q/V1)*A(2) + (Q/V2)*A(1)

$ERROR
  Y = F*(1 + EPS(1))

$THETA (0, 5)    ; TVCL
$THETA (0, 20)   ; TVV1
$THETA (0, 8)    ; TVQ
$THETA (0, 60)   ; TVV2
$OMEGA 0.09      ; ETA_CL
$SIGMA 0.04

$ESTIMATION METHOD=1 INTER MAXEVAL=9999
