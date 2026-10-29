$PROBLEM 2-cpt IV, $INPUT names CMT, CSV header agrees
;
; Fixture for issue #36. Negative control: the header names the column CMT, so NONMEM and ferx read the same column.
;
; The two twins differ in exactly one thing, the header of the CSV beside them
; (cmt_header_match.csv): the $INPUT is identical and names a CMT item in both. NONMEM
; binds columns by $INPUT position and ferx by CSV header name, so only the
; header decides whether ferx finds the CMT column. DO NOT EDIT THE HEADER
; OR $INPUT: either edit makes the pair stop discriminating.
;
$INPUT ID TIME AMT DV MDV CMT
$DATA cmt_header_match.csv IGNORE=@

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
