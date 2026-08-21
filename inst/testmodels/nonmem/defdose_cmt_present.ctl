$PROBLEM 2-cpt IV, DEFDOSE on compartment 2, dataset HAS a CMT column
;
; Negative control for issue #27, and the twin of defdose_no_cmt.ctl. Identical
; in every respect except the CMT data item on $INPUT.
;
; With a CMT column the dataset says which compartment each dose reaches and
; NM-TRAN never consults DEFDOSE, so NONMEM and ferx agree and the translator
; must say nothing. This file is what shows the #27 guard is not simply firing
; on DEFDOSE-not-1: widen it to drop the CMT half and this model starts
; reporting a divergence that does not exist.
;
; DO NOT REMOVE CMT FROM $INPUT. That turns this into defdose_no_cmt.ctl and
; leaves the guard with no negative control at all.
;
$INPUT ID TIME AMT DV MDV CMT
$DATA defdose_cmt_present.csv IGNORE=@

$SUBROUTINE ADVAN6 TOL=6

$MODEL
  COMP=(PERIPH)
  COMP=(CENTRAL, DEFDOSE, DEFOBS)

$PK
  TVCL = THETA(1)
  TVV1 = THETA(2)
  CL = TVCL*EXP(ETA(1))
  V1 = TVV1
  Q  = THETA(3)
  V2 = THETA(4)
  S2 = V1

$DES
  DADT(1) = (Q/V1)*A(2) - (Q/V2)*A(1)
  DADT(2) = -(CL/V1)*A(2) - (Q/V1)*A(2) + (Q/V2)*A(1)

$ERROR
  Y = F*(1 + EPS(1))

$THETA (0, 5)    ; CL
$THETA (0, 20)   ; V1
$THETA (0, 8)    ; Q
$THETA (0, 60)   ; V2
$OMEGA 0.09      ; ETA_CL
$SIGMA 0.04

$ESTIMATION METHOD=1 INTER MAXEVAL=9999
