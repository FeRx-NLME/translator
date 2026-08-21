$PROBLEM $ERROR reads a bare A(n), which NONMEM does not scale by S<n>
; Regression fixture for issue #32. The source declares `S2 = V` AND writes
; `Y = A(2)*(1+EPS(1))`. Those look like they belong together and do not:
; measured against NONMEM 7.6.0 (tests/nonmem-anchor/), `S<n>` is applied to `F`
; and NOT to a bare `A(n)` -- the ratio between the two spellings is exactly V at
; every timepoint. So the correct translation of this file emits NO [scaling]
; block, and the prediction is an amount, exactly as it is in the source.
;
; Before the fix this emitted `obs_scale = V`, dividing by a scale NONMEM never
; applied and putting every prediction low by a factor of 50, silently, in a file
; that validated and fit.
;
; DO NOT "fix" this model by changing $ERROR to `Y = F` or `Y = A(2)/S2`. Both of
; those ARE scaled and both are already covered -- by pk_1cmt_oral.mod and
; defobs_not_last.ctl respectively. The bare spelling is the whole point of this
; file, and it is the only one of the three that can show the defect.
;
; DO NOT add a Tier-4 concordance fit for this model either. Unscaled, CL and V
; are not separately identifiable: V enters the ODE only through CL/V, so amount
; data pins the ratio and nothing else, and a fit sits at its starting values
; regardless of whether the translation is right. It is checked by emitted text.
$INPUT ID TIME AMT DV MDV CMT
$DATA amount_readout.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9
$MODEL
  COMP=(DEPOT, DEFDOSE)
  COMP=(CENTRAL, DEFOBS)
$PK
  KA = THETA(1)
  CL = THETA(2)
  V  = THETA(3)
  S2 = V
$DES
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1) - (CL/V)*A(2)
$ERROR
  Y = A(2)*(1 + EPS(1))
$THETA
  (0, 1.0)   ; KA
  (0, 3.0)   ; CL
  (0, 50.0)  ; V
$OMEGA 0.09
$SIGMA 0.04
$EST METHOD=1
