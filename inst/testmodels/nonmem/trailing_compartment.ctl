; DO NOT REORDER. The compartment with no DADT must be the LAST one declared in
; $MODEL, and $DES must write DADT(2) before DADT(1). A gap anywhere else keeps
; its names in nonmem2rx and exercises cmt_order_gap.ctl instead; a $DES in
; compartment order cannot tell the DADT index from the statement position.
$PROBLEM trailing $MODEL compartment with no DADT
$INPUT ID TIME DV AMT EVID MDV CMT
$DATA trailing_compartment.csv IGNORE=@
$SUBROUTINE ADVAN13 TOL=6
$MODEL
  COMP=(DEPOT, DEFDOSE)
  COMP=(CENTRAL, DEFOBS)
  COMP=(SPARE)
$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
  KA = THETA(3)
  S2 = V
$DES
  DADT(2) =  KA*A(1) - (CL/V)*A(2)
  DADT(1) = -KA*A(1)
$ERROR
  IPRED = A(2)/S2
  Y = IPRED*(1 + EPS(1))
$THETA (0,5) (0,50) (0,1)
$OMEGA 0.09
$SIGMA 0.04
$ESTIMATION METHOD=1
