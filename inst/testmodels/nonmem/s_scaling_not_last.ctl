$PROBLEM 2-cpt oral, observed compartment named only by $PK's S2
; The observed compartment is CENTRAL (compartment 2), and the ONLY thing in
; this control stream that says so is `S2 = V`. There is no DEFOBS attribute and
; $ERROR goes through NONMEM's bare F, which names no compartment.
;
; DO NOT REORDER the COMP lines and do not add a DEFOBS. PERIPH is declared last
; on purpose: it is what makes this fixture able to fail. With the observed
; compartment declared last, `tail(states)` and the right answer coincide and
; the fixture cannot tell the S<n> tier from the positional guess it replaced --
; which is true of every other bundled ODE model, and is why this one exists.
$INPUT  ID TIME DV AMT EVID MDV
$DATA   d.csv IGNORE=@

$SUBROUTINE ADVAN6 TOL=6

$MODEL
  COMP=(DEPOT)
  COMP=(CENTRAL)
  COMP=(PERIPH)

$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
  KA = THETA(3)
  Q  = THETA(4)
  V3 = THETA(5)
  S2 = V

$DES
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1) - (CL/V)*A(2) - (Q/V)*A(2) + (Q/V3)*A(3)
  DADT(3) =  (Q/V)*A(2) - (Q/V3)*A(3)

$ERROR
  IPRED = F
  Y     = IPRED*(1 + EPS(1))

$THETA
  (0, 5)    ; CL
  (0, 50)   ; V
  (0, 1)    ; KA
  (0, 8)    ; Q
  (0, 60)   ; V3
$OMEGA 0.09
$SIGMA 0.04
$EST METHOD=1 INTERACTION MAXEVAL=9999
