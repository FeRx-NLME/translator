$PROBLEM $MODEL declares a compartment with no DADT, so d/dt order is not COMP order
; Regression fixture for issue #25. ferx numbers compartments by POSITION in
; `states=[...]` -- measured on 0.3.0, a dose row's CMT is applied as
; `u[CMT - 1] += F*AMT` straight into the state vector. NONMEM numbers by $MODEL
; COMP order. nonmem2rx hands back $DES statement order, which is neither.
;
; DUMMY is the whole point of this file. It is declared as compartment 2 and has
; no DADT, so nonmem2rx materialises it as `d/dt(DUMMY) = 0` and places it
; FIRST -- yielding states [DUMMY, DEPOT, CENTRAL] against COMP
; [DEPOT, DUMMY, CENTRAL]. NONMEM's CMT=1 dose then landed in DUMMY, whose
; derivative is zero, and every prediction came back 0.0 from a file that
; validated clean and warned about nothing.
;
; DO NOT "tidy" DUMMY away, and do not give it a DADT. Note that $DES is written
; in ascending index order: this fixture needs NO reordered $DES, which is what
; makes it the realistic trigger rather than a curiosity. A $MODEL compartment
; without a differential equation is ordinary NONMEM.
;
; DO NOT move DUMMY to the end either. A TRAILING gap is a different defect
; (issue #26): nonmem2rx drops the compartment and loses every state name, so
; the COMP list cannot be matched by name at all and the translator declines
; instead of renumbering. That path is covered by an inline fixture in
; test-integration.R, not here.
$INPUT ID TIME AMT DV MDV CMT
$DATA cmt_order_gap.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9
$MODEL
  COMP=(DEPOT, DEFDOSE)
  COMP=(DUMMY)
  COMP=(CENTRAL, DEFOBS)
$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
  KA = THETA(3)
; S3, not S2: the scaled compartment is CENTRAL, which is $MODEL compartment 3
; and d/dt position 3 only after the renumbering. Reading it as a position gave
; the third state of [DUMMY, DEPOT, CENTRAL], which is CENTRAL by coincidence
; here -- so this file alone cannot show the obs_scale half of #25. The
; reordered-$DES fixture in test-integration.R is the one that separates those.
  S3 = V
$DES
  DADT(1) = -KA*A(1)
  DADT(3) =  KA*A(1) - (CL/V)*A(3)
$ERROR
; `F`, not a bare `A(3)`, and that is load-bearing twice over.
;
; NONMEM applies `S<n>` to F and NOT to a bare A(n) -- anchored against NONMEM
; 7.6.0 in tests/nonmem-anchor/, ratio exactly V at every timepoint. So with
; `Y = A(3)` this model would be unscaled, and then CL and V are not separately
; identifiable: V enters the ODE only as CL/V, so amount data determines the
; ratio and nothing more, and the concordance fit sits at its starting values.
; The scaled prediction A(3)/V is what makes V identifiable at all.
;
; Change this line and the Tier-4 test below it stops meaning anything. The bare
; A(n) case has its own fixture -- see amount_readout.ctl -- which is checked by
; simulation rather than by fitting, for exactly this reason.
  Y = F*(1 + EPS(1))
$THETA
  (0, 3)     ; CL
  (0, 50)    ; V
  (0, 1.2)   ; KA
$OMEGA 0.09
$SIGMA 0.04
$EST METHOD=1
