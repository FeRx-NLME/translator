$PROBLEM compartment name collides with an $ERROR variable
; Regression fixture for issue #6 defect 1. $MODEL names a compartment RTOT and
; $ERROR also assigns a variable RTOT, so nonmem2rx renames the compartment to
; `c.RTOT` -- a name ferx cannot parse at any reference site. The translator has
; to sanitise it, and the corpus sweep in test-integration.R asserts that every
; emitted identifier is legal, which without this model would hold vacuously.
$INPUT ID TIME AMT DV MDV
$DATA dotted_state.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9
$MODEL
  COMP=(CENT, DEFDOSE)
  COMP=(RTOT, DEFOBS)
$PK
  KEL  = THETA(1)*EXP(ETA(1))
  KDEG = THETA(2)
$DES
; the second state is referenced from the FIRST equation on purpose: the
; original defect renamed the declaration but left this inlined cross-state
; reference behind, so a fixture where each equation touches only itself would
; not exercise it.
  DADT(1) = -KEL*A(1) - A(2)
  DADT(2) = -KDEG*A(2)
$ERROR
  RTOT = A(2)
  Y = A(2)*(1+EPS(1))
$THETA (0,0.05) (0,0.3)
$OMEGA 0.09
$SIGMA 0.04
$EST METHOD=1
