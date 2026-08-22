;; Purpose: FIX on every kind of variance parameter (issue #31).
;;
;; Derived from pk_1cmt_oral.mod so it fits the SAME bundled dataset
;; (inst/testdata/ode_1cpt_oral_concordance.csv), which lets the concordance
;; test assert what a fixed parameter DOES rather than only what is emitted.
;;
;; Every initial below is chosen so that "held at its start" and "estimated"
;; are different observations:
;;
;;   ETA_KA  .01  free   -- simulation truth; the control that must MOVE
;;   ETA_CL  .05  FIX    -- truth is .02, and a free fit of this dataset lands
;;                          on .0197 (measured, ferx 0.3.0). Held, it must read
;;                          exactly .05. An initial AT the truth could not tell
;;                          "fixed" from "recovered it" -- see CLAUDE.md's
;;                          discriminating-fixture rule.
;;   SIGMA   .3   FIX    -- truth is .1. Emitted on the SD scale as sqrt(.3).
;;   V       1    FIX    -- theta FIX, which already worked; here as the control
;;                          that the theta path is untouched by this change.
$PROB 1 CMT oral, fixed variance parameters
$INPUT ID TIME DV EVID AMT CMT MDV
$DATA ode_1cpt_oral_concordance.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=3
$MODEL
COMP = (ABS)
COMP = (CENTRAL)

$PK
KA = THETA(1) * EXP(ETA(1))
CL = THETA(2) * EXP(ETA(2))
V  = THETA(3)
S2  = V
K20 = CL/V

$DES
DADT(1) = - KA*A(1)
DADT(2) =   KA*A(1) - K20*A(2)

$THETA
(0,.1)    ; KA (1/h)
(0,2)     ; CL (l/h)
(0,1) FIX ; V (l), fixed

$OMEGA
.01     ; ETA KA, estimated
.05 FIX ; ETA CL, fixed away from the simulation truth of .02

$ERROR
Y = F * (1 + EPS(1))

$SIGMA
.3 FIX ; proportional error, fixed away from the simulation truth of .1

$EST METHOD=1 INTERACTION MAXEVAL=9999 NOABORT POSTHOC
