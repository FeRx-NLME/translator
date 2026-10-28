;; Purpose: $SIGMA BLOCK with a non-zero residual covariance (issue #39).
;;
;; Derived from pk_1cmt_oral.mod so it fits the SAME bundled dataset
;; (inst/testdata/ode_1cpt_oral_concordance.csv). The covariance is non-zero so
;; that "block kept" and "block flattened to its diagonal" are different models:
;; the two give different OFVs on this dataset.
$PROB 1 CMT oral, combined error with correlated residuals
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
(0,1)     ; V (l)

$OMEGA
.01     ; ETA KA
.02     ; ETA CL

$ERROR
Y = F * (1 + EPS(1)) + EPS(2)

$SIGMA BLOCK(2)
.04            ; proportional
.001 .01       ; covariance, additive

$EST METHOD=1 INTERACTION MAXEVAL=9999 NOABORT POSTHOC
