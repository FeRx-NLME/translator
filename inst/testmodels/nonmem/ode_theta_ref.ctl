$PROBLEM $DES names a THETA directly, with no $PK assignment to carry it
; Regression fixture for issue #6 defect 2. A theta is not in scope in a ferx
; [odes] block -- ferx resolves a d/dt right-hand side against states, individual
; parameters and covariates only -- so a theta named there has to be carried in by
; an individual parameter.
;
; KTP is referenced from $DES and never assigned in $PK. nonmem2rx does not bind a
; theta referenced by its $THETA label: it leaves a free `KTP` in the expression
; and records the theta in an `rxmissingvars` placeholder, so without the carrier
; the emitted [odes] block names an identifier declared nowhere. Measured against
; ferx 0.3.0, that IS caught -- "[odes]: RHS references undefined name(s): KTP",
; with or without a dataset -- so this fixture fails the concordance corpus sweep
; and aborts to_ferx(strict = TRUE) if the carrier ever stops being emitted.
;
; KSS is the negative case and must be left alone: it is read from $PK, which
; becomes [individual_parameters], where a theta IS in scope.
$INPUT ID TIME AMT DV MDV
$DATA ode_theta_ref.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9
$MODEL
  COMP=(CENT, DEFDOSE)
  COMP=(PERIPH, DEFOBS)
$PK
  KEL = THETA(1)*EXP(ETA(1))*KSS
$DES
; The intermediate touches a state, so it is inlined into the equation rather than
; emitted as a parameter of its own -- the path that bypassed the [odes] scope and
; left a bare theta in the output even once the carrier existed.
  FLUX    = KTP*A(1)
  DADT(1) = -KEL*A(1) - FLUX
  DADT(2) =  FLUX
$ERROR
  Y = A(2)*(1+EPS(1))
$THETA
  (0,0.05)   ; KEL
  (0,0.2)    ; KTP
  (0,2.0)    ; KSS
$OMEGA 0.09
$SIGMA 0.04
$EST METHOD=1
