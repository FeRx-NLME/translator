;; Source: amp.sim package (LeidenAdvancedPKPD/amp.sim), inst/example_models/PKPD.IR.mod
;; Description: 2-cpt PK + indirect response PD (ADVAN6, 4-state ODE)
$PROBLEM  indirect response model Effect compartment initialized with unity dose
$DATA dat IGNORE=#
$INPUT ID TIME DV MDV AMT BW DOSE CMT iCL
$SUBROUTINES ADVAN6 TOL=5
$MODEL
COMP=(DOSE, DEFDOSE)
COMP=(CENTRAL)
COMP=(PERIPH)
COMP=(EFFECT,DEFOBS)

$PK
K12= THETA(1)
CL = THETA(2)
V2 = THETA(3)
V3 = THETA(4)
Q  = THETA(5)
K20=CL/V2
K23=Q/V2
K32=Q/V3

BL  = THETA(6)*EXP(ETA(1))
Kout= THETA(7)
Rin = BL*Kout
A_0(4)=BL

EMAX= THETA(8)
EC50= THETA(9)

$DES
DADT(1) = -K12*A(1)
DADT(2) =  K12*A(1) - K20*A(2) - K23*A(2) + K32*A(3)
DADT(3) =                        K23*A(2) - K32*A(3)

C2  = A(2)/V2
EFF =  1 - C2*EMAX/( C2 + EC50 )
DADT(4)= EFF*Rin - Kout*A(4)

$ERROR
IPRE = F
Y    = IPRE + ERR(1)

$THETA
(1.95   )   ; TH1 K12
(1      )   ; TH2 CL
(1.31   )   ; TH3 V2
(4.15   )   ; TH4 V3
(0.904  )   ; TH5 Q

(0 100  )   ; TH6 BL
(0 3.94 )   ; TH7 Kout
(0 0.9 1)   ; TH8 EMAX
(0 2    )   ; TH9 EC50

$OMEGA
0.1         ; ETA1 = BL

$SIGMA
30          ; ERR1

$EST PRINT=5 MAX=9999 METHOD=0 POSTHOC
$COV COMP
$TABLE ID TIME DOSE CMT BL EC50 C2 EFF IPRE
ONEHEADER NOPRINT FILE=par
