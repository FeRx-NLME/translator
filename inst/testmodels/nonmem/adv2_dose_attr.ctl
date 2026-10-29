$PROBLEM  One-compartment oral analytical model with bioavailability and lag time
$INPUT  ID TIME DV AMT EVID MDV CMT
$DATA dose_attr.csv IGNORE=@
$SUBROUTINE ADVAN2 TRANS2
$PK
KA = THETA(1)*EXP(ETA(1))
CL = THETA(2)
V = THETA(3)
S2 = V
F1 = THETA(4)
ALAG1 = THETA(5)
$ERROR
Y = F*(1+EPS(1))
$THETA (0,1) (0,2) (0,10) (0,0.6,1) (0,0.5)
$OMEGA 0.1
$SIGMA 0.04
$EST METHOD=1 INTER
