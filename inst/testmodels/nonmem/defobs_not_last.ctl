$PROBLEM 2-cpt IV, DEFOBS declared FIRST
;
; Regression fixture for review finding 7. Every other bundled model declares
; its observed compartment LAST, so tail(state_names, 1) guessed correctly and
; the translator's inability to read DEFOBS was invisible across the whole
; corpus. Here DEFOBS is compartment 1 and PERIPH is compartment 2.
;
; DO NOT REORDER THE $MODEL BLOCK. The declaration order IS the subject under
; test: with DEFOBS first, the old positional guess picked PERIPH, which both
; attributed observations to the peripheral compartment AND made the S1 = V
; scaling below unreachable -- the compartment number for the scaling lookup is
; derived from the same guess, so `S1 = V` parsed correctly and was then
; silently discarded. That is the S2=V silent-divergence class CLAUDE.md warns
; about, reached by a different route.
;
$INPUT ID TIME AMT DV MDV
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=6
$MODEL
  COMP=(CENT, DEFDOSE, DEFOBS)
  COMP=(PERIPH)
$PK
  CL = THETA(1)*EXP(ETA(1))
  V  = THETA(2)
  Q  = THETA(3)
  V2 = THETA(4)
  S1 = V
$DES
  DADT(1) = -CL/V*A(1) - Q/V*A(1) + Q/V2*A(2)
  DADT(2) =  Q/V*A(1) - Q/V2*A(2)
$ERROR
  Y = A(1)/S1*(1+EPS(1))
$THETA (0,5) (0,20) (0,8) (0,60)
$OMEGA 0.09
$SIGMA 0.04
$ESTIMATION METHOD=1 INTERACTION
