$PROBLEM 2-cpt IV, DEFDOSE on compartment 2, dataset has NO CMT column
;
; Fixture for issue #27. NM-TRAN resolves a dose row with no CMT data item to
; the $MODEL DEFDOSE compartment; ferx-core resolves it to compartment 1
; (`cmt_col.and_then(...).unwrap_or(1)` in src/io/datareader.rs, three sites),
; and has no model-side dose-compartment binding on the fitting path at all --
; `dose_cmt` exists only on SimulationSpec, parsed from [simulation]. So the
; two disagree by a whole compartment and the emitted file cannot say so.
;
; Measured on this model: doses land in PERIPH instead of CENTRAL, a factor of
; ~6 at t=1, with no warning from anything.
;
; BOTH halves are the subject under test and NEITHER alone shows it:
;
;   - Add a CMT column to $INPUT and the data decides, DEFDOSE is irrelevant,
;     and the code path is dead.
;   - Move DEFDOSE to compartment 1 and the two rules agree, so the wrong code
;     gives the right answer.
;
; DO NOT ADD CMT TO $INPUT AND DO NOT MOVE DEFDOSE. Either edit makes this file
; pass whether the guard works or not. See defdose_cmt_present.ctl for the
; negative control that holds DEFDOSE at 2 and adds the column back.
;
; The $MODEL order here already matches the $DES DADT order, so this is not
; issue #25: the compartment NUMBERING is correct and the dose still goes to
; the wrong place.
;
$INPUT ID TIME AMT DV MDV
$DATA defdose_no_cmt.csv IGNORE=@

$SUBROUTINE ADVAN6 TOL=6

$MODEL
  COMP=(PERIPH)
  COMP=(CENTRAL, DEFDOSE, DEFOBS)

$PK
  TVCL = THETA(1)
  TVV1 = THETA(2)
  CL = TVCL*EXP(ETA(1))
  V1 = TVV1
  Q  = THETA(3)
  V2 = THETA(4)
  S2 = V1

$DES
  DADT(1) = (Q/V1)*A(2) - (Q/V2)*A(1)
  DADT(2) = -(CL/V1)*A(2) - (Q/V1)*A(2) + (Q/V2)*A(1)

$ERROR
  Y = F*(1 + EPS(1))

$THETA (0, 5)    ; CL
$THETA (0, 20)   ; V1
$THETA (0, 8)    ; Q
$THETA (0, 60)   ; V2
$OMEGA 0.09      ; ETA_CL
$SIGMA 0.04

$ESTIMATION METHOD=1 INTER MAXEVAL=9999
