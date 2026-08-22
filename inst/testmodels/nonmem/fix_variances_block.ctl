;; Purpose: FIX on a $OMEGA BLOCK, end to end (issue #31, review of PR #38).
;;
;; The diagonal-omega and sigma cases live in fix_variances.ctl, which is also
;; fitted by the Tier-4 concordance test. This one is emitted-text only: adding a
;; correlated block to that model would change what it fits and cost the
;; discriminating fit its meaning.
;;
;; TWO blocks, one FIX and one not, because the discriminator is the pair. A
;; single fixed block cannot distinguish "FIX emitted from the source flag" from
;; "FIX appended to every block", and ferx fixes a block_omega all or nothing --
;; so the free block is what proves the flag is being read rather than assumed.
$PROB 2 blocks, one fixed and one free
$INPUT ID TIME DV EVID AMT CMT MDV
$DATA nonexistent.csv IGNORE=@
$SUBROUTINES ADVAN3 TRANS4

$PK
TVCL = THETA(1)
TVV1 = THETA(2)
TVQ  = THETA(3)
TVV2 = THETA(4)
CL = TVCL * EXP(ETA(1))
V1 = TVV1 * EXP(ETA(2))
Q  = TVQ  * EXP(ETA(3))
V2 = TVV2 * EXP(ETA(4))

$ERROR
Y = F * (1 + EPS(1))

$THETA
(0,5)   ; TVCL
(0,20)  ; TVV1
(0,8)   ; TVQ
(0,60)  ; TVV2

$OMEGA BLOCK(2) FIX
  0.09
  0.01 0.04
$OMEGA BLOCK(2)
  0.02
  0.005 0.03

$SIGMA
0.04

$EST METHOD=1 MAXEVAL=9999 NOABORT
