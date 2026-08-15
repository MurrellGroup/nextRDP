# Session 26 handoff — RDP-only cyclic parity audit

Version: `0.26.0-session-26`
Project schema: `org.rdp-web.project/v1alpha19`

This checkpoint starts from the Session 25 source archive and compares the browser core against ten
simulated, RDP-only saved analyses. Inputs were passed opaquely to the scanner; comparison code read
only WebRDP result JSON and the numeric saved-analysis tail after the RDP5 `Recombination Data`
marker. The supplied desktop source remained read-only and was not compiled. No alternate RDP
implementation was consulted.

## Source semantics restored

- RDP candidate breakpoints now follow the supplied `CentreBP` midpoint rule. The detector retains
  the first/last informative coordinates for scoring, then stores the midpoint between adjacent
  informative sites. After cyclic erasure, a midpoint that lands on missing data is shifted toward
  the detected tract exactly as in the supplied routine.
- The durable BestXOList-style action cache now covers every committed event, not only restored or
  manually fixed events. Ordinary cyclic reuse is deliberately narrow: the same original triplet
  and more than 80% tract overlap. This prevents the same acted-on signal from creating fragments
  indefinitely without hiding neighbouring events that merely share a pair.
- The supplied `AddjustCXO`/`MakePairsP` scheduling rule is active. Current XOverList signals that
  touch `RList(WinPP)` enable their three original-sequence pairs; those permissions are propagated
  to every selected RList member, and a dirty next-round triplet runs only when all three pairs are
  enabled. A sparse pair-key set replaces the native dense `DoPairs` matrix without changing the
  decision. Unchanged signals remain in the carried XOverList shortlist.
- `ModSeqNumY` erasure uses the selected role's final `ConsensusOK`/`FinalTrim` distance list. The
  diagnostic two-of-three evidence union remains visible but is no longer erased as though it were
  `RList(WinPP)`.
- The active `MakeConsensusC` subset quantizes PhPr-family and distance statistics to the supplied
  precisions before tie decisions, applies the exact post-rounding ±1 family sentinels and supplied
  tree-distance cap. The TreePhPr second-place contribution is the supplied 14 points, not a generic
  half of its 18-point win.
- Automated cyclic scanning has a fixed 64-round safety ceiling. Reaching it preserves all completed
  events, records `cycle-limit-reached`, and proceeds to Review/export instead of iterating without a
  diagnostic bound.

## Validation result and boundary

The midpoint rule is strongly confirmed: whenever an ordered event role matched in the focused
opening rounds, stored boundaries frequently became exact numeric matches, and the probability did
not change because centering is post-score. The saved source membership arrays also confirmed that
the selected role's final distance list—not the two-of-three union—is the erasure group for every
same-role event examined.

The ten-case audit is not full cyclic parity. Original RDP reports 336 events. Before the source
pair gate, this checkpoint reported 509 events, matched 82 ordered-role events and 152 unordered
triplets, reached the 64-round ceiling in five cases, and took 125.5 seconds. With the source
`DoPairs` gate it reports 453 events, matches 95 roles and 171 triplets, reaches the ceiling in one
case, and takes 72.7 seconds. Heavy triplet-kernel entries fall from 1,399,142 to 307,911 (78%),
while total audit time falls 42% because late event grouping/tree work remains substantial. The
gate improves all three aggregate parity measures, so it is retained as a general source invariant
rather than a fixture-specific adjustment. Exact `MakeConsensusC` role selection and source event
consolidation remain the next causal deficits.

An experimental partial `TestMoveInTreeAlt` acceptance gate was rejected from active behavior. It
reduced total events but worsened early ordered parity because source RDP's mutable
`PermValid`/`PermDiffs`, `vQuickDist`, and final tree stages must be ported as one state machine. The
helper trace remains non-authoritative and is not called by event selection.

## Recommended next checkpoint

1. Port the remaining active `MakeConsensusC` inputs (`dMax`, RCompat families, list correlations,
   outlier checks, and combined prizes) with the exact source role order and rounding.
2. Complete source XOverList/BestXOList event consolidation after the now-active `MakePairsP`
   scheduler, including how alternate retained records are redirected into existing events.
3. Port `TestMoveInTreeAlt` only with its mutable pair totals, `UFDist`, `vQuickDist`, and final
   Clearcut tree decisions together; then rerun all ten opaque comparisons.
4. Keep the 64-round cap during this work so a parity regression terminates with inspectable results.

Full native parity remains unclaimed.
