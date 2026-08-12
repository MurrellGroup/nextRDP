# Supplied-source cyclic shortlist and correction trace

## Scope

This trace uses only the supplied RDP5 VB and DLL source. It covers the signal summaries and
shortlists Darren identified (`XOverList` / `XOverDefine`, `BestXOList`, `Worthwhilescan`, and
`StoreLPV`) and the multiple-comparison factor used during cyclic rescanning. It does not use an
alternate RDP implementation. No alternate RDP implementation was consulted.

## What the supplied code retains

- `Module31.bas::XOverDefine` stores the method, provisional roles, beginning/end, probability,
  event number, acceptance state, and auxiliary breakpoint/probability fields. `XOverList` is
  therefore a compact detected-signal catalogue rather than a copy of each triplet profile.
- `Module2.bas::MakeBestXOList` and `MakeNextBestXOLists` copy actionable `CollectEvents` records
  into per-sequence `BestXOList` rows. They retain the main signal, grow in bounded increments, and
  prefer replacing a weaker same-method non-main entry when storage cannot grow.
- `MathFuncsDll.cpp::MakeAListIS*` consults `Worthwhilescan` plus method bits in `ProgBinRead` before
  adding a triplet to a method analysis list. `AlistRDP4`, `AlistGC2`, `AlistMC2`, and `AlistChi`
  combine the shortlist with `StoreLPV` thresholds and mark only selected entries in `RedoL3` for
  the detailed VB method call.
- `Module3.bas::FindBetterRecSignal` starts from the three members and probability of an existing
  `XOverList` signal. It marks replacements already present in `AnalysisListX` or `BanTriplet`, then
  creates only new triplets that retain two anchor members and substitute the third. Its threshold
  is the probability of the signal it is trying to improve.

These structures have two related jobs: retain summaries of useful signals and avoid repeating
method work that has already been shown irrelevant or acted upon.

## Project correction lifetime

`Module2.bas::BuildFirstXOList` resets `NextNo` to `PermNextno` and calls `MakeMCCorrection` before
the event-reconstruction loop. For an ordinary exploratory analysis, the active factor is the
initial unmasked `choose(NS, 3)` count, capped at `(255^4)/2`. The query/reference branch obtains its
group-pair × query factor from `MakeAnalysisListQvR`. It is not recalculated from each later
fragment-expanded working schedule.

WebRDP now keeps two counts deliberately separate:

- `totalTriplets` is the current round's real working schedule and can change after erasure or
  fragment re-entry;
- `correctionTests` is the initial source scan-plan factor and remains fixed for every cyclic
  detection and post-group recheck.

This also makes cached probabilities safe to reuse: the sequence bytes and the probability factor
must both be unchanged.

## Browser-safe shortlist implementation

Materializing a desktop-sized entry for every empty triplet would waste browser memory. The WASM
scanner instead retains only signal-bearing per-working-triplet summaries and a compact dirty-row
set:

1. Before mapping fragment rows back to original sequence identities, each emitted signal records
   the exact three mutable working rows that produced it.
2. After an event is selected and erased, every row whose bytes changed and every newly created
   fragment row is marked dirty.
3. A next-round triplet containing a dirty row runs all enabled discovery kernels from fresh
   sequence state.
4. A triplet containing no dirty row replays its retained signal summary. The absence of a retained
   summary is also definitive: the same triplet and same fixed correction factor previously
   produced no threshold-passing signal, so its stable method kernels can be skipped.
5. Support-bearing summaries remain eligible when their exact working rows were not modified. The
   selected event's durable support records remain separate; replay simply reproduces what a fresh
   scan of the unchanged bytes would emit. Replayed summaries re-enter in current schedule order,
   preserving strongest-first deterministic tie behavior.

The one method boundary is 3SEQ. Its first post-erasure pass enables the supplied
`FindSubSeqTS2` / `CheckSplit3Seq` behavior, so first-round 3SEQ summaries are not reused. Unchanged
triplets rerun 3SEQ once while reusing the other stable methods. Later rounds may reuse 3SEQ too,
because the post-erasure mode is then unchanged.

This is an exact state-based reuse rule, not a probabilistic cache. Manual correction/rejection
rebuilds clear the shortlist and begin with a full fresh pass.

## Diagnostics and regression boundary

Progress now reports `tripletKernelEvaluations`, `tripletSummariesReused`, `cachedSignalsReused`,
and `methodScansSkipped` separately from scheduled/cumulative triplets. The production Pages
verification includes a deterministic two-event WASM scan that asserts:

- at least two cyclic events and three rounds;
- an unchanged initial correction factor despite fragment-expanded later schedules;
- cached signal replay and skipped method work; and
- fewer kernel evaluations than scheduled cumulative triplets.

On that fixture, 476 scheduled triplets require 280 triplet-kernel evaluations; 196 unchanged
triplet summaries are reused and 196 RDP method scans are skipped. With all five current discovery
families enabled on the same fixture, 131,488 scheduled triplets require 32,911 triplet-kernel
evaluations, with 98,661 unchanged summary reuses, 1,199 cached signals replayed, and 493,221
individual method scans skipped.

A temporary full-rescan build ran all 476 RDP kernels on the two-event fixture. The optimized and
full-rescan builds produced identical event roles, breakpoints, detection rounds, support-signal
order, signal triplets, and p-values. The complete selected-result digest matched across all 22
events and 422 retained signals. This comparison is a scheduler regression, not a substitute for
authorized native RDP5 golden output.

## Remaining fidelity boundary

WebRDP still enumerates the lightweight current-round triplet schedule so progress, cancellation,
query/reference constraints, same-origin fragment exclusion, and deterministic order remain exact.
It does not yet replace that enumeration with the narrower `FindBetterRecSignal` two-anchor
substitution list. The expensive alignment/profile/probability kernels are skipped safely; removing
the remaining `O(T)` scheduler walk should wait for native golden fixtures covering event ties,
masked rows, query/reference groups, fragment aliases, and `BanTriplet` behavior.
