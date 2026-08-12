# Session 4 handoff

This is a source-only checkpoint. Nothing in the supplied desktop source or browser port was
compiled, bundled, served, previewed, type-checked, or executed.

## Delivered in this phase

1. Added `phylogeny.cpp`/`.hpp` as a standalone browser-core layer. It computes pairwise
   Jukes–Cantor distances from valid aligned states, represents insufficient/saturated pairs with
   distance 10, constructs deterministic neighbour-joining trees, and calculates leaf patristic
   matrices.
2. Added deterministic column bootstrap with the event-path native count of ten replicates.
   Internal branches are matched by canonical bipartition and branches with support below 50% are
   collapsed to zero length before the second patristic matrix is generated.
3. Added the manual's six event subalignments: two regions around each breakpoint plus the complete
   outside and inside tract partitions. Tree diagnostics retain site count, sequence count,
   replicate count, supported/total internal branches, and usability for every region.
4. Added the source-shaped paired-tree membership test for every presumed-recombinant role. A
   candidate-anchor pair must be closer than all four anchor/candidate-to-parent comparisons in
   both trees of at least one region pair. Raw and collapsed margins are exported separately.
5. Bounded the cubic NJ panel at 100 active sequences, always retaining all event anchors and then
   the most similar active sequences. Larger active sets use cached Jukes–Cantor affinity for
   omitted candidates; every fallback row is explicit in UI/JSON. Masked sequences stay out of
   phylogenetic construction.
6. Expanded detectable sets with the supplied `FindSets` cross-role closure: membership in the two
   other role sets adds the candidate to the remaining role, repeated to stability.
7. Combined detectable, distance, and phylogenetic membership with the manual's general
   at-least-two-of-three rule. Each role hypothesis and event now carries a complete
   co-recombinant group instead of the session-3 lower-bound intersection.
8. Ported the three category swaps and two cyclic relabellings from `CalCR`. Direct and selected
   correlations, native inversion class, and exact six-value P are retained.
9. Ported the dominant-category `RCorrWarn` gate and the positive warning-adjusted aggregate score
   from `MakeRList`. `MakeGoodC` overlap and inverse `r > 0.83` support remain visible per candidate.
10. Added an auditable role recommendation. `PhPr`, `TreePhPr`, `SubDist`, `TreeSubDist`, triplet
    topology change, and three-set support each contribute an equal vote when informative. The UI
    shows every winner and vote margin and requires an explicit apply action.
11. Added final primary-use-case exports. The tract-masked FASTA replaces accepted event tracts with
    gaps in the current co-recombinant group. The mosaic FASTA performs the same operation in event
    order and appends aligned fragment-only records. Buttons remain locked until every event is
    reviewed and no correction awaits downstream reconciliation.
12. Advanced project JSON to `v1alpha4`, retained `v1alpha1`–`v1alpha3` import, expanded CSV/JSON
    diagnostics, and updated the responsive review/export interface.

## Performance shape

- The six tree families are built once per event and reused by all three roles.
- JC matrices load each alignment column once and update all valid pairs from compact byte states.
- Candidate-to-anchor reference distances are cached once per event region.
- NJ work is bounded at 100 sequences; no unbounded `O(N³)` tree is attempted for large inputs.
- Primary scanning remains in 512-triplet yielding batches; reconciliation/tree work stays in the
  dedicated worker and never moves to the DOM thread.

## Fidelity boundary

- The role vote is source shaped but equally weighted. The native application uses a larger method
  battery and native weights, so `nativeWeightParity` is false.
- The positive `MakeRList` aggregate and inversion/warning paths are mapped, but exact
  `AcceptableCoR`, contradictory-inversion cleanup, `StripDupInv`, `ConsensusOK`, and `FinalTrim`
  parity remains.
- The >100-sequence JC fallback is a declared browser performance adaptation rather than native
  phylogenetic placement.
- Ordered event rebuilding is implemented, but has not been proven identical to native sequence
  tract erasure, fragment insertion, and repeated scan termination.
- The alignment and tree evidence are numeric/diagnostic; the desktop's full graphical inspection
  surface is not reproduced.
- No execution was permitted, so syntax, ABI, numerical results, UI behavior, and native parity are
  not runtime-validated.

## Next phase

Trace the remaining native late filters and weighted role-identification path from
`MakeACOR`/`ConsensusOK`/`FinalTrim`, preserve each method as separately inspectable evidence, and
map native event erasure/re-scan state where it changes later event discovery. Do not replace the
current transparent vote until exact native weights are known. When execution is eventually
authorized, start with tiny fixed alignments and the cases in `validation-plan.md`; do not begin
with a production dataset.
