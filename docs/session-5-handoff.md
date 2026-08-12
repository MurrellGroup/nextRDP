# Session 5 handoff

This is a source-only checkpoint. Nothing in the supplied desktop source or browser port was
compiled, bundled, served, previewed, type-checked, or executed.

## Delivered in this phase

1. Replaced the prior single-pass event list with strongest-first cyclic discovery. A complete
   eligible-triplet pass selects one event, builds its evidence/group, erases that group’s tract,
   and starts a fresh pass. Unselected hits do not leak into the next modified-alignment round.
2. Added source-shaped synthetic fragment re-entry below the supplied 100,000-site cutoff. Each
   retained tract is gap padded and tagged with original sequence and source-event provenance.
   Same-origin copies cannot occupy one triplet.
3. Added practical browser bounds without concealing them: fragments shorter than
   `max(5, window, ceil(1% alignment length))` are omitted, exact same-origin copies are deduplicated,
   and at most 256 working fragments are retained. JSON and review notices report cutoff/cap state.
4. Recomputed the Bonferroni opportunity count for every round over valid distinct-origin triplets.
   Each signal now stores its own correction factor, and project replay restores that factor before
   recalculating its event evidence.
5. Added fragment-assisted evidence calculation. Event identities remain original sequence IDs,
   while correlation profiles, representative distances, event trees, and role metrics can use the
   saved working fragment that generated the anchor.
6. Ported active `MakeINList`/`MakeACOR` topology-affinity gating, including raw-patristic preference
   and bounded JC fallback; added the source distance-triangle warning and exact XOR warning behavior.
7. Preserved the active `MakeRList` `corc == 2` quirk over its first two correlations, separated
   positive and inverse support, and implemented active `StripDupInv` removal of inverse-only rows.
8. Expanded recombinant identification from the session-4 equal vote to nine displayed metrics.
   Eight carry mapped native full/half weights: PhPr, raw/collapsed-tree PhPr, two leave-one-role-out
   scores, two displacement scores, and TrpScore. Three-set size remains zero-weight context.
9. Closed downstream invalidation for rejected events. Rejecting a detected event now requires a
   rebuild that restores its tract, preserves the rejected record for audit, fixes earlier calls,
   and rediscovers all later events from a clean cyclic pass.
10. Tightened the review UI around the manual’s order: later undecided calls remain inspectable but
    cannot be decided or edited before the next event; a correction/rejection blocks downstream
    work until its rebuild finishes.
11. Added manual co-recombinant group correction. The searchable editor locks in the current
    recombinant, excludes the current parents, preserves the automatic two-of-three baseline, and
    applies the corrected group to later erasure and accepted-event alignment exports.
12. Enforced ordered review below the UI. The C++ boundary rejects out-of-order decisions/edits,
    unsolicited rebuilds, and final FASTA exports while an event is undecided or downstream state
    is stale.
13. Made mid-repair project reload resumable. A pending snapshot restores only the valid event
    prefix and its signal evidence, compacts/remaps anchor IDs, preserves the changed event and
    marker, and leaves later events to the required fresh cyclic scan. Review now exposes a direct
    project-checkpoint action so this state can be saved without leaving the event workspace.
14. Ported the active opening `FinalTrim` duplicate-correlation pass. Direct-polarity pairs repeated
    across competing role lists are marked in JSON, CSV, and the review ledger without prematurely
    applying the still-unported late pruning stack.
15. Advanced project JSON to `org.rdp-web.project/v1alpha6`. Imports accept `v1alpha1` through
    `v1alpha6`; signals carry correction/provenance fields and events carry original/working erasure,
    fragment count, fragment-assisted status, current/automatic groups, and manual group state. CSV
    contains the same high-level audit fields.

## Performance shape

- Primary and cyclic passes yield every 512 working combinations in the dedicated worker.
- Same-origin synthetic combinations are rejected before profile construction and excluded from the
  reported correction opportunity count.
- Round signals use a canonical original-triplet hash index with exact collision checks instead of
  scanning the accumulated round vector; reordered fragment aliases deduplicate correctly.
- Short-copy filtering, exact-copy deduplication, and the 256-fragment cap bound re-entry growth.
- The six event tree families are built once and reused by every role; NJ remains capped at 100
  sequences, with explicit candidate-level fallback outside the panel.
- Fragment representatives are used only when provenance requires one. Other aliases do not inflate
  the tree panel.
- Rebuilding working history copies compact sequence state but omits the immutable `O(N²)` original
  identity matrix and parser diagnostics.

## Fidelity boundary

- `MakeACOR`, the active positive `MakeRList` path, its dual-r override, `StripDupInv`, and the first
  `FinalTrim` duplicate-correlation cleanup are now represented. Exact native distance-matrix
  preprocessing can still change marginal calls.
- Complete `ConsensusOK` and the remaining `FinalTrim` duplicate/pattern/tree-boundary/distance
  pruning stack are not implemented.
- Eight mapped role metrics use native weights, but the full native method battery is missing;
  `nativeWeightParity` therefore remains false.
- The fragment cap, tree cap, and omission of non-anchor fragment aliases from event trees are
  explicit browser performance adaptations.
- Erasure endpoint behavior and breakpoint uncertainty beside a previously deleted tract require
  supplied-desktop golden fixtures.
- The alignment and tree evidence surface is diagnostic rather than a full desktop graphical editor.
- No execution was permitted, so syntax, ABI, numerical results, UI behavior, and native parity are
  not runtime validated.

## Source-only audit

- C/C++ and changed TypeScript/TSX files passed comment/string-aware delimiter scans.
- All 30 public C ABI names and parameter counts match between `rdp_api.h`, `rdp_api.cpp`, the
  worker module interface, and the Emscripten export list.
- No conflict markers, trailing whitespace, generated `.wasm`/`.mjs`/object files, dependency
  directory, build directory, or bundled site is included.
- These checks do not replace compilation, type checking, browser execution, or numerical tests.

## Next phase

Port `ConsensusOK` incrementally by first retaining every needed evidence component in project JSON,
then continue the active `FinalTrim` branches after the duplicate-correlation stage. Keep each late
decision auditable. Add breakpoint uncertainty for deleted-neighbor cases and the remaining
role-method families. When execution is explicitly authorized, begin with tiny fixed fixtures from
`validation-plan.md`, never a production alignment.
