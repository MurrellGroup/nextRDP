# Changelog

## 0.5.0-session-5

- Replaced one-shot ordered support rebuilding with strongest-first cyclic discovery: after each
  event, the current co-recombinant tract is erased and every eligible distinct-origin triplet is
  screened again from the beginning.
- Added source-shaped fragment re-entry below the supplied 100,000-site cutoff. Working fragments
  are gap padded, retain original/event provenance, cannot form same-origin triplets, omit short or
  duplicate copies, and are bounded by an explicit 256-fragment browser cap.
- Made the Bonferroni opportunity count round specific, excluding same-origin synthetic
  combinations, and retained each signal's correction factor for faithful project replay.
- Replaced linear round-signal duplicate scans with a canonical original-triplet hash index and
  removed repeated copies of the immutable `O(N²)` identity matrix during history rebuilds.
- Added fragment-assisted event analysis: saved anchors recover the appropriate working fragment
  for distance profiles, role metrics, and the bounded event-tree panel while all reported sequence
  identities remain original alignment records.
- Ported the active `MakeACOR` topology-affinity inequalities, the exact `MakeRList` two-correlation
  override quirk, the distance-triangle warning branch, and active `StripDupInv` removal of
  inverse-only candidates. All pass/block reasons remain visible in JSON and review rows.
- Expanded role identification to nine displayed metrics. Eight use the mapped native full/half
  contribution weights, including raw/collapsed-tree PhylPro, leave-one-role-out, displacement, and
  weighted triplet-ordering paths; `nativeWeightParity` remains false while method families are missing.
- Closed the rejected-event workflow hole. Rejecting an event now invalidates the downstream chain,
  restores its tract during replay, preserves the rejected call for audit, and rediscovers later
  events in fresh cyclic passes. The UI gates decisions to analysis order.
- Added manual complete co-recombinant group editing with a preserved automatic two-of-three
  baseline. Current group membership now drives cyclic erasure and accepted-event FASTA variants.
- Enforced event order and rebuild readiness in the C++/WASM boundary, including rejection of final
  FASTA export while reviews or downstream reconciliation remain incomplete.
- Made pending project snapshots resume from a valid compact prefix: stale later signals/events are
  omitted on import, retained signal IDs and event anchors are remapped, and the rebuild marker remains.
- Added a one-click project checkpoint directly in ordered review so the manual's save-often loop
  remains practical during accepted, unreviewed, and pending-rebuild states.
- Ported the active opening `FinalTrim` duplicate-positive-correlation cleanup as non-pruning,
  per-pair JSON/UI/CSV diagnostics.
- Advanced reloadable project JSON to `v1alpha6`, accepting `v1alpha1`–`v1alpha6`, and expanded
  JSON/CSV with round correction factors, fragment provenance, erasure counts, late-correlation
  gates, and weighted role contributions.
- Kept the checkpoint source-only and deliberately uncompiled per project instruction.

## 0.4.0-session-4

- Added a standalone deterministic phylogenetics core with pairwise Jukes–Cantor distances,
  neighbour joining, canonical split matching, ten-replicate column bootstrap, and 50% branch collapse.
- Added the six event subalignments and all three paired-tree phylogenetic membership checks, with a
  100-sequence tree panel and explicitly labelled Jukes–Cantor fallback for larger active sets.
- Added iterative detectable-set closure and complete at-least-two-of-three co-recombinant groups for
  each of the three presumed-recombinant hypotheses.
- Ported `CalCR`'s five inverse category relabellings, the dominant-category `RCorrWarn` gate, and the
  positive aggregate score/target portion of `MakeRList`.
- Added an auditable six-metric recombinant-role vote using distance `PhPr`/`SubDist`, tree
  `TreePhPr`/`TreeSubDist`, topology change, and three-set support; exact native weighting remains labelled.
- Added explicit application of role recommendations in the ordered review workflow.
- Added accepted-event tract-masked and event-ordered mosaic-fragment FASTA exports, gated in the UI
  until every event is reviewed and downstream reconciliation is clear.
- Expanded project JSON and CSV with tree regions, bootstrap summaries, all evidence sets, per-sequence
  phylogenetic diagnostics, correlation relabelling/warning data, role metrics, and final group membership.
- Advanced reloadable projects to schema `v1alpha4` while retaining `v1alpha1`–`v1alpha3` import support.
- Kept the checkpoint source-only and deliberately uncompiled per project instruction.

## 0.3.0-session-3

- Added source-shaped five-region extraction around both breakpoints, including the native
  60-informative-site boundary walk and boundary inclusion/exclusion behavior.
- Added three paired six-value distance-pattern correlations with exact six-observation Pearson
  probabilities, native-style overlap eligibility, and uncorrected `P < 0.05` membership.
- Added all three presumed-recombinant role hypotheses for every reconciled event.
- Added per-sequence correlation coefficients, p-values, site diagnostics, eligibility, and
  detectable-support flags to project JSON.
- Added the conservative intersection of detectable-signal and distance-correlation evidence as a
  two-set consensus subset; it is explicitly not labelled a complete co-recombinant group.
- Expanded the review UI and CSV export with hypothesis counts, detailed correlation evidence, and
  the pending phylogenetic boundary.
- Advanced reloadable projects to schema `v1alpha3` while retaining `v1alpha1` and `v1alpha2`
  migration support.
- Kept the checkpoint source-only and deliberately uncompiled per project instruction.

## 0.2.0-session-2

- Added strongest-first reconciliation of repeated RDP triplet hits into detectable-signal event sets.
- Added the supplied greater-than-30% symmetric tract-overlap rule and two-shared-sequence indexing.
- Added relaxed RDP trace-profile checks for masked sequences, with significant and trace-only evidence separated.
- Reworked review around event hypotheses, grouped evidence, and masked trace evidence.
- Added manual role/breakpoint correction and ordered re-identification of later events.
- Added direct project import/resume, session-1 migration, and project schema `v1alpha2`.
- Changed CSV output from primary-signal rows to reconciled-event rows.
- Kept the checkpoint source-only and deliberately uncompiled per project instruction.

## 0.1.0-session-1

- Mapped the supplied RDP5 manual workflow into a five-stage browser interface.
- Added local alignment parsing, diagnostics, and the supplied RDP5 auto-mask routine.
- Added the first C++/WASM port of the primary RDP triplet screen.
- Added worker batching, progress, cancellation, provisional-signal review, and exports.
- Added responsive UI styling and static-host-relative asset configuration.
- Documented fidelity boundaries, performance design, and the native parity plan.
- Deliberately left the checkpoint uncompiled per project instruction.
