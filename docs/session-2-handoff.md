# Session 2 handoff

This checkpoint advances the port through the detectable-signal portion of RDP5's secondary
analysis phase. It remains source-only and has not been compiled or executed.

## Completed in this checkpoint

1. Primary signals are ordered by corrected p-value, local p-value, and stable signal ID.
2. The strongest unexplained signal anchors an event candidate.
3. Signals sharing two sequence IDs are found through an unordered pair index.
4. Candidate tracts must pass the supplied strict symmetric overlap rule:
   `2 × intersection / (first tract length + second tract length) > 0.3`.
5. Compatible signals are assigned to the strongest matching anchor; a signal cannot bridge two
   events merely by sharing a different pair with another supporting signal.
6. Every masked sequence is checked against the event's current parent representatives with a
   relaxed RDP profile pass; corrected significance is recorded without discarding trace-only
   profiles.
7. Event roles and breakpoints can be corrected. A correction marks later events stale; after the
   corrected event is accepted, the engine preserves earlier events and rebuilds the corrected and
   later event sets.
8. Project schema `v1alpha2` reloads the alignment without a lossy intermediate format, restores
   primary signals, replays earlier corrections, and preserves review and pending-rescan state.

## Deliberate scientific boundary

RDP5 combines the detectable-signal set with phylogenetic-correlation and distance-correlation
sets. A sequence must occur in at least two of those three sets before it is treated as a
co-recombinant. This checkpoint implements only the detectable-signal set. It therefore exports
event candidates and provisional role candidates, not final co-recombinant groups or an
RDP5-equivalent minimum-event history.

## Recommended next development phase

Port the evidence used to complete the secondary phase before adding more UI breadth:

1. Jukes-Cantor distance lists inside and outside each event tract.
2. Pearson distance-correlation significance for each presumed recombinant in the anchor triplet.
3. The six breakpoint-flanking sub-alignments and bootstrapped neighbor-joining tree movement
   evidence used for the phylogenetic-correlation set.
4. The two-of-three evidence-set consensus that forms co-recombinant groups.
5. The supplied recombinant-identification score consensus, with each component exposed in the
   review UI rather than collapsed into an unexplained label.

Only after those steps and native parity validation should recombination-free alignment exports be
unlocked.
