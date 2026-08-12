# RDP5 workflow mapping

The manual treats recombination analysis as hypothesis construction and ordered review, not a
single “run” operation. RDP Web uses the same staged shape.

| Manual workflow | Browser implementation | Session 9 boundary |
| --- | --- | --- |
| Load an aligned nucleotide dataset | Drag/drop local alignment; validate format and length | GenBank/ORFMap, PDB, and `.rdp5` import remain pending |
| Inspect diversity and alignment quality | Identity, missing-data, variable-site, and informative-site summaries | Full alignment editor/checker remains pending |
| Mask or disable sequences | Supplied RDP5 iterative auto-mask plus explicit enabled/masked/disabled states; individual selectors and auto-mask/enable-all/mask-all/disable-all actions remain usable above the 500-row rendering cap; masked rows skip primary triplets but enter secondary checks, while disabled rows are tree-context only | Source-only and not runtime validated; the 100-row event-tree cap remains an explicit browser bound |
| Select topology, methods, window, p cutoff, correction, and breakpoint polishing | Circular RDP defaults, window 30, p = 0.05, Bonferroni, supplied BURT flag enabled; the method panel distinguishes discovery from confirmation | RDP discovery is active; MaxChi is confirmation-only |
| Primary scan across eligible triplets | Incremental worker scan with progress and cancellation | Primary RDP only; MaxChi exploratory discovery and other methods remain locked |
| Select the strongest event and rescan | Reconcile one strongest unexplained event, erase its inferred group tract, re-enter bounded fragments, then run a fresh complete pass | Fragment/tree caps are explicit browser adaptations |
| Reconcile signals into an event | Two shared original sequences plus strict >30% symmetric tract overlap; source-shaped MaxChi strongest-peak evidence corroborates the representative triplet | RDP authors the event and coordinates; MaxChi multi-peak event construction remains pending |
| Test all three recombinant roles | Each anchor rotates through recombinant position | Current roles are preserved until the weighted recommendation is applied |
| Build detectable sets | Direct supporting signals plus iterative cross-role closure | Native event-network parity remains a golden-test item |
| Build distance-correlation sets | Five regions; six-value Pearson; five relabellings; warnings; overlap; `MakeACOR`; `MakeRList`; `StripDupInv`; duplicate cleanup; nearest-nonrecombinant `OKSeq 6`; both expansions; active matrix `OKSeq` 7–14; selected pruning; `OKSeq 15`; `CalcMatchY` 17/18; complete CScore; all three `ConsensusOK` passes; shared selected-tree cleanup; primary-RDP finalized-list rechecks; and MaxChi strongest-peak corroboration for each finalized nonrepresentative | Source-only and unvalidated; full MaxChi and other method-family rechecks remain |
| Build phylogenetic-correlation sets | Six JC matrices; ten bootstrap NJ replicates; 50% branch collapse; three tree pairs | Above 100 panel sequences, labelled JC fallback controls cost |
| Form co-recombinant groups | Any original sequence present in at least two of the three evidence sets | Implemented, not native validated |
| Identify recombinant and parents | Nine visible metrics; eight use mapped native full/half weights | Missing native method families keep full-weight parity false |
| Examine masked relatives | Recheck every masked sequence against event parent representatives | Full native trace score remains pending |
| Review events in order | Profile, weighted role evidence, three sets, traces, bounded original-alignment windows, and lazy paired graphical views of all six saved NJ topologies before accept/reject | Alignment/tree inspectors are source-only and not yet runtime validated |
| Correct an early event | Edit/apply recommendation, record decision, re-erase corrected group, rediscover later calls | The manual's RDP erased-tract one-window uncertainty rule and separate BURT/`PolishBP` statistical 95%/99% confidence/repositioning evidence are exposed; manual breakpoint edits stay authoritative |
| Correct over/under-grouping | Search and edit current co-recombinant membership while retaining the automatic post-ConsensusOK two-of-three baseline | Manual group remains authoritative for downstream erasure/export |
| Reject an event | Preserve rejection for audit, restore its tract, suppress that fixed call, rediscover later calls | Runtime/native replay parity remains unvalidated |
| Save often | Review-header checkpoint plus export-page project download; a dirty/current indicator, tab-exit warning, and destructive-replacement guard track whether the latest scan/review state was downloaded; `v1alpha9` retains alignment, round factors, fragment provenance, RDP/MaxChi rechecks, final membership stages, groups, edits, and decisions | Alpha schema, not an `.rdp5` replacement |
| Export tables and alignment variants | CSV; full, enabled-only, and masked/disabled-only row subsets available immediately after load; plus the common accepted-event variants: remove recombinant sequences, remove recombinant columns, tract-mask each current group, or split aligned mosaic fragments | Event variants use accepted current roles/groups and all remain unvalidated |

## Defaults carried over

- Circular sequences: enabled.
- Highest acceptable p-value: `0.05`.
- Multiple comparisons: Bonferroni over valid active distinct-origin triplets, recomputed each
  cyclic round and capped at the supplied native limit.
- RDP window: `30`; the source computes `2 × floor(30 / 2) + 1 = 31` information-rich sites.
- Polish detected breakpoints with BURT/BenHMM: enabled; disabling it preserves primary RDP or
  manually edited coordinates and skips the HMM pass.
- Event-specific tree bootstrap: ten replicates with branches below 50% support collapsed.
- Synthetic fragment condition: alignment length below 100,000 sites; this port additionally
  exposes a 256-fragment retention cap.
- Method panel: RDP, GENECONV, and MAXCHI occupy the manual’s primary positions; RDP performs
  discovery, while MaxChi's strongest-peak kernel is enabled only for event/list confirmation.

## Review contract

The scanner does not hide how the selected roles or group were formed. The review screen shows the
current roles, weighted recommendation, every underlying role metric/contribution, all three sets,
the complete two-of-three group, fragment-assisted status, and whether phylogenetic evidence came
from the tree panel or bounded fallback. The graphical breakpoint inspector loads only when opened,
uses the immutable original alignment, and places the current roles before co-group, trace, and
supporting-evidence rows so a browser row cap cannot hide the three representatives.
Its closest expected-state brackets remain manual review intervals rather than statistical CIs.
The separate supplied BURT/`PolishBP` HMM path now reports its signed source-labelled 99%/95%
ranges, HMM positions, coordinate movement, missing/gap adjustments, and revert state.
The graphical tree inspector is also lazy: it reuses the six edge lists built for reconciliation,
shows outside/inside pairs for the whole event and both boundaries, labels retained fragments, and
states that its internal-node root is chosen only for display.
Mapped late matrix scores are shown separately with both the readable raw subtotal and exact native
declared-`Long` `NS` multiplier. The screen exposes raw and topology-filtered breakpoint classes,
`OKSeq 15`, completed CScore, the ConsensusOK admission stage, and empty-role restoration. These
values now determine the automatic distance set.
MaxChi is displayed as a separate corroboration layer: its variable-site count, pair track, peak,
window growth, χ², raw tail, within-triplet probability, corrected probability, and ban reasons are
visible, but it neither shifts the RDP coordinates nor claims an independently discovered event.

Undecided events are handled in analysis order. Applying a recommendation, editing breakpoints, or
correcting group membership marks that event for downstream reconciliation. Rejecting a call does
the same because later calls were discovered after its tract had been erased. Both the interface
and WASM API reject final FASTA generation until every event is decided and no downstream chain is
stale.
