# RDP5 workflow mapping

The manual treats recombination analysis as hypothesis construction and ordered review, not a
single “run” operation. RDP Web uses the same staged shape.

| Manual workflow | Browser implementation | Session 5 boundary |
| --- | --- | --- |
| Load an aligned nucleotide dataset | Drag/drop local alignment; validate format and length | GenBank/ORFMap, PDB, and `.rdp5` import remain pending |
| Inspect diversity and alignment quality | Identity, missing-data, variable-site, and informative-site summaries | Full alignment editor/checker remains pending |
| Mask overly similar sequences | Supplied RDP5 iterative auto-mask; every sequence remains manually toggleable | Disabled and masked states are not yet separate |
| Select topology, methods, window, p cutoff, and correction | Circular RDP defaults, window 30, p = 0.05, Bonferroni | Only RDP is active |
| Primary scan across eligible triplets | Incremental worker scan with progress and cancellation | Primary RDP only; other methods are locked |
| Select the strongest event and rescan | Reconcile one strongest unexplained event, erase its inferred group tract, re-enter bounded fragments, then run a fresh complete pass | Fragment/tree caps are explicit browser adaptations |
| Reconcile signals into an event | Two shared original sequences plus strict >30% symmetric tract overlap | Remaining late native pruning can still alter marginal membership |
| Test all three recombinant roles | Each anchor rotates through recombinant position | Current roles are preserved until the weighted recommendation is applied |
| Build detectable sets | Direct supporting signals plus iterative cross-role closure | Native event-network parity remains a golden-test item |
| Build distance-correlation sets | Five regions; six-value Pearson; five relabellings; warnings; overlap; `MakeACOR`; `MakeRList`; `StripDupInv`; opening `FinalTrim` duplicate cleanup | Remaining `ConsensusOK`/`FinalTrim` refinement does not yet prune the group |
| Build phylogenetic-correlation sets | Six JC matrices; ten bootstrap NJ replicates; 50% branch collapse; three tree pairs | Above 100 panel sequences, labelled JC fallback controls cost |
| Form co-recombinant groups | Any original sequence present in at least two of the three evidence sets | Implemented, not native validated |
| Identify recombinant and parents | Nine visible metrics; eight use mapped native full/half weights | Missing native method families keep full-weight parity false |
| Examine masked relatives | Recheck every masked sequence against event parent representatives | Full native trace score remains pending |
| Review events in order | Profile, weighted role evidence, three sets, tree diagnostics, traces, accept/reject | Numeric rather than graphical tree/alignment inspection |
| Correct an early event | Edit/apply recommendation, record decision, re-erase corrected group, rediscover later calls | Breakpoint uncertainty beside deleted tracts remains pending |
| Correct over/under-grouping | Search and edit current co-recombinant membership while retaining the automatic two-of-three baseline | Manual group is authoritative for downstream erasure/export; native late consensus can still alter the automatic baseline |
| Reject an event | Preserve rejection for audit, restore its tract, suppress that fixed call, rediscover later calls | Runtime/native replay parity remains unvalidated |
| Save often | Review-header checkpoint plus export-page project download; `v1alpha6` retains alignment, round factors, fragment provenance, groups, edits, and decisions; a pending repair reloads only its valid prefix | Alpha schema, not an `.rdp5` replacement |
| Export tables and alignment variants | CSV, tract-masked FASTA, and event-ordered mosaic-fragment FASTA | Variants use accepted current roles/groups and remain unvalidated |

## Defaults carried over

- Circular sequences: enabled.
- Highest acceptable p-value: `0.05`.
- Multiple comparisons: Bonferroni over valid active distinct-origin triplets, recomputed each
  cyclic round and capped at the supplied native limit.
- RDP window: `30`; the source computes `2 × floor(30 / 2) + 1 = 31` information-rich sites.
- Event-specific tree bootstrap: ten replicates with branches below 50% support collapsed.
- Synthetic fragment condition: alignment length below 100,000 sites; this port additionally
  exposes a 256-fragment retention cap.
- Method panel: RDP, GENECONV, and MAXCHI occupy the manual’s primary positions; RDP alone is enabled.

## Review contract

The scanner does not hide how the selected roles or group were formed. The review screen shows the
current roles, weighted recommendation, every underlying role metric/contribution, all three sets,
the complete two-of-three group, fragment-assisted status, and whether phylogenetic evidence came
from the tree panel or bounded fallback.

Undecided events are handled in analysis order. Applying a recommendation, editing breakpoints, or
correcting group membership marks that event for downstream reconciliation. Rejecting a call does
the same because later calls were discovered after its tract had been erased. Both the interface
and WASM API reject final FASTA generation until every event is decided and no downstream chain is
stale.
