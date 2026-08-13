# RDP5 workflow mapping

The manual treats recombination analysis as hypothesis construction and ordered review, not a
single “run” operation. RDP Web uses the same staged shape.

| Manual workflow | Browser implementation | Session 20 boundary |
| --- | --- | --- |
| Load an aligned nucleotide dataset | Drag/drop local alignment; validate format and length | GenBank/ORFMap, PDB, and `.rdp5` import remain pending |
| Inspect diversity and alignment quality | Identity, missing-data, variable-site, and informative-site summaries | Full alignment editor/checker remains pending |
| Mask or disable sequences | Supplied RDP5 iterative auto-mask plus explicit enabled/masked/disabled states; individual selectors and auto-mask/enable-all/mask-all/disable-all actions remain usable above the 500-row rendering cap; masked rows skip primary triplets but enter secondary checks, while disabled rows are tree-context only | Source-only and not runtime validated; the 100-row event-tree cap remains an explicit browser bound |
| Choose fully exploratory or automated query-vs-reference analysis (§4.2) | Per-row query/reference group editor, filter-aware bulk assignment beyond the render cap, first-appearance group compaction, documented `REF-A<name>` detection, and a scheme selector. Constrained triplets contain one query and two references from different groups; reference roles remain eligible to be called recombinant and are visibly flagged | `MakeAnalysisListQvR` scheduling and correction are source-shaped but not native golden validated; name detection is editable and deliberately preserves visible curation rather than silently clearing masks as the desktop prompt does |
| Select topology, methods, window, p cutoff, correction, and breakpoint polishing | Circular defaults, RDP window 30, enabled MaxChi window 70, enabled CHIMAERA window 60, enabled GENECONV G-scale 1/one overlap, opt-in primary/secondary BootScan window 200/step 20/100 replicates/70% support/seed 3, enabled no-window 3SEQ, p = 0.05, project correction, supplied BURT flag enabled; the method panel distinguishes discovery from confirmation | Primary distance-mode BootScan is active but defaults off while unvalidated; its tree/similarity/permutation/manual modes and the remaining method modes stay explicitly pending |
| Primary scan across eligible triplets | Incremental combined RDP/GENECONV/BootScan/MaxChi/CHIMAERA/3SEQ worker scan with method-aware progress and cancellation; BootScan reuses bounded shared-pair profiles, while a lazy cursor selects either every exploratory triplet or only one-query/two-cross-group-reference triplets | Source-shaped and unvalidated; primary BootScan/cache has linked host and production-WASM regression coverage, while authorized native golden comparison remains separate |
| Select the strongest event and rescan | Reconcile one strongest unexplained event, erase its inferred group tract, re-enter bounded fragments, then run a fresh complete pass | Fragment/tree caps are explicit browser adaptations |
| Reconcile signals into an event | Two shared original sequences plus strict >30% symmetric tract overlap across RDP, GENECONV, BootScan, MaxChi, CHIMAERA, and 3SEQ signals; any enabled method may anchor the strongest event and coordinates | BootScan/MaxChi/CHIMAERA/GENECONV/3SEQ preliminary roles feed the same late role arbitration; cross-method native ordering still needs golden validation |
| Test all three recombinant roles | Each anchor rotates through recombinant position | Current roles are preserved until the weighted recommendation is applied |
| Build detectable sets | Direct supporting signals plus iterative cross-role closure | Native event-network parity remains a golden-test item |
| Build distance-correlation sets | Five regions; six-value Pearson; five relabellings; warnings; overlap; `MakeACOR`; `MakeRList`; `StripDupInv`; duplicate cleanup; nearest-nonrecombinant `OKSeq 6`; both expansions; active matrix `OKSeq` 7–14; selected pruning; `OKSeq 15`; `CalcMatchY` 17/18; complete CScore; all three `ConsensusOK` passes; shared selected-tree cleanup; primary-RDP finalized-list rechecks; MaxChi/CHIMAERA strongest-peak rechecks; ordinary six-track GENECONV; and `TSXOver(1)` 3SEQ Findall for each finalized nonrepresentative | Source-only and unvalidated; full native method-stack event reconstruction and other method-family rechecks remain |
| Build phylogenetic-correlation sets | Six JC families; supplied float Clearcut NJ; Microsoft-CRT `SEQBOOT2`; retained-base pseudocount; `TreeMidP`/`UltraTreeDistP`; below-50% parent-rank promotion; raw/collapsed `MakeTreeArrayXP2` rank distances; three tree pairs | Above 100 panel sequences, labelled JC fallback controls cost; source-shaped path and explicit writer/parser repair are host-regression-tested but still await native golden matrices |
| Form co-recombinant groups | Any original sequence present in at least two of the three evidence sets | Implemented, not native validated |
| Identify recombinant and parents | Nine visible metrics; eight use mapped native full/half weights | Missing native method families keep full-weight parity false |
| Examine masked relatives | Recheck every masked sequence against event parent representatives | Full native trace score remains pending |
| Review events in order | Profile, weighted role evidence, three sets, traces, bounded original-alignment windows, and lazy paired graphical views of all six saved NJ topologies before accept/reject | Alignment/tree inspectors are source-only and not yet runtime validated |
| Correct an early event | Edit/apply recommendation, record decision, re-erase corrected group, rediscover later calls | The manual's RDP erased-tract one-window uncertainty rule and separate BURT/`PolishBP` statistical 95%/99% confidence/repositioning evidence are exposed; manual breakpoint edits stay authoritative |
| Correct over/under-grouping | Search and edit current co-recombinant membership while retaining the automatic post-ConsensusOK two-of-three baseline | Manual group remains authoritative for downstream erasure/export |
| Reject an event | Preserve rejection for audit, restore its tract, suppress that fixed call, rediscover later calls | Runtime/native replay parity remains unvalidated |
| Save often | Review-header checkpoint plus export-page project download; a dirty/current indicator, tab-exit warning, and destructive-replacement guard track whether the latest scan/review state was downloaded; `v1alpha18` retains alignment, event-tree provenance, analysis scheme/reference groups, every discovery trace, BootScan cache/work counters, fragment provenance, late rechecks, final membership stages, groups, edits, and decisions | Alpha schema, not an `.rdp5` replacement; pre-v17 imports keep primary BootScan disabled and all earlier method gates remain conservative |
| Export tables and alignment variants | CSV; full, enabled-only, and masked/disabled-only row subsets available immediately after load; plus the common accepted-event variants: remove recombinant sequences, remove recombinant columns, tract-mask each current group, or split aligned mosaic fragments | Event variants use accepted current roles/groups and all remain unvalidated |

## Defaults carried over

- Circular sequences: enabled.
- Highest acceptable p-value: `0.05`.
- Multiple comparisons: one active opportunity count over valid distinct-origin triplets,
  recomputed each cyclic round and capped at the supplied native limit. RDP/MaxChi/CHIMAERA/
  GENECONV use Bonferroni multiplication; 3SEQ uses the supplied Dunn–Šidák form.
- Automated query-vs-reference: `0`/blank is a query and each positive integer is a reference group;
  progress counts every scheduled cross-group reference-record/query triplet, while the supplied
  correction is `choose(active reference groups, 2) × active query origins` before the native cap.
- RDP window: `30`; the source computes `2 × floor(30 / 2) + 1 = 31` information-rich sites.
- MaxChi discovery: enabled, with a fixed `70`-site requested window; its source window fallback is
  applied when a triplet has fewer usable variable sites.
- CHIMAERA discovery: enabled, with the supplied fixed `60`-site requested window independently
  applied to each candidate-recombinant information-rich binary profile.
- GENECONV discovery: enabled for the supplied ordinary automated triplet route, with ignored
  indels, mismatch scale `G=1`, one overlapping fragment, strict KA critical score, and the three
  desktop minimum-fragment predicates explicitly inactive because they are commented out in
  `GCXoverD`.
- 3SEQ discovery: enabled; it has no automated window control and uses three target rotations,
  exact bounded hypergeometric random-walk tails, then the supplied `SiegmundDiscrete` fallback.
- Primary and secondary BootScan: available but disabled by default while unvalidated. The supplied
  distance-mode defaults are window `200`, step `20`, 100 counted profiles (including the
  unresampled replicate zero), 70% support, and Microsoft-CRT seed `3`.
- Polish detected breakpoints with BURT/BenHMM: enabled; disabling it preserves detected or
  manually edited coordinates and skips the HMM pass.
- Event-specific tree bootstrap: ten replicates with branches below 50% support collapsed.
- Synthetic fragment condition: alignment length below 100,000 sites; this port additionally
  exposes a 256-fragment retention cap.
- Method panel: RDP, GENECONV, BOOTSCAN, MAXCHI, CHIMAERA, and 3SEQ retain their source method-major order;
  all six ordinary triplet discovery routes are active. 3SEQ split and late Findall are active;
  its manual permutation/full event-catalogue modes and
  GENECONV permutations/manual pairs/alternative indel handling and full late event reconstruction
  remain separate pending modes. The MaxChi
  `FastRecCheckMC2` event/list confirmation path remains separate from `MCXoverF` discovery.

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
For a MaxChi-authored anchor, the screen exposes its discovery attempt number, tract side, variable-
site pair and peak, initial/grown windows, flank χ² values, raw/within/project probabilities, and
missing/linear filters. Method badges distinguish RDP and MaxChi support. The separately labelled
MaxChi representative/finalized-list recheck remains corroborative and never silently replaces the
event's detected or manually edited coordinates. Plot thinning retains the event bounds and peak;
when a cyclic signal depended on prior erasure or a retained fragment, the graph is labelled as an
original-alignment reconstruction rather than the unavailable historical working trace.
For a CHIMAERA-authored anchor, the screen separately exposes the candidate-recombinant rotation,
target/parent-one trace, information-rich count, parent-one inside/outside contrast, raw peak,
window growth, flank choice, breakpoint result, all probability scopes, filters, and retry attempt.
It never presents the one target profile as MaxChi's three pairwise profiles.
For a GENECONV-authored anchor, the screen exposes the inner/outer track and provisional role,
polymorphic/positive/discordant counts, mismatch penalty, fragment/critical scores, lambda/K, raw
and project-corrected KA probabilities, ignored-indel state, configured overlap rule, and inactive
minimum filters. Its three-colour plot pairs the supplied inner and outer tracks and does not
mislabel fragment envelopes as sliding-window or χ² profiles.
When MaxChi and CHIMAERA are the only methods supporting an event, a persisted review/export flag
states that they are closely related and must not be interpreted as independent confirmation.
For a constrained analysis, every current role, repair choice, breakpoint-alignment row, and saved
tree leaf also retains its query/reference input role. A reference selected as recombinant receives
the manual's distinct amber warning without being rejected or relabelled automatically.

Undecided events are handled in analysis order. Applying a recommendation, editing breakpoints, or
correcting group membership marks that event for downstream reconciliation. Rejecting a call does
the same because later calls were discovered after its tract had been erased. Both the interface
and WASM API reject final FASTA generation until every event is decided and no downstream chain is
stale.
