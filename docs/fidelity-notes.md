# Primary RDP fidelity notes

## Supplied-source mapping

The scanner follows the active RDP5 path reached from the supplied VB analysis-list loop and
`AlistRDP4`/`FastRecCheckP`. No alternate implementation was used.

| Supplied routine or behavior | Port location | Treatment |
| --- | --- | --- |
| `FillFSSRDP` and compressed equality lookup | `alignment.cpp`, `RdpScanner::build_profile_on` | Direct byte-state equality preserves the three information-rich-site categories |
| `FindSubSeqPB3` / `FindSubSeqPB6` | `build_profile_on` | Excludes missing/ambiguous triplet sites and records pair support plus alignment coordinates |
| `XOHomologyP2` | `compute_rolling_counts` | Uses `2 × floor(window/2) + 1`, including circular padding |
| Average-homology ordering and cycle rotation | `ranked_pairs`, `scan_triplet` | Scans middle/low alternatives and applies the additional source rotation below 0.7 |
| `FindNextP` | `append_candidate_signals` | Candidate support must exceed both alternative tracks |
| `DefineEventP2` | `append_candidate_signals` | Maps runs to tracts, respects circular/linear termination, and retains the primary `EN != BE` acceptance gate |
| `ProbCalcP2` | `binomial_tail`, `rdp_probability` | Log-space binomial tail, opportunity factor, 169-site scaling, and `1e-300` floor |
| `MakeMCCorrection` | `refresh_active_sequences`, per-signal `correction_tests` | Recomputed for every working round over valid distinct-origin triplets, with supplied `(255^4)/2` cap |
| `AutoMaskmnu_Click` | `make_suggested_mask` | Iterative closest-pair pruning and diversity-preserving representative choice |
| `EnableAllMnu_Click` / `MaskAllMnu_Click` / `DisableAllMnu_Click` | dataset bulk curation controls | Replaces every row state in one operation and recomputes the primary opportunity count |
| `MaskSeq` states 0/1/2 | sequence-state UI, `ScanOptions::mask`/`disabled`, active/trace/late/tree paths | Enabled enters every screen; masked skips primary only; disabled skips primary and secondary evidence; all remain bounded tree context |
| `SaveAlign` modes 0/5/6 | loaded-alignment curation FASTA routes | Full, masked-or-disabled, and enabled row subsets are available before scanning and preserve aligned source rows |

## Cyclic discovery and fragment state

The supplied workflow does not retain every first-pass hit as an independent final event. The port
now performs complete strongest-first cycles:

1. Screen every eligible working triplet and retain significant signals for that round.
2. Select the lowest corrected probability, then local probability, then stable signal order.
3. Group its detectable support using two shared original identities and strict `> 0.3` symmetric
   tract overlap.
4. Build all three role hypotheses and select the current complete co-recombinant group.
5. Erase the event tract from that group and begin a fresh complete triplet pass.
6. Stop when no unexplained signal remains, no new working sites can be erased, or fewer than three
   distinct active origins remain.

Below the supplied `< 100000` alignment-length condition, each erased tract also becomes a
gap-padded working fragment. `working_origins_` and `working_fragment_events_` preserve identity and
event provenance. Two working copies of one original can never occupy the same triplet. Copies with
fewer than `max(5, window, ceil(length/100))` usable states and exact same-origin duplicates are
discarded. A visible 256-fragment cap bounds browser memory and combinatorial growth; reaching it is
reported rather than hidden.

Fragment-assisted signals are mapped back to original sequence IDs for event reporting. During
downstream evidence calculation, a saved fragment context selects the matching working copy with
the most usable event-tract sites. Other fragment aliases are not all placed in each tree panel;
that is a declared performance adaptation.

## Secondary event reconciliation

| Supplied routine or rule | Port location | Treatment |
| --- | --- | --- |
| Strongest remaining probability | `finish_detection_round` | Corrected p, local p, then stable signal ID choose the next unexplained anchor |
| Two shared triplet sequences | `triplet_pair_keys`, round pair index | Each signal is indexed by its three unordered original-sequence pairs |
| `FillSetsP*` overlap | `tract_overlap` | Symmetric `2 × intersection / (size1 + size2)`, strict `> 0.3`, including origin-spanning tracts |
| `FindSets` role closure | `refresh_role_hypotheses` | Builds three direct detectable sets and repeatedly adds a sequence to the remaining role when the other two contain it |
| Masked follow-up | `refresh_trace_evidence` | Rechecks each masked sequence; trace-only and corrected-significant calls remain distinct |
| `MakeBPosLR`, `VSN = 60` | `correlation_boundary`, `build_correlation_regions` | Four circular walks with source boundary inclusion/exclusion |
| `MakeSDMP2` | `correlation_region_profile` | Parent-one, parent-two, and parent-agree/candidate-differs fractions over five regions |
| `FillRmat` | `correlation_vector` | Two breakpoint pairs plus outside-average/tract six-value pair |
| `CalCR` Pearson and relabelling | `pearson_six_two_sided_p`, `best_category_correlation` | Direct correlation, three pair swaps, and both cyclic relabellings; cycles share native class four |
| `RCorrWarn` | `correlation_warnings`, distance-triangle guard | Dominant-category warnings plus the source distance-triangle branch and XOR tract/outside warning behavior |
| `MakeGoodC` | `breakpoint_overlap_sites` | Strictly more than ten non-gap candidate sites across either breakpoint span |
| `MakeINList` / `MakeACOR` | `source_in_list`, `acceptable_affinity` | Closest outside/inside representative pair maps to the active topology-affinity inequalities; raw patristic distance is preferred with JC fallback |
| `MakeRList` | distance evidence finalisation | Preserves the exact first-two-correlation `corc == 2` override quirk and warning-adjusted positive P-score aggregate |
| `StripDupInv` | `stripped_inverse_only` final membership | Inverse `r > 0.83` remains diagnostic, while active inverse-only rows are removed from the distance set |
| Opening `FinalTrim` duplicate pass | `duplicate_filtered`, `duplicate_cleaned_support` | Warning/inversion-cleared direct pairs above `0.83` are counted across all three role lists; once a sequence/pair count exceeds one, that pair is suppressed for every occurrence in every role list, matching the source's pair-wide second pass |
| `FinalTrim` `OKSeq 6` | nearest-membership fixed-point block | Reconstructs ascending `MakeRList` order plus `StripDupInv` swap-last deletion, then repeats raw/collapsed nearest-nonrecombinant limits, strict `0.83`/`0.95`/`0.99` correlation tickets, unfound-event deletion, and paired four-breakpoint JC veto until the source last-index stop condition holds |
| `FinalTrim` `OKSeq` 7/8 | `FinalTrimMatrixEvidence` tree-position scores | Applies the source's four anchor/parent and candidate/parent comparisons to collapsed and raw outside/inside patristic matrices, including strict/tied rewards and role modifiers |
| `FinalTrim` `OKSeq` 9/12/13 | `FinalTrimMatrixEvidence` distance scores | Ports whole-tract and both breakpoint-pair JC branches, native warning/`<3` gates, positive downweighting, weaker-positive suppression, and negative quarter/half/full penalties; scores feed active CScore |
| `FindActualEvents` / `MakeMatchMatX2P`, `OKSeq 14` | `SourceDetectedTractGrid`, match-matrix prefix path | Uses the active direct-event catalogue, pre-inversion-strip membership, inverse-parent gate, strict `>1/3` overlap and best tract, circular/intersected JC distance with native saturation, all signed branches, and the source's bare-`CompMat` index lookup |
| `CheckPatternX`, `OKSeq 3` | `source_pattern_scores` | Counts source `1`/`0.5` informative-state patterns across both breakpoint spans and the tract, then emits each role's share with its representative sentinel |
| `FinalTrim` ordered expansions / `OKSeq 15` | active final-list block | Retains the fixed-point nearest-nonrecombinant thresholds, appends correlation-gated candidates in ascending order, applies strict raw-tree parent bounds in the second pass, then preserves selected-role swap-last pruning and its inherited third-list loop index before recording `OKSeq 15` |
| `ConsensusOK` completed score | `ConsensusScoreEvidence` | Maps corrected correlation P, event overlap, set/pattern/pre-strip/duplicate/nearest/final membership, `RCorrX`, representative sentinel/zero cells, and the active matrix multiplier into the score used for membership |
| `ConsensusOK` `NS As Long` | `source_vb_long`, completed score | Applies half-to-even narrowing after the averaged 7/8 term and every nonzero 9–14 addition; the native three-assignment sub-one branch narrows its final `2^-1` value to zero. The UI separately shows the readable raw subtotal |
| `CalcMatchY` `OKSeq` 17/18 | `source_calc_match`, `CalcMatchEvidence` | Preserves source flank reconstruction/bounds, VB half-to-even window selection, signed match map, circular single-precision smoothing, regional product, six checkpoint samples, standard breakpoint thresholds, and the opening `ConsensusOK` raw-tree topology filter; available values are active grouping inputs |
| `ConsensusOK` list rebuild | primary/equivalent/straggler membership block | Clears all three lists, applies score/regional/class thresholds and collapsed/raw topology fallbacks, snapshots each widening pass, preserves exact six-distance equality and direct/raw straggler comparisons, and restores all old lists if any rebuilt role is empty |
| Shared selected-role conservative cleanup | post-rebuild tree cleanup block | Accumulates direct whole-region movement sums as Single, ranks all three representatives, applies native raw/direct outlier and clade-equality constraints with swap-last removal, then adds bounded or strict four-matrix inliers; preserves the `FAMatSmall <= HDF` and always-true `Or x = x` branches |
| Six event subalignments | `build_phylogenetic_regions` | Both breakpoint flanks plus whole outside/inside tract partitions |
| Jukes–Cantor matrices | `phylogeny.cpp` | Pairwise valid sites; `-0.75 ln(1-4p/3)`; saturated/insufficient pairs use native-style distance 10 |
| Event tree path | `build_tree_region_evidence` | Deterministic NJ, ten column-bootstrap replicates, canonical split support, and collapse below 50% |
| Manual tree displays | `event_trees_json`, `EventTreeInspector` | Reuses exact saved edges for all six reconciliation trees; compares each region pair, exposes branch support/collapse and fragment provenance, and arbitrarily roots only the SVG representation |
| Large-set tree bound | `select_tree_sequences` | All three anchors plus the closest active originals, capped at 100; omitted candidates use labelled JC affinity |
| `CheckBSTree` grouping shape | phylogenetic evidence loop | Candidate-anchor distance must beat both anchor-parent and candidate-parent distances in both trees of a pair |
| Manual “two of three” rule | `complete_two_of_three_set` | Detectable, distance, and phylogenetic membership are counted independently; two memberships include the sequence |
| Ordered correction/rejection | `update_event`, `set_event_review_state`, `reconcile_after` | Accepted corrections are re-erased; rejected tracts are restored; fixed earlier events remain while all later events are rediscovered |
| Manual group repair | `update_event_group` | Keeps automatic and current groups separate; forces the recombinant, excludes current parents, and invalidates the downstream chain |
| Mid-repair save/reload | worker project restore | Drops stale later events/signals, compacts anchor IDs, restores the valid prefix, and preserves the first-changed-event marker |
| Manual sequence display | `event_alignment_json`, `EventAlignmentInspector` | On-demand immutable-alignment windows around both breakpoints prioritize current roles and relevant group/evidence rows, wrap circular coordinates, bracket expected major/minor informative-state transitions, and colour parent/recombinant patterns without moving full sequences to the UI |
| Native RDP `CheckEnds` uncertainty | `breakpoint_uncertainty`, `uncertain_erasure_event_ids`, event/alignment JSON, CSV, and UI warning | Reconstructs current triplet information positions after prior erasures, the source's ten-character input `MissingData` baseline, distinct beginning/ending DLL ranges, strict linear edge gates, and literal wrap comparison; erased-event attribution and immediate contact remain separate |
| BURT / `PolishBP` statistical confidence | `burt_confidence.cpp`, settings, result/alignment JSON, CSV, and review card | Full supplied three-state Viterbi/posterior/range/matching/repositioning chain is active with a local seeded Microsoft-CRT random adapter; the supplied flag defaults on, disabled runs preserve coordinates, and parent-transition brackets remain separate review-only evidence |
| MaxChi `FastRecCheckMC2` strongest-peak recheck | `maxchi.cpp`, event/list result evidence, CSV, and review cards | Three variable-site pair profiles, native missing/erasure and linear-edge bans, strict critical-difference gate, rolling χ² maximum, `ChiPVal2P`, source-shaped growth, and raw/within/project probability scopes are active as corroboration only; MaxChi does not discover or reposition events |
| Manual alignment outcomes | FASTA export methods | Accepted group tracts are gapped; fragment export appends aligned event fragments in event order |

Project JSON identifies this as `reconciliationTier: "detectable-distance-phylogenetic"` and
retains individual sets, complete groups, round-specific correction factors, fragment provenance,
tree-panel/fallback status, bootstrap summaries, correlation gates, and role metrics.

## Role recommendation boundary

The role panel exposes nine inspectable metrics:

| Metric | Native contribution represented here |
| --- | --- |
| `PhPr` | JC profile correlation, weight 8 |
| `TreePhPr` | Raw-tree profile correlation, weight 18 |
| `CollapsedTreePhPr` | Collapsed-tree correlation with the native `SubPhPr` condition, weight 20 |
| `SubPhPr` | Leave-one-role-out JC correlation, weight 10 |
| `TreeSubPhPr` | Leave-one-role-out raw-tree correlation, weight 8 |
| `SubDist` | JC outside/inside displacement, weight 2 |
| `TreeSubDist` | Raw-tree outside/inside displacement, weight 10 |
| `TrpScore` | Weighted triplet-ordering change, weight 8 |
| `ThreeSetSupport` | Display-only group-size context, weight 0 |

The mapped voting metrics use the supplied full contribution for a best role and half contribution
for a role better than one alternative, with the collapsed-tree special condition handled
separately. Confidence is the top-minus-second contribution margin divided by all contributions.
The full desktop method battery is not present, so JSON deliberately keeps
`nativeWeightParity: false`.

## Representation and performance choices

- Internal vectors are zero based; emitted alignment coordinates remain one based.
- A/C/G/T/U are unambiguous byte states. Gaps, ambiguity, and missing markers are excluded.
- Probability tables become direct log-space evaluation.
- Per-triplet buffers, breakpoint layouts, overlap counts, and reference distances are reused.
- Scans yield after 512 working combinations; same-origin combinations skip the numerical path.
- Sequence curation preserves the manual's three states: only enabled rows enter primary triplets,
  masked rows can enter secondary evidence, disabled rows cannot enter event roles/groups, and all
  three states remain candidates for the bounded phylogenetic panel.
- Fragment retention is bounded at 256 and event tree construction at 100 sequences.
- Tree distance matrices load a column's states once before updating all pairs.
- Six tree families are built once per event and reused across all role hypotheses and metrics.
- Graphical breakpoint context is generated only when opened and is capped at 64 rows, 100 sites
  on either side, and two windows; no full sequence row crosses the worker boundary.
- `CalcMatchY` reconstructs no more than three alignment lengths, retains no more than 160 variable
  sites, and reuses its match/smoothing scratch vectors for every ordinary candidate.
- `MakeMatchMatX2P` builds three valid/difference prefix pairs in `O(NL)` and answers all role/tract
  intersections without the source's repeated whole-alignment rescans; `CheckPatternX` is `O(NL)`.
- `OKSeq 6` reuses saved tree/JC matrices. Its order-sensitive fixed point can be `O(N²)` distance
  checks in the one-removal-per-pass case but performs no new alignment scan or tree construction.
- MaxChi keeps rolling left/right totals for all three pair tracks, so the strongest-peak pass is
  `O(V)` rather than `O(V × H)`; only the selected peak is grown and no full χ² profile is retained.

## Known parity risks requiring golden comparison

1. Informative-site upper-bound conventions can create one-site boundary or probability deltas.
2. Average-homology tie perturbations use stable pair identities rather than every native lookup table.
3. Circular tracts at coordinate one require focused `DefineEventP2` comparison.
4. The erasure helper now follows active `ModSeqNumY`'s inclusive linear/wrapped loops, including
   the source's full-alignment traversal when both endpoints are equal; exact behavior still
   requires a native golden fixture.
5. Text readers normalize metadata and do not cover every desktop format.
6. Dot-gap and missing-data arrays are separate natively; compact state zero treats them uniformly.
7. Exact FAMat/SMat categorical normalization around `MakeACOR` can change marginal affinity calls.
8. Full contradictory-inversion handling beyond active `StripDupInv` remains incomplete.
9. The active `FinalTrim`/`ConsensusOK` grouping path through shared selected-tree cleanup and the
   primary-RDP post-group recheck are ported but unvalidated. Exact floating-point equality branches,
   literal sequence-zero fallback, collapsed-matrix availability boundary, Single movement ranks,
   list ordering, empty-role restoration, and the browser mapping of the native widened `LowestProb`
   threshold need focused native fixtures.
10. The 100-sequence tree panel, 256-fragment cap, and omission of non-anchor fragment aliases from
    event trees are explicit browser adaptations.
11. NJ tie/negative-branch behavior and split support require native tree fixtures.
12. Several native recombinant-identification method families remain absent despite the mapped
    weights for the implemented subset.
13. The manual's RDP-specific uncertainty flag, supplied `CheckEnds` range construction, active
    inclusive `ModSeqNumY` erasure, and distinct BURT/`PolishBP` path are implemented but unvalidated.
    Fragment-to-original role mapping, the local Microsoft-CRT random adapter, odd-length circular
    expansion, signed ranges, accumulated missing-data repositioning, and final active source quirks
    require native fixtures. Other breakpoint probability-distribution families remain outside this
    primary-RDP slice.
14. Project replay (including pending-prefix compaction), rejected-event tract restoration, manual
    group repair, final-export gating, and cyclic termination need native fixtures.
15. Enabled/masked/disabled replay and the exact influence of disabled tree-context rows on native
    late distance pruning need focused fixtures, especially when the 100-row browser tree cap binds.
16. MaxChi confirmation is source-shaped but unvalidated. Native lookup-table rounding, zero/length
    `BanWin` aliases, `MaxX = 0`, missing-boundary growth ordering, the browser's nonoverlapping-
    half-window growth bound, and correction factors need golden comparison. Smoothing,
    `DestroyPeakP`, retry/event construction, and `AlistMC3` discovery remain unported.

Until those cases are validated, the project provides a complete review/export path for primary
RDP use cases plus MaxChi triplet corroboration, but does not claim native result parity.
