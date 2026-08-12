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
| `FindNextP` | `candidate_signals` | Candidate support must exceed both alternative tracks |
| `DefineEventP2` | `candidate_signals` | Maps runs to tracts and respects circular/linear termination |
| `ProbCalcP2` | `binomial_tail`, `rdp_probability` | Log-space binomial tail, opportunity factor, 169-site scaling, and `1e-300` floor |
| `MakeMCCorrection` | `refresh_active_sequences`, per-signal `correction_tests` | Recomputed for every working round over valid distinct-origin triplets, with supplied `(255^4)/2` cap |
| `AutoMaskmnu_Click` | `make_suggested_mask` | Iterative closest-pair pruning and diversity-preserving representative choice |

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
| Opening `FinalTrim` duplicate pass | `duplicate_filtered`, `duplicate_cleaned_support` | Warning/inversion-cleared direct pairs above `0.83` are counted across all three role lists; repeated pairs are suppressed in the diagnostic copy exactly as the active source stage does |
| Six event subalignments | `build_phylogenetic_regions` | Both breakpoint flanks plus whole outside/inside tract partitions |
| Jukes–Cantor matrices | `phylogeny.cpp` | Pairwise valid sites; `-0.75 ln(1-4p/3)`; saturated/insufficient pairs use native-style distance 10 |
| Event tree path | `build_tree_region_evidence` | Deterministic NJ, ten column-bootstrap replicates, canonical split support, and collapse below 50% |
| Large-set tree bound | `select_tree_sequences` | All three anchors plus the closest active originals, capped at 100; omitted candidates use labelled JC affinity |
| `CheckBSTree` grouping shape | phylogenetic evidence loop | Candidate-anchor distance must beat both anchor-parent and candidate-parent distances in both trees of a pair |
| Manual “two of three” rule | `complete_two_of_three_set` | Detectable, distance, and phylogenetic membership are counted independently; two memberships include the sequence |
| Ordered correction/rejection | `update_event`, `set_event_review_state`, `reconcile_after` | Accepted corrections are re-erased; rejected tracts are restored; fixed earlier events remain while all later events are rediscovered |
| Manual group repair | `update_event_group` | Keeps automatic and current groups separate; forces the recombinant, excludes current parents, and invalidates the downstream chain |
| Mid-repair save/reload | worker project restore | Drops stale later events/signals, compacts anchor IDs, restores the valid prefix, and preserves the first-changed-event marker |
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
- Fragment retention is bounded at 256 and event tree construction at 100 sequences.
- Tree distance matrices load a column's states once before updating all pairs.
- Six tree families are built once per event and reused across all role hypotheses and metrics.

## Known parity risks requiring golden comparison

1. Informative-site upper-bound conventions can create one-site boundary or probability deltas.
2. Average-homology tie perturbations use stable pair identities rather than every native lookup table.
3. Circular tracts at coordinate one require focused `DefineEventP2` comparison.
4. The erasure helper currently preserves the prior endpoint-exclusive interpretation; supplied
   `ModSN` boundary behavior requires a golden inclusive/exclusive fixture.
5. Text readers normalize metadata and do not cover every desktop format.
6. Dot-gap and missing-data arrays are separate natively; compact state zero treats them uniformly.
7. Exact FAMat/SMat categorical normalization around `MakeACOR` can change marginal affinity calls.
8. Full contradictory-inversion handling beyond active `StripDupInv` remains incomplete.
9. `ConsensusOK` and the remaining `FinalTrim` score, pattern, tree-boundary, and
   breakpoint-distance branches after the opening duplicate cleanup are not yet ported.
10. The 100-sequence tree panel, 256-fragment cap, and omission of non-anchor fragment aliases from
    event trees are explicit browser adaptations.
11. NJ tie/negative-branch behavior and split support require native tree fixtures.
12. Several native recombinant-identification method families remain absent despite the mapped
    weights for the implemented subset.
13. Breakpoint uncertainty next to a previously erased tract is not yet reported.
14. Project replay (including pending-prefix compaction), rejected-event tract restoration, manual
    group repair, final-export gating, and cyclic termination need native fixtures.

Until those cases are validated, the project provides a complete review/export path for primary
RDP use cases but does not claim native result parity.
