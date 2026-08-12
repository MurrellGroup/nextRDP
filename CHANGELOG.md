# Changelog

## 0.9.0-session-9

- Ported a bounded MaxChi confirmation kernel directly from the supplied DNA5/VB paths, without
  consulting alternate implementations: three variable-site match profiles, the source
  `NormalZ`/`ChiPVal2P` tail, `GetCriticalDiff`-shaped match screen, and fixed 70-site default.
- Reproduced `MakeBanWinP`'s distinct zero/length boundary map for native input-missing runs,
  accumulated prior-event erasures, and linear sequence ends. Missing windows are excluded from
  the first peak pass and stop the source-shaped growth path at the corresponding outer edge.
- Replaced repeated per-boundary window summation with capacity-light rolling left/right match
  totals for all three pair tracks. The strongest-peak scan is now `O(variable sites)` with a small
  constant; only the selected peak enters the bounded growth pass.
- Preserved MaxChi's three probability scopes in the result contract: raw one-degree-of-freedom
  tail, within-triplet `length / min(initial,grown window) * 3`, and the current project-wide
  correction. A source recheck hit is reported separately from profile availability.
- Added representative-triplet MaxChi evidence to each reconciled event and reran the same kernel
  for every finalized nonrepresentative distance-list candidate beside the primary-RDP recheck.
  MaxChi remains corroboration-only and does not change RDP-authored event coordinates or groups.
- Added a responsive MaxChi review card and finalized-list column, including pair, peak position,
  χ², window sizes, all probability levels, missing/edge filters, and explicit event-discovery
  boundary. Settings now distinguish full RDP discovery from available MaxChi confirmation.
- Extended JSON, project JSON, CSV, late-consensus capability metadata, TypeScript contracts, and
  the static source-contract checker for the MaxChi evidence path.
- Advanced engine/source labelling to `0.9.0-session-9` and reloadable projects to `v1alpha9`, while
  accepting `v1alpha1`–`v1alpha9`.
- Kept MaxChi smoothing, ordered multi-peak destruction/retry, exploratory event construction, and
  native golden comparison explicit as the next fidelity boundary.
- Kept the checkpoint source-only and deliberately uncompiled per project instruction.

## 0.8.0-session-8

- Revalidated the reattached RDP4/RDP5 DLL and VB source archives for safe paths and traced the
  manual's dataset → preliminary scan → hypothesis refinement → ordered review/export workflow.
- Completed the active `FinalTrim` RFF=0 list path after `OKSeq 6`: retained the final fixed-point
  nearest-nonrecombinant thresholds, added the correlation-gated first expansion and strict raw-tree
  parent-bounded second expansion, and preserved ascending candidate admission.
- Ported selected-role list pruning for all three `INList` positions, retaining swap-last removals,
  the collapsed/raw matrix branch, correlation free tickets, the source's asymmetric index spelling,
  and the inherited third-list loop index rather than normalizing it.
- Made final `OKSeq 15` membership active and completed `ConsensusOK` CScore with the exact ordering
  of the final-membership multiplier, `RCorrX`, and declared-`Long` matrix accumulator semantics.
- Ported `ConsensusOK`'s primary list rebuild, score/regional/class thresholds, raw/collapsed topology
  gates, exact position fallback, six-distance equivalence pass, low-score direct widening,
  straggler pass, and all-three-list restoration when any rebuilt role is empty.
- Ported the shared selected-role conservative cleanup after the RFF guard: native Single-precision
  movement sums and representative ranks, raw/direct outlier constraints, swap-last removals,
  bounded-list admission, and the always-true strict four-matrix inlier branch.
- Ported the primary-RDP branch of the post-group signal recheck: every finalized nonrepresentative
  candidate is rescanned against its role's two representatives using the supplied `LowP * 100000`
  threshold lift. Emitted, candidate-recombinant, event-overlapping, and ordinarily significant
  outcomes plus the best tract and probabilities remain separate in JSON, CSV, and review UI.
- Extended erased-tract breakpoint context from immediate contact to the manual's RDP-specific
  uncertainty rule and the supplied `CheckEnds` implementation. The port rebuilds the current
  erased-triplet information map, applies the distinct beginning/ending DLL ranges, source-shaped
  input `MissingData`, literal wrap comparison, and linear edge gates. Immediate adjacency,
  qualifying event IDs, nearest information-rich count, native range, and each warning reason remain
  separate in result JSON, on-demand alignment JSON, CSV, and review UI.
- Added a line-level supplied-source trace for native RDP breakpoint uncertainty and kept its manual
  review flag distinct from statistical breakpoint probability/confidence evidence.
- Ported the complete active BURT/`PolishBP` path as an isolated capacity-reusing unit: supplied
  three-symbol recoding and circular expansion, DNA5's 20-argument/21-start Viterbi training,
  pseudocount updates, forward/reverse posteriors, strict `0.995`/`0.999` bounds, signed nearest-
  interval matching, breakpoint selection, source-shaped missing-data adjustment, three-usable-site
  revert, and final gap relocation. A local 15-bit Microsoft-CRT `rand()` adapter preserves the
  supplied Windows DLL's seeded emission starts without depending on Emscripten libc.
- Made BURT evidence active before cyclic tract erasure and after representative-only corrections.
  Statistical 99%/95% ranges, HMM positions, signed containment state, applied/reverted movement,
  training diagnostics, and unavailable reasons are retained in result/project JSON, on-demand
  alignment JSON, CSV, and the responsive review UI; manual breakpoint edits remain authoritative.
- Corrected cyclic erasure to the active supplied `ModSeqNumY` path: both breakpoint coordinates are
  erased for linear and circular tracts, and the same inclusive definition now drives fragment
  re-entry, erased-tract review context, recombination-free output, and fragmented FASTA output.
  Equal endpoints follow the source's wrapped branch and therefore represent the full alignment,
  including in event-overlap grouping and circular result labelling.
- Reconstructed the current triplet's accumulated native `MissingData` union before each BURT pass.
  It combines source-shaped input missing runs with prior inclusive erased tracts that affected any
  current representative, matching `PolishBP`'s position in the native cycle without copying an
  alignment-sized mask for every sequence.
- Completed the manual's four common accepted-event recombination-free FASTA variants by adding
  current-group sequence removal and union-tract column removal beside tract masking and aligned
  mosaic fragments. All four share ordered-review/downstream-reconciliation gating; empty sequence
  or column results fail with an explicit explanation. Column removal uses a linear difference mask.
- Added the supplied default-enabled “polish breakpoints” setting across UI, C ABI, worker, project
  save/restore, and result metadata. Disabling it skips BenHMM work and explicitly preserves the RDP
  or manually edited coordinates.
- Fed finalized distance membership back into every role's manual two-of-three evidence rule and the
  automatic co-recombinant group used by cyclic tract erasure. Manual group overrides remain separate.
- Added memoized outside/inside direct, raw-patristic, and collapsed-patristic candidate-pair lookups
  so the new widening passes do not repeatedly rescan alignment columns.
- Removed heap allocation from the frequently called symmetric tract-overlap test by using fixed
  two-segment stack storage for linear/wrapped spans, while retaining inclusive and equal-endpoint
  source semantics.
- Removed temporary-vector allocation and concatenation from every primary triplet candidate pass;
  the scanner now appends into one capacity-retaining signal buffer without changing call order.
- Split the manual's enabled, masked, and disabled sequence states across the UI, worker/C ABI,
  core, results, and project replay. Masked rows now correctly re-enter secondary RDP, distance, and
  phylogenetic evidence; disabled rows remain tree context but cannot enter scans, event roles,
  co-recombinant groups, or finalized evidence lists.
- Added full, enabled-only, and masked-or-disabled-only FASTA exports around sequence curation.
  The two complementary partitions preserve original aligned rows, use the current UI state rather
  than requiring a scanner instance, work directly after alignment loading, bypass unrelated
  event-decision gating, and reject an empty selection explicitly.
- Restored the manual's large-dataset curation commands as modern controls: reapply the supplied
  closest-pair auto-mask, enable all, mask all, or disable all. Each operation updates the primary
  triplet count, keeps masked/disabled sets disjoint, and uses the same unsaved-analysis guard as an
  individual row edit; this also makes every row controllable when rendering is capped at 500 rows.
- Restored `FastRecCheckP`'s `EN != BE` detected-event gate so a full circular candidate run cannot
  create an equal-endpoint primary call; later manual/restored equal-endpoint events retain the
  active `ModSeqNumY` full-alignment semantics.
- Reused alignment-sized CheckEnds and BURT symbol/lattice/posterior buffers across events, avoiding
  repeated capacity growth in both new breakpoint paths; HMM lattice/backpointer/path buffers now
  follow the DLL's one-clear-per-event lifetime instead of being filled per Viterbi iteration. The
  compact per-triplet accumulated missing mask retains capacity across events as well.
- Expanded JSON, CSV, and review UI evidence with both expansion stages, selected-role removal,
  `OKSeq 15`, completed CScore, ConsensusOK admission stage, final membership, fallback status, and
  the primary-RDP post-group recheck.
- Advanced engine/source labelling to `0.8.0-session-8` and reloadable projects to `v1alpha8`, while
  accepting `v1alpha1`–`v1alpha8`.
- Refreshed the Pages workflow to the current Node-24-generation official actions (`checkout@v7`,
  `setup-node@v7`, `configure-pages@v6`, `upload-pages-artifact@v5.0.0`, and `deploy-pages@v5`) and
  explicitly includes the verified `.nojekyll` file in the Pages artifact. The application build
  remains on Node 20 and pinned Emscripten 5.0.1; the Emscripten setup wrapper is the current v16.
- Added a line-level supplied-source trace and updated the parity corpus for order, exact-equality,
  class-threshold, widening, and fallback cases.
- Kept the checkpoint source-only and deliberately uncompiled per project instruction.

## 0.7.0-session-7

- Added a first-party GitHub Pages Actions workflow. Default-branch pushes and manual dispatches
  install the locked npm graph under Node 20, provision Emscripten 5.0.1, check TypeScript, build
  the single-worker WASM/Vite application, upload `dist/`, and deploy through the `github-pages`
  environment with the required least-privilege token permissions.
- Added a Pages artifact verifier that requires the entry point, `.nojekyll`, Emscripten module,
  and a binary with the WebAssembly magic header; rejects development and root-relative asset URLs;
  and rejects symbolic/hard links before upload.
- Added a pre-build source-contract gate that compares every public C header function with its
  keepalive definition, the CMake export list, and worker calls, then aligns engine/package/lock
  versions and the emitted/imported project schema.
- Added `package-lock.json` and `npm ci` deployment instructions so Actions resolves the same
  dependency graph on each build. Lock metadata was generated with lifecycle scripts disabled;
  dependencies were not installed and project code was not executed in this checkpoint.
- Preserved repository-root, repository-subpath, and custom-domain deployment through the existing
  relative Vite base and document-relative worker WASM URL.
- Versioned the dynamically imported Emscripten loader and its `locateFile` requests with the same
  package checkpoint, preventing a newly hashed UI bundle from reusing stale stable-name WASM. The
  worker also rejects and destroys an engine whose exported version differs from the UI package.
- Kept Pages on the single-worker compatibility module because Pages cannot provide the COOP/COEP
  headers needed for pthread WASM. The numerical core remains off the DOM thread.
- Added explicit unsaved-project state after scans and every review/repair mutation. The UI now
  shows dirty/current/saving checkpoint state, warns before tab exit, and confirms destructive
  dataset, mask/settings, or re-scan replacement until a project checkpoint is downloaded.
- Advanced the engine/source checkpoint to `0.7.0-session-7` without changing project schema
  `v1alpha7` or the documented native-consensus boundary.
- Kept the checkpoint source-only and deliberately uncompiled per project instruction.

## 0.6.0-session-6

- Added a bounded C++/C-ABI/worker endpoint for inspecting the original alignment around both
  event breakpoints without transferring complete alignment rows to the main thread.
- Added an on-demand graphical breakpoint inspector to ordered review. Current recombinant and
  parent rows are always first, followed by current/automatic co-recombinant groups, masked traces,
  and distance/tree supporting evidence; rows are deduplicated and explicitly capped.
- Added 1-based coordinate rulers, circular-origin wrapping, linear-end clipping, selectable
  ±15/30/60-site context, role/group labels, and parent-aware nucleotide colouring.
- Added persistent per-event warnings when a beginning or ending boundary touches a prior erased
  tract, plus closest expected major→minor/minor→major informative-state brackets in the inspector
  and corresponding JSON/CSV fields.
- Corrected the active opening `FinalTrim` duplicate pass: after a candidate/pair is counted in more
  than one RList, the diagnostic correlation copy is cleared for every occurrence of that candidate
  in every role list, including occurrences not independently above the duplicate threshold.
- Ported the source branches that populate `FinalTrim` `OKSeq` elements 7, 8, 9, 12, and 13:
  collapsed/raw whole-event tree position, whole-tract relative JC distance, and both paired
  breakpoint JC distance scores. Native strict/tied rewards, negative penalties, warning/saturation
  gates, closest-pair role modifiers, and deliberately asymmetric cases are preserved.
- Added each mapped matrix score, availability/fallback state, and the raw active matrix subtotal to
  JSON, CSV, and a dedicated review table. These diagnostics do not change membership.
- Ported `CalcMatchY` `OKSeq` 17/18 from the supplied VB/DLL sources: exact four-flank reconstruction,
  VB banker's rounding, signed variable-site map, circular rolling window, regional product, six
  breakpoint samples, standard grouping thresholds, and source bounds. Scratch vectors are reused.
- Ported the opening `ConsensusOK` raw-tree topology-consistency rejection for `OKSeq` 18. JSON,
  CSV, and review retain both the raw and filtered class, checkpoint values, breakpoint-presence
  flags, and labelled JC fallback; this evidence remains deliberately non-pruning.
- Corrected every role-indexed comparison to the supplied `CompMat` ordering, including the role-1
  parent order, while retaining the active bare-sequence-index quirk in the positive `OKSeq 14`
  branch.
- Confirmed from the active VB path that `OKSeq` 10 is explicitly cleared and 11 is never assigned;
  both are now represented as intentional zero slots rather than missing functionality.
- Ported active `FindActualEvents` and `MakeMatchMatX2P` behavior for `OKSeq 14`: best direct event
  overlap, inverse-parent rejection, circular/equal-endpoint tract traversal, candidate/reference
  tract intersection, JC saturation, exact rewards/penalties, and the source index quirk. Prefix
  counts replace repeated full-alignment rescans without changing interval results.
- Ported `CheckPatternX`'s three-region informative-state shares and the `ConsensusOK` base inputs
  `OKSeq` 0–5 plus `RCorrX`, including representative sentinel/zero behavior.
- Ported `FinalTrim` `OKSeq 6` nearest-nonrecombinant fixed-point membership with native post-
  `StripDupInv` swap order, correlation gates, unfound-event removals, raw/collapsed tree limits,
  and paired breakpoint-distance veto.
- Added a labelled non-pruning CScore preview through `OKSeq` 0–6. It preserves the VB `NS As Long`
  half-to-even narrowing after each matrix contribution and leaves only pending final membership
  `OKSeq 15` neutral. JSON, CSV, and the review table expose every input and the raw versus source-
  rounded matrix values.
- Added machine-readable `lateConsensus` and `breakpointInspection` capability metadata to result
  JSON so consumers can distinguish complete early evidence stages from the non-pruning late stack.
- Preserved the exact edge lists, branch lengths, bootstrap support, and working-fragment leaf
  provenance from all six event NJ trees, then exposed them through a bounded on-demand C ABI and
  worker request without enlarging ordinary result or project JSON.
- Added a responsive graphical tree reviewer with whole-tract/5′/3′ paired comparisons, arbitrary-
  root disclosure, current role/group colouring, retained-fragment labels, bootstrap percentages,
  unavailable-region fallback messaging, and a reversible below-50% branch-collapse view.
- Added machine-readable `treeInspection` capability metadata alongside the breakpoint inspector.
- Advanced reloadable projects to `v1alpha7`, accepting `v1alpha1`–`v1alpha7`, and advanced the
  engine/source checkpoint to `0.6.0-session-6`.
- Kept the checkpoint source-only and deliberately uncompiled per project instruction.

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
