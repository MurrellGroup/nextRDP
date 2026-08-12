# Changelog

## 0.18.0-session-18

- Re-skinned the complete application as a Windows 95 workstation without changing analytical
  behavior: teal desktop, navy active title bar, classic menu and status strips, gray system face,
  square raised/recessed bevels, Tahoma/MS Sans Serif fallbacks, inset tables and plots, and
  segmented navy progress blocks.
- Applied the same visual system to workflow navigation, upload/curation controls, settings,
  notices, scan phases, review lists, breakpoint/tree/alignment inspectors, evidence cards, and
  export panels, including a compact raised drawer on narrow screens. No external visual asset or
  network dependency was added.
- Fixed stop handling during cyclic discovery. At the next 512-triplet boundary, the engine now
  discards only transient signals from the unfinished round, retains every event committed by
  prior complete rounds, clears incomplete shortlist summaries, records `user-stopped`, and runs
  normal final reconciliation so the UI moves to Review with usable results.
- Added an explicit stopped-cleanly notice and changed the control label to “Stop and review
  completed events.” Stopping the opening round safely returns a zero-event result when no event
  has yet reached a complete round boundary.
- Added a production-WASM graceful-stop regression that enters a second cyclic round, evaluates a
  partial batch, stops, reconciles, and verifies the exact completed-event prefix with no orphaned
  signals. The existing FASTA-upload and exact cyclic-shortlist digest regressions remain active.
- Advanced the engine/package version to `0.18.0-session-18`; project schema remains
  `org.rdp-web.project/v1alpha16` because the saved analysis shape did not change.

## 0.17.0-session-17

- Traced the supplied GENECONV probability lifecycle through `GCCalcPValP2`, `GCXoverD`,
  `MakeMCCorrection`, and later BURT polishing. WebRDP's raw KA probability is calculated from the
  detected fragment before BURT; RDP5's stored `XOverList.Probability` corresponds to WebRDP's
  project-corrected value. Review UI now states both scopes and the exact initial opportunity
  multiplier explicitly.
- Fixed cyclic multiple-comparison fidelity by freezing the source `MakeMCCorrection` factor from
  the initial exploratory or query/reference scan plan. Later fragment-expanded schedules still
  report their real triplet workload, but no longer silently change stored probability scope.
- Added source-shaped cyclic shortlist reuse based on `XOverList`/`XOverDefine`, `BestXOList`,
  `Worthwhilescan`, and `StoreLPV`: significant summaries for unchanged exact working triplets are
  replayed, empty unchanged summaries skip their method kernels, and only triplets touching erased
  rows or new fragments run fresh full discovery.
- Preserved 3SEQ's semantic boundary by refreshing it once on every unchanged triplet when the scan
  first enters the supplied post-erasure `CheckSplit3Seq` mode; later rounds can safely reuse it.
  Manual event repair/rejection rebuilds clear the shortlist and start with a full pass.
- Added live counters for triplet kernel evaluations, summaries/signals reused, and method scans
  skipped. A deterministic two-event production-WASM Pages regression asserts fixed correction,
  cached signal replay, skipped work, and fewer kernel evaluations than scheduled triplets.
- Compared the optimized and temporary full-rescan builds on an all-five-family cyclic fixture;
  event roles, breakpoints, rounds, support order, signal triplets, and p-values produced the same
  digest across 22 events/422 signals while kernel evaluations fell from 131,488 to 32,911.
- Recorded two deliberate GENECONV safety fixes: a bounded root-bracketing fallback is used only
  when the supplied Newton iteration becomes unstable (and is counted), and the apparent
  `NDiff[3] == 1` typo is treated as the intended assignment to avoid degenerate division by zero.
- Completed WASM and web production builds plus FASTA-upload, cyclic-shortlist, source-contract,
  TypeScript, and GitHub Pages artifact verification.

## 0.16.1-session-16

- Fixed FASTA and project uploads under Emscripten 5.0.1 by exporting the `HEAPU8` view used by
  the worker to copy file bytes into WASM memory.
- Added a GitHub Pages artifact smoke test that instantiates the production WASM module, copies a
  three-sequence FASTA through that heap view, loads it, and validates the alignment summary.
- Completed a clean Node 20 / Emscripten 5.0.1 production build and Pages artifact verification.

## 0.15.0-session-15

- Ported the supplied late `TSXOver(1)` 3SEQ Findall shape for the representative event triplet and
  every finalized nonrepresentative distance-list row. It scans all three source target rotations,
  preserves the initial lower-P/low-information/corrected gate, and then evaluates both walk
  orientations rather than allowing the ordinary strictly-better-orientation shortcut.
- Retained forced post-erasure `CheckSplit3Seq`/`SubPVal` evaluation for both Findall orientations,
  the distinct post-`SwapRound` correction branch, and the extra inverse-parent/inverse-interval
  `XOverList` copy for each accepted orientation. Recheck evidence records profile/evaluation counts,
  qualifying orientations/list entries, best target/direction, bounds, excursions, probability
  route, split state, raw/corrected P, and cutoff result without moving event coordinates.
- Added event and post-group 3SEQ recheck JSON, review cards/badges/table cells, detailed CSV columns,
  and late-consensus diagnostics. The late evidence status now names five active families while
  retaining `nativeThreeSeqFullRecheckComplete=false` until golden and full catalogue parity exist.
- Ported the supplied later-round 3SEQ `FindSubSeqTS2` → `CheckSplit3Seq`/`SubPVal` path. The
  compact `upper_bound` mapping preserves `XPosDiff[x]` as the retained-site count at or before an
  alignment coordinate without allocating a length-sized integer map per target.
- Preserved the two-stage `TSXOver` gate: only an initially significant call enters split handling;
  beginning- and ending-side missing runs are tested, the reverse orientation is retried only when
  its unsplit probability beats the trimmed selection, strict lower-P swapping is retained, and
  project correction is recalculated before emission.
- Reused the fused triplet missing/erasure byte map and compact walk, so later-round split handling
  adds bounded coordinate scans and probability evaluations but no alignment-byte pass. Split state,
  revised bounds, raw/corrected probability, route, and excursion survive JSON, CSV, review, and
  hardened restore.
- Changed the compact exact DP state from extended precision to supplied `Single` precision. This
  narrows every accumulation like the desktop lookup and materially lowers bounded WASM memory and
  arithmetic cost.
- Corrected the scaled `GetTSPVal` fallback to use the supplied post-truncation `onM/nM` or
  `onN/nN` exponent and the `1e-300` positive-underflow floor. Also separated the source's pre-wrap
  probability excursion from the post-`CheckwrapC` boundary excursion throughout the contract.
- Aligned the visible and serialized kernel order with the actual source call sequence:
  `FindSubSeqTS` → `Seq3PVals`/`GetTSPVal` → `CheckwrapC` → `TSXOver`.
- Advanced the engine to `0.15.0-session-15` and project schema to `v1alpha15`. Imports accept
  `v1alpha1`–`v1alpha15`; v14 retains ordinary 3SEQ support but cannot claim v15 split evidence.
- Kept manual 3SEQ permutation envelopes and full late event-catalogue reconstruction explicit
  rather than claiming the new non-coordinate-changing corroboration is complete dispatcher parity.

## 0.14.0-session-14

- Ported ordinary automated 3SEQ directly from the supplied manual and active `FindSubSeqTS` /
  `FindSubSeqTS2`, `Seq3PVals`/`Get3SeqPvalC`, `GetTSPVal`, `CheckwrapC`,
  `SiegmundDiscrete`, `SwapRound`, and `TSXOver` routes, without consulting alternate
  implementations.
- Added the three source-order candidate-recombinant rotations and corrected target-to-parent
  equality-slot mappings `{0,2,1}` / `{1,0,2}`. Parent-equal, all-different, and missing sites are
  excluded, while the DLL's effective four-information-rich-site minimum is retained.
- Preserved maximum descent/ascent discovery, the source's probability-before-`CheckwrapC` call
  order, origin extension and beginning advance, linear interval conversion, strict lower-P
  orientation choice, provisional parent swap, and both supplied low-information exits. JSON,
  restore, CSV, and review retain pre-wrap probability and post-wrap boundary excursions separately.
- Replaced the desktop four-dimensional exact lookup with a compact finite hypergeometric
  maximum-drawdown dynamic program. Results narrow to source `Single`, cache by `(m,n,k)`, and obey
  an eight-million-transition / 8192-entry browser bound.
- Added the supplied large-profile `SiegmundDiscrete`/`ApproxNu` route and bounded scaled-exact
  fallback, including its post-truncation exponent ratio and positive-underflow floor, with exact
  and approximate evaluations retained separately. Literal zero-P emission,
  the `p > 10^-15` Dunn–Šidák branch, and the smaller-tail product route are preserved.
- Reused the existing fused RDP/MaxChi/CHIMAERA/GENECONV equality preparation so 3SEQ adds no full
  alignment-byte pass. Method code 4 joins the source-major RDP → GENECONV → MaxChi → CHIMAERA →
  3SEQ strongest-first cyclic scheduler and the existing grouping, erasure, fragment, consensus,
  BURT, ordered review/repair, and export workflow.
- Added enable/default state, correction explanation, target/exact/approximate/candidate progress,
  complete JSON/restore/CSV anchor evidence, a signed three-target random-walk plot, an anchor
  evidence card, and visible approximation disclosure.
- Added `rdp_restore_threeseq_discovery`, four authoritative restore counters, CMake source/export
  coverage, and hardened import checks for counts, excursion bounds, role rotations, pair IDs,
  low-information calls, post-erasure split claims, and unknown method labels.
- Advanced the engine to `0.14.0-session-14` and project schema to `v1alpha14`. Imports accept
  `v1alpha1`–`v1alpha14`; all pre-v14 projects restore 3SEQ disabled so saved discovery semantics
  cannot silently change.
- Added the supplied-source 3SEQ trace, Session 14 handoff, status/workflow/architecture/validation
  updates, and static contracts for pair mappings, source gates, bounded exact/fallback work, ABI,
  restore, UI, plots, exports, version, schema, and Pages deployment.
- Kept the ambiguous `FindSubSeqTS2` alias, `CheckSplit3Seq`/`SubPVal` post-erasure re-probability,
  manual permutation envelopes, late-list corroboration, authorized native golden comparison, and
  runtime validation explicit as unresolved boundaries. No compiler, bundler, preview server,
  dependency install, or project/runtime test was invoked.

## 0.13.0-session-13

- Ported ordinary automated triplet GENECONV directly from the supplied manual and active
  `FindSubSeqGCAP6`, `GetFragsP`, `GetMaxFragScoreP`, `CalcKMaxP`,
  `GCCalcPValP2`, lowest-P/overlap helpers, and `Module31.bas::GCXoverD`, without consulting
  alternate implementations.
- Added the three inner pair-match and three outer discordant-sequence signed tracks with the
  supplied six-way provisional recombinant/minor/major role mapping, ignored-indel default, skew
  rejection, inactive automated minimum-fragment predicates, and configurable one-overlap default.
- Reused the existing fused RDP/MaxChi/CHIMAERA variable-site/equality preparation, so enabling
  GENECONV does not add another complete triplet alignment-byte scan.
- Preserved source mismatch penalties, Newton initialization, lambda/K calculation, strict
  critical-score gate, single-precision narrowing points, and Karlin–Altschul low/high-score
  probability branches. Added a bounded bisection fallback only when the native unbounded Newton
  iteration becomes numerically unsafe, with an exported diagnostic counter.
- Restored literal critical/tail rounding and native zero-underflow handling, preventing an exact
  zero KA probability from being promoted to an artificial positive display-floor value.
- Replaced the direct quadratic positive-fragment extension with signed prefix scores,
  next-strictly-lower bounds, and rightmost range maxima. This preserves the native first-negative
  stop and latest-ending tie while reducing each track to `O(R log R)`.
- Corrected circular construction to preserve `GetFragsP`'s initial category run plus terminal
  origin-wrapped copy without repeating the complete signed list. One category pass now builds all
  six tracks and coalesces only the positive outer fragments coalesced by the supplied helper.
- Preserved stable track/fragment lowest-raw-P ordering and implemented the overlap allowance with
  a lazy range-add/range-maximum tree in `O(F log V)`.
- Joined independently labelled GENECONV signals to the strongest-first cyclic scheduler, shared
  cross-method grouping, three-set reconciliation, BURT polishing, event correction/rebuild,
  fragment re-entry, and final alignment exports. Exact corrected/local-P ties use the supplied
  RDP → GENECONV → MaxChi → CHIMAERA method-major order.
- Added non-coordinate-changing ordinary-kernel GENECONV corroboration for event representatives
  and every finalized distance-list row. It reuses MaxChi's prepared triplet workspace and retains
  best-track/tract/role, raw/corrected probability, skew, overlap, workload, and fallback evidence
  in project JSON, the review matrix, and CSV.
- Added GENECONV settings, validation, live counters, method-aware scan state, complete per-signal
  score/KA evidence, three-colour inner/outer `-log10(raw P)` review plot, anchor card, event
  method badges, CSV diagnostics, and interpretation boundaries.
- Added ABI/worker restore for method code 3, `rdp_restore_geneconv_discovery`, three scan
  settings, five authoritative progress counters, CMake source/export coverage, and cross-layer
  parameter-arity checks.
- Advanced the engine to `0.13.0-session-13` and project schema to `v1alpha13`. Imports accept
  `v1alpha1`–`v1alpha13`; pre-v13 projects restore GENECONV disabled so saved discovery
  semantics cannot silently change.
- Added the supplied-source GENECONV trace, Session 13 handoff, architecture/workflow/fidelity/
  validation documentation, and static contracts across kernel, fused preparation, ABI, restore,
  types, settings, progress, review, plots, exports, version, schema, and Pages deployment.
- Kept permutation probabilities, manual pair scans, alternative indel modes, full late native
  event reconstruction, and all native-versus-WASM golden validation explicit as unresolved boundaries. No
  compiler, bundler, preview server, dependency install, or project/runtime test was invoked.

## 0.12.0-session-12

- Ported ordinary automated CHIMAERA discovery directly from the supplied manual §8.5.1,
  `FindSubSeqDP3/6`, `FastRecCheckChim`, `AlistChi`, and active VB `CXoverA` dispatch without
  consulting alternate implementations.
- Added all three candidate-recombinant target rotations and their target-specific information-rich
  binary strings: monomorphic and all-different sites are excluded, target/parent-one matches encode
  `1`, and target/parent-two matches encode `0`.
- Reused the combined RDP/MaxChi alignment-byte pass for CHIMAERA input. Target profiles are filtered
  from cached equality tracks, avoiding three additional complete triplet-row reads while preserving
  independently capacity-retained workspaces.
- Added fixed/fallback CHIMAERA window selection with supplied default 60, missing/earlier-erasure and
  linear-end bans, rolling raw χ² construction, source-shaped critical screen, and a per-target lazy
  peak heap with strict raw-score/coordinate ordering.
- Added source-shaped window growth, far-side tract choice, optimized second breakpoint, alignment-
  coordinate mapping, completed/rejected smoothed-basin destruction, three-consecutive-wasted exit,
  and the supplied 100-attempt bound independently for every target rotation.
- Preserved raw, within-triplet position/rotation-corrected, and optional project-corrected
  probabilities separately, plus initial/grown windows, critical difference, flank χ², filter flags,
  and inside/outside parent-one match rates.
- Joined method-labelled CHIMAERA signals to the strongest-first RDP/MaxChi cyclic scheduler. The
  selected target and provisional parent order remain auditable before shared late role consensus,
  BURT confidence, tract erasure, fragment re-entry, and fresh complete passes.
- Added CHIMAERA enable/window settings, target-profile/peak/candidate/limit progress counters,
  method-aware null and retry-limit notices, a dedicated evidence card, and a single-trace target/
  parent-one χ² plot rather than misrepresenting the profile as MaxChi's three pair traces.
- Expanded result/project JSON and event CSV with CHIMAERA methods and complete anchor diagnostics.
  Added ABI/worker restore for method code 2, all discovery evidence, four authoritative workload
  counters, and the exported `rdp_restore_chimaera_discovery` symbol.
- Added secondary `FastRecCheckChim` representative/finalized-list corroboration across all three
  target rotations. It reuses MaxChi's prepared triplet bytes, retains raw/within/project-corrected
  strongest-target evidence in project JSON/review/CSV, and never moves an event's coordinates.
- Persisted a review and CSV warning when MaxChi and CHIMAERA are the only supporting methods, since
  the manual treats those methods as closely related rather than independent confirmation.
- Advanced the engine to `0.12.0-session-12` and project schema to `v1alpha12`. Imports accept
  `v1alpha1`–`v1alpha12`; v10/v11 projects retain MaxChi but restore CHIMAERA disabled, while projects
  through v9 retain their RDP-only discovery semantics.
- Added filter-wide reference assignment, first-appearance group compaction, exact dataset
  eligibility feedback, per-signal query/reference role metadata, live cyclic working-role counts,
  and honest zero-work terminal progress retained from the opening audit of this session.
- Added a supplied-source CHIMAERA trace, Session 12 handoff, validation cases, and static contracts
  across the kernel, fused preparation, C ABI/CMake/worker bridge, restore gates, types, settings,
  progress, review, plots, exports, version, and schema.
- Kept CHIMAERA permutations, full late-list event reconstruction, native lookup rounding, and all
  native-versus-WASM golden comparison explicit as unresolved boundaries. No compiler, bundler,
  preview server, dependency install, or project/runtime test was invoked.

## 0.11.0-session-11

- Ported the supplied automated query-vs-reference analysis-list workflow from manual §4.2 and
  `Module3.bas::MakeAnalysisListQvR`, without consulting alternate implementations. Each scheduled
  primary triplet contains one query and two references from different reference groups; later role
  inference remains unconstrained so a reference can still be identified as recombinant.
- Added dataset-level query/reference assignment with positive numeric group IDs, one-click parsing
  of documented `REF-A<name>`-style names, an all-queries reset, per-row status, filter-aware
  multi-row assignment beyond the 500-row rendering cap, first-appearance group compaction, and
  validation that at least one eligible query plus references from two groups exist before the scan
  plan proceeds.
- Added fully exploratory/query-vs-reference scheme selection, exact initial workload counts,
  scheme-aware scan/review/export language, live round-specific correction counts, result summaries,
  CSV role/assignment columns, and full reference-group persistence in project checkpoints.
- Made cyclic progress report the current active working-row, query-record, reference-record, and
  represented-group counts after every erasure/re-entry round; legitimate terminal zero-work rounds
  are no longer visually replaced by the initial workload.
- Added the manual's visible reference-as-recombinant distinction: event navigation and detail use
  an amber flag, all three current role cards/edit choices expose their query/reference input role,
  event JSON stores the selected recombinant's input role/group, and export summaries count these
  calls without forcing a query-based interpretation.
- Replaced the desktop routine's materialized triplet array with an `O(N)` active query/reference
  catalogue and a lazy cross-group reference-pair/query cursor. Source pair/query order is retained;
  bounded worker batches, cancellation, cyclic tract erasure, and fragment re-entry use the same
  scheduler without allocating `O(T)` triplet records.
- Kept the exact scheduled reference-record workload distinct from the supplied
  `choose(reference groups, 2) × query origins` multiple-testing rule, retained the native correction
  cap, added saturating 64-bit arithmetic, and exposed the distinction in settings and saved results.
- Made query/reference identity follow original-sequence provenance through working fragments,
  excluded same-origin working copies, recomputed groups/correction/workload after every erasure or
  repair, and added a dedicated no-eligible-constrained-triplets terminal state.
- Fixed schedule rebuilding for corrections at event zero and all-rejected saved prefixes by
  unconditionally refreshing active sequence, role, cursor, total, and correction state after replay.
- Hardened v11 restore by rejecting repeated or constraint-breaking saved triplets, rebuilding a
  missing grouped correction from the scheme rather than raw workload, filling legacy per-signal
  factors consistently, and validating every C ABI role-buffer length instead of silently defaulting.
- Completed method-correct MaxChi plot reconstruction by feeding the selected triplet's supplied
  ten-character input-missing map into banned-window generation, while retaining forced breakpoint,
  peak, and pair-maximum samples and the later-round reconstruction disclosure.
- Advanced the engine to `0.11.0-session-11` and projects to `v1alpha11`. Imports accept
  `v1alpha1`–`v1alpha11`; older projects remain exploratory, and projects through `v1alpha9` retain
  MaxChi discovery disabled so saved detection semantics do not change.
- Extended the source-only contract checker across the constrained scheduler, C ABI/worker pointer
  bridge, restore compatibility boundary, options/results/types, dataset/settings/scan/review/export
  UI, version/schema markers, and delimiter balance. No compiler, bundler, preview, or project
  runtime was invoked.

## 0.10.0-session-10

- Ported the supplied ordinary-triplet MaxChi `MCXoverF` discovery path as an independently
  labelled method stream beside RDP, without consulting alternate implementations.
- Added all three variable-site match profiles, fixed/fallback window selection, native
  critical-difference screening, `MissingData`/prior-erasure and linear-edge exclusions, and the
  source's strict raw χ² peak ordering across boundary first and pair second.
- Added source-shaped peak growth, `FindSide` left/right tract selection, `OptLeftBPMC` and
  `OptRightBPMC` breakpoint refinement, missing-boundary adjustment, alignment-coordinate mapping,
  and preliminary roles that remain subject to the shared late role consensus.
- Reproduced the literal `SmoothChiValsP` eleven-divisor/twelve-term off-by-one behavior and padded
  terminal cell used by the supplied destruction routines.
- Ported accepted-region `DestroyPeaks` and rejected-peak `DestroyPeakP` behavior with bounded circular
  traversal, the source's numerical peak-inclusion adjustment, accepted-hit `WasteOfTime` reset,
  three-consecutive-failure exit, and 100-attempt cap.
- Replaced repeated full-profile maximum scans with a linear-time heap build and lazy invalidation of
  destroyed raw peaks. Rolling left/right match totals still make all three initial χ² profiles
  linear in variable sites rather than window width. The combined path prepares RDP categories,
  MaxChi matches/variable prefixes, and the triplet missing/erasure map in one alignment-byte pass.
- Joined RDP and enabled MaxChi signals in the strongest-first cyclic scheduler. Either method can
  anchor an event; cross-method supporting signals share the existing two-identities/strict-30%-
  overlap grouping, late role/group evidence, tract erasure, bounded fragment re-entry, and fresh
  complete-pass lifecycle.
- Kept `FastRecCheckMC2` representative/finalized-list MaxChi confirmation separate from exploratory
  discovery, so its evidence cannot silently replace detected or manually edited event coordinates.
- Added method-aware progress counters, signal/event badges, RDP pair-identity versus MaxChi χ²
  plots, and a responsive MaxChi anchor card exposing attempt, side, pair, peak, windows, flank χ²,
  raw/within/project probabilities, and missing/linear filters. Scan and review now warn when any
  triplet exhausts the supplied 100-peak retry bound with positive raw peaks still present.
- Moved the native global-maximum significance gate ahead of smoothing and retry accounting, so
  rejected triplets avoid unnecessary `SmoothChiValsP` work and report zero attempted peaks just as
  the supplied `MCXoverF` path does.
- Made plot thinning retain both inferred breakpoints, the selected MaxChi peak, and each pair
  profile's maximum. Later-round/fragment-assisted plots are now visibly labelled as original-
  alignment reconstructions instead of being mistaken for the historical erased working profile.
- Expanded CSV with anchor/detection methods and full MaxChi discovery diagnostics; expanded result
  and project JSON with discovery method, method set, settings, counters, and the per-signal trace.
- Added C ABI and worker restore support for method-labelled signals and MaxChi discovery state.
  Completed-project reload also retains cumulative rounds and all four MaxChi workload/limit
  counters, the exact terminal reason, and final-round processed/total triplets without inflating
  cumulative work while saved events are replayed. Advanced projects to `v1alpha10`; imports accept
  `v1alpha1`–`v1alpha10`, and older
  projects restore MaxChi discovery disabled so their saved RDP-only detection semantics do not
  change.
- Extended the source-only contract checker across the MaxChi kernel, C ABI/CMake/worker bridge,
  schema/version, settings, result types, review UI, and legacy-project compatibility boundary.
- Added a supplied-source MaxChi discovery trace, a Session 10 handoff, and a golden corpus for
  smoothing/index aliases, peak ties, growth/side optimization, destruction basins, retry bounds,
  method-aware plots/exports, and combined cyclic ordering.
- Kept native lookup-table rounding, preliminary role selection, permutation/manual-doublet modes,
  and all native-versus-WASM golden comparison explicit as unresolved parity boundaries.
- Kept the checkpoint source-only and deliberately uncompiled per project instruction.

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
