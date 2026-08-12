# Browser architecture and performance

## Runtime boundary

The React interface owns workflow state and presentation. A module worker owns the WASM instance,
alignment, cyclic working alignment, scan cursor, event list, and exports. The main thread transfers
the input `ArrayBuffer` rather than copying it and receives compact JSON summaries.

The worker advances 512 working combinations per batch and yields to its event loop between
batches. This keeps progress and cancellation deterministic without placing numerical work on the
DOM thread. At the end of a pass, the C++ core selects/reconciles one strongest event, calculates its
secondary evidence, erases its tract, and reports a new-pass status. The worker yields once before
starting that pass. Event reconciliation/tree construction is worker-isolated but atomic, so
cancellation is available between scan batches rather than inside one event’s finishing stage.

## State model

The immutable alignment retains normalized text, compact states, diagnostics, and the original
pair-identity matrix. The cyclic alignment contains only names, sequences, compact states, basic
sequence summaries, and length; it deliberately does not copy the `O(N²)` identity matrix whenever
a correction or rejection rebuilds history.

Each working record has:

- an original sequence ID;
- a fragment-source event ID, or `-1` for an original record;
- an aligned compact-state row in which non-fragment coordinates are missing.

Primary-active rows must be enabled (neither masked nor disabled) and have at least
`max(5, window)` usable states. Masked originals are retained for secondary RDP evidence and every
original—including disabled rows—remains eligible for bounded tree context. A triplet is eligible
only when its three original IDs are distinct. Erasure below 100,000 alignment sites can append
fragments; short/duplicate copies are omitted and at most 256 fragments are retained.

## Complexity and memory

- Alignment encoding: `O(NL)` time, one compact byte per sequence/site plus normalized source text.
- Original pair identities: `O(N²L)` once, stored as an `N × N` float matrix.
- Auto-masking: initial `O(N²)` nearest-neighbor pass; affected rows are refreshed as representatives
  are masked.
- Cyclic RDP scans: `O(Σ T_r L)` where `T_r` is the valid distinct-origin triplet count in round
  `r`. Same-origin synthetic combinations are skipped before profile construction and excluded from
  the correction factor.
- Probability: at most 169 log-space binomial recurrence steps after long-tract scaling.
- Round signal deduplication: expected constant-time hashed bucket lookup over canonical original
  triplet, recombinant, and breakpoints, followed by exact collision comparison.
- Detectable-signal grouping: pair-index construction is `O(S_r)` for retained signals in a round;
  lookup uses three unordered original-sequence pairs plus tract overlap.
- Erasure/re-entry: `O(EWL)` across event history in the conservative case, with working copies
  bounded by `W ≤ N + 256`.
- Masked trace follow-up: `O(EML)` for `E` events and `M` masked, non-disabled sequences with one
  parent pair per event role currently rechecked. Disabled rows incur tree-matrix context cost only.
- Role/correlation evidence: `O(ENL)` after event grouping. Five coordinate lists are shared; each
  candidate is reduced to fixed category counts and three six-value correlations.
- Event phylogenetics: six regions and eleven matrices per region (base plus ten bootstraps). NJ is
  `O(K³)` for `K = min(panel candidates, 100)`. Column-state loading makes each JC matrix
  `O(K²L_region)`; omitted original candidates use cached anchor-relative JC affinity.
- Role consensus: nine fixed metrics over three roles after distance/tree matrices exist.
- Mapped late-matrix scoring: constant score work per retained candidate after the same six
  distance/tree matrix families exist. Detected-region matching builds valid/difference prefixes in
  `O(NL)` and queries tract intersections in constant time; pattern shares add one `O(NL)` pass.
- `FinalTrim` nearest-nonrecombinant membership: saved raw/collapsed tree and breakpoint JC lookups
  only. The native order-sensitive fixed point is worst-case `O(N²)` if one list entry is removed per
  iteration; no alignment is rescanned and no tree is rebuilt.
- `ConsensusOK` rebuilding: primary and two widening passes are worst-case `O(N²)`. Symmetric
  candidate-pair direct/raw/collapsed distances are memoized per outside/inside region, so a fallback
  pair scans alignment columns at most once per region instead of once per branch comparison.
- Selected-tree cleanup: the native all-pairs direct movement ranks are `O(N²)` and reuse the same
  caches; list constraints and strict-inlier admission add `O(N)` work per role afterward.
- Primary-RDP post-group rechecks: `O(GL)` across the three finalized role-list memberships `G`.
  The information-rich profile scratch buffers retain their capacity across candidates, and only
  compact counts plus the best overlapping tract are kept in ordinary results.
- MaxChi confirmation: `O(GL)` across the same finalized memberships plus one representative
  triplet per event. Each triplet makes one `O(L)` variable-site/missing-boundary pass and scans all
  three pair tracks with rolling left/right totals in `O(V)`, where `V ≤ L`; only the strongest
  peak enters bounded window growth. No per-breakpoint half-window rescan or full χ² profile is
  retained.
- Native `CheckEnds` uncertainty: immutable input `MissingData` is `O(NL)` once. Per-event erased-
  triplet map/range reconstruction is conservatively `O(EL)` and therefore `O(E²L)` over history;
  only reason flags, native ranges, prior event IDs, and nearest information-rich counts are retained.
- `CalcMatchY` scoring: four source-bounded flank walks and at most 160 retained variable sites;
  rolling smoothing is `O(NV)` for `V ≤ 160`, with two reusable scratch vectors per role.
- Alignment exports: tract mask/fragments are `O(ELG + NL)`, where `G` is the accepted current
  co-recombinant group size; sequence removal is `O(EG + N)`, and a difference mask makes accepted-
  tract column removal `O(E + NL)` rather than rescanning all event tracts for every row. Full,
  enabled-only, and masked/disabled-only row exports are single `O(NL)` passes over the loaded
  alignment, use curation masks supplied with the request, and need neither a scanner instance nor
  an event-readiness gate.
- Plot output: capped near 2,000 profile points per selected signal.
- Breakpoint inspection: requested on demand as two original-alignment windows only, bounded to
  64 prioritized rows and 100 flanking sites per side; deleted-tract uncertainty is precomputed as
  compact event IDs/counts, and the complete alignment stays in the worker.
- BURT/BenHMM confidence uses fixed three-state contiguous matrices plus capacity-retaining symbol,
  lattice, posterior, path, and coordinate buffers. Each event preserves the supplied 21 seeded
  starts and at most 100 Viterbi updates per start; only compact interval and movement evidence is
  retained after the worker completes the linear passes. Lattice, backpointer, and path storage is
  cleared once per event, matching the DLL allocation lifetime rather than per training iteration;
  disabling the supplied setting bypasses the HMM altogether. A capacity-retaining `O(L)` triplet
  mask merges input missing runs with relevant prior inclusive erasures before repositioning.
- Tree inspection: six event topologies are already constructed for reconciliation; only their
  `O(K)` edge lists and common leaf metadata are retained and transferred on demand. Ordinary
  results/projects keep bootstrap summaries rather than topology payloads or `O(K²)` matrices.

## Implemented hot-path choices

- Compact integer nucleotide states rather than string comparisons.
- Constant-time rolling-window updates plus reused profile/count and primary-scan candidate buffers;
  candidate passes append directly instead of allocating and concatenating temporary vectors.
- Round-local hashed signal deduplication, including reordered fragment aliases.
- Allocation-free fixed two-segment tract overlap for signal grouping and post-group rechecks.
- Log-space binomial evaluation and significance rejection before UI retention.
- Round-specific correction opportunities computed from distinct original identities.
- Unordered pair indexing for event support rather than all-signal pairwise comparison.
- Shared breakpoint layouts, overlap counts, representative profiles, and reference distances.
- Prefix-counted detected-region comparisons instead of the DLL's repeated whole-alignment
  `MakeMatchMatX2P` scans.
- Event-scoped symmetric distance caches for the active `ConsensusOK` candidate-pair passes.
- Reused information-rich profile buffers for masked traces and finalized-list RDP rechecks.
- Rolling MaxChi totals across all three variable-site pair tracks, with compact boundary-ban and
  match arrays; strongest-peak-only confirmation avoids the pending smoothing/destroy/retry profile.
- Six tree families built once per event and reused by all role hypotheses and metrics.
- Canonical split keys for bootstrap support rather than tree-string comparison.
- A lightweight mutable alignment reset that omits the original pair-identity matrix and parser
  diagnostics.
- Dedicated worker isolation and bounded scheduling.
- Event-scoped alignment slices rather than main-thread transfer of complete sequence rows.
- Event-scoped topology edge lists rather than main-thread transfer of tree distance matrices or
  review-time tree reconstruction.

## Project restoration

Project import does not round-trip the saved alignment through FASTA or another lossy wrapper. The
worker transfers each saved name/normalized sequence through a record-oriented C ABI, restores
signals with their own correction factors and fragment provenance, then replays events in order.
Every tract/fragment state needed by a later saved anchor is rebuilt before that anchor’s evidence.
Manual role/breakpoint/group edits, rejected calls, decisions, and any pending downstream
invalidation marker are retained. The breakpoint-polishing flag is also restored, with older
schemas taking the supplied enabled default. If a snapshot was saved during a pending repair, only events
through the changed call and their supporting signals are restored. Signal IDs are compacted,
event anchors are remapped, and the stale tail is deliberately left for the next cyclic scan to
rediscover. This prevents pre-correction evidence from being replayed as current state.

Ordered review is independently enforced by the C++ boundary. Only the first undecided event (or
the event holding a pending rebuild marker) can be changed; revisiting an earlier accepted event is
still allowed and creates a new downstream invalidation point. The same boundary rejects final
alignment variants unless the scan is complete, every event is decided, and no rebuild is pending.

The emitted schema is `org.rdp-web.project/v1alpha9`; the import path accepts `v1alpha1` through
`v1alpha9` and supplies conservative defaults for fields absent from an older checkpoint.

## Checkpoint lifecycle

The UI treats a completed scan as unsaved until project JSON is downloaded. Accept/reject actions,
role or breakpoint edits, co-group edits, and downstream re-identification mark it unsaved again.
Importing a project or successfully initiating a project download marks the displayed state
current. While state is unsaved—or a scan/reconciliation is active—the browser's `beforeunload`
guard protects tab close/reload. Replacing the dataset, mask/settings, or an existing result with a
new scan also requires explicit confirmation when the latest state has not been checkpointed.

This is deliberately a loss-prevention layer rather than implicit persistence: no alignment or
analysis is written to browser storage, and no sequence data leaves the tab.

## Static deployment

The checked-in Pages workflow runs only for a manual dispatch or a push to the repository's actual
default branch. It restores the locked npm graph under Node 20, provisions Emscripten 5.0.1, checks
the strict TypeScript contract, builds the single-worker C++/WASM target and Vite application, then
uploads only `dist/` through the official Pages artifact/deployment path.
The wrappers use the current Node-24-generation action releases; hidden-file inclusion is explicit
so the verified `.nojekyll` marker is not dropped by the artifact action.

`vite.config.ts` uses `base: "./"`; the worker resolves `wasm/` from `document.baseURI`. Therefore
the same artifact works at an account root, repository subpath, or custom domain without knowing
the repository name at build time. A post-build verifier rejects development-entry remnants,
root-relative asset URLs, missing/empty loader output, invalid WASM magic, and file links that are
not accepted in a Pages artifact. The `.nojekyll` marker is included defensively.

Vite content-hashes the UI bundles, while Emscripten emits stable loader/WASM filenames. The
client therefore sends the package version to the worker and applies the same version query to
both the dynamic `.mjs` import and every `locateFile` request. A deployment cannot pair the new UI
with a cached engine binary from an older checkpoint. After instantiation, the worker also compares
the native `_rdp_version()` result with that package version and destroys/refuses a mismatched
context rather than analysing data across an ABI-skewed UI/engine pair.

Pages cannot attach COOP/COEP headers, so its workflow intentionally omits the optional pthread
module. The non-threaded WASM instance remains isolated in the module worker; a different static
host can separately build the optional thread target when it supplies cross-origin isolation.

## Later performance work

- Golden profiling before changing source-equivalent behavior.
- Pattern-compressed bootstrap matrices if column resampling dominates.
- Triangular original pair-identity storage for large `N`.
- Chunked/compressed project export to avoid duplicating very large normalized alignments.
- A source-supported alternative to the 256-fragment cap if real datasets reach it routinely.
- Real parallel triplet partitions in the optional pthread build; the current target is build
  plumbing only and does not divide scan work.
- Resumable in-progress round checkpoints. Completed analyses/review state are resumable now.
- If MaxChi becomes a discovery method, preserve rolling profiles while adding the supplied
  smoothing and ordered peak-destruction scheduler; do not fall back to per-window rescans.

The non-threaded worker is the compatibility target because ordinary static hosting does not always
supply the COOP/COEP headers required by `SharedArrayBuffer`.
