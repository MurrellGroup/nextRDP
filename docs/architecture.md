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

Active rows must be unmasked and have at least `max(5, window)` usable states. A triplet is eligible
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
- Masked trace follow-up: `O(EML)` for `E` events and `M` masked sequences with one parent pair per
  event role currently rechecked.
- Role/correlation evidence: `O(ENL)` after event grouping. Five coordinate lists are shared; each
  candidate is reduced to fixed category counts and three six-value correlations.
- Event phylogenetics: six regions and eleven matrices per region (base plus ten bootstraps). NJ is
  `O(K³)` for `K = min(panel candidates, 100)`. Column-state loading makes each JC matrix
  `O(K²L_region)`; omitted original candidates use cached anchor-relative JC affinity.
- Role consensus: nine fixed metrics over three roles after distance/tree matrices exist.
- Alignment exports: `O(ELG + NL)`, where `G` is the accepted co-recombinant group size.
- Plot output: capped near 2,000 profile points per selected signal.

## Implemented hot-path choices

- Compact integer nucleotide states rather than string comparisons.
- Constant-time rolling-window updates and reused profile/count buffers.
- Round-local hashed signal deduplication, including reordered fragment aliases.
- Log-space binomial evaluation and significance rejection before UI retention.
- Round-specific correction opportunities computed from distinct original identities.
- Unordered pair indexing for event support rather than all-signal pairwise comparison.
- Shared breakpoint layouts, overlap counts, representative profiles, and reference distances.
- Six tree families built once per event and reused by all role hypotheses and metrics.
- Canonical split keys for bootstrap support rather than tree-string comparison.
- A lightweight mutable alignment reset that omits the original pair-identity matrix and parser
  diagnostics.
- Dedicated worker isolation and bounded scheduling.

## Project restoration

Project import does not round-trip the saved alignment through FASTA or another lossy wrapper. The
worker transfers each saved name/normalized sequence through a record-oriented C ABI, restores
signals with their own correction factors and fragment provenance, then replays events in order.
Every tract/fragment state needed by a later saved anchor is rebuilt before that anchor’s evidence.
Manual role/breakpoint/group edits, rejected calls, decisions, and any pending downstream
invalidation marker are retained. If a snapshot was saved during a pending repair, only events
through the changed call and their supporting signals are restored. Signal IDs are compacted,
event anchors are remapped, and the stale tail is deliberately left for the next cyclic scan to
rediscover. This prevents pre-correction evidence from being replayed as current state.

Ordered review is independently enforced by the C++ boundary. Only the first undecided event (or
the event holding a pending rebuild marker) can be changed; revisiting an earlier accepted event is
still allowed and creates a new downstream invalidation point. The same boundary rejects final
alignment variants unless the scan is complete, every event is decided, and no rebuild is pending.

The emitted schema is `org.rdp-web.project/v1alpha6`; the import path accepts `v1alpha1` through
`v1alpha6` and supplies conservative defaults for fields absent from an older checkpoint.

## Later performance work

- Golden profiling before changing source-equivalent behavior.
- Pattern-compressed bootstrap matrices if column resampling dominates.
- Triangular original pair-identity storage for large `N`.
- Chunked/compressed project export to avoid duplicating very large normalized alignments.
- A source-supported alternative to the 256-fragment cap if real datasets reach it routinely.
- Real parallel triplet partitions in the optional pthread build; the current target is build
  plumbing only and does not divide scan work.
- Resumable in-progress round checkpoints. Completed analyses/review state are resumable now.

The non-threaded worker is the compatibility target because ordinary static hosting does not always
supply the COOP/COEP headers required by `SharedArrayBuffer`.
