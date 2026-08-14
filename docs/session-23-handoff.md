# Session 23 handoff — cyclic pruning and browser throughput

Version: `0.23.0-session-23`
Project schema: `org.rdp-web.project/v1alpha19`

Session 23 is a behavior-preserving performance checkpoint based on Darren's pruning guidance and
the supplied `DoRDP`/`DropSeqs`, `XOverList`/`XOverDefine`, `BestXOList`, and `Worthwhilescan`
paths. No alternate RDP implementation was consulted, and no supplied native source was compiled.

## Added in this checkpoint

- Made an unchanged triplet that was clean on its first complete screen a permanent negative
  shortlist entry. Later inner/outer rounds still enumerate it for deterministic schedule and stop
  semantics, but invoke no discovery method. The special first post-erasure 3SEQ split refresh is
  now limited to unchanged triplets that previously carried a 3SEQ-capable signal; it never
  overrides the clean-triplet rule.
- Added `DropSeqs`-shaped removal for fragment rows created by the preceding event. After exactly
  one complete follow-up round, an event-free fragment is removed before the next event can add
  more rows. Removal swaps the last working row into the vacant slot and remaps exact working-row
  provenance in live signals, the current `XOverList`-equivalent map, the carried
  `BestXOList`-equivalent shortlist, and dirty-row state.
- Kept same-origin fragment combinations out of the analytical batch budget. They are still
  traversed in deterministic order but no longer cause otherwise empty worker/WASM round trips.
- Maintained valid/missing-site summaries and deterministic state fingerprints during tract
  erasure. Active-row refresh is now `O(N)` instead of rescanning `O(NL)` state bytes, fragment data
  is collected sparsely before allocation, and duplicate fragments use a fingerprint prefilter
  followed by an exact byte comparison.
- Replaced prior-event tract marking in breakpoint polishing and `CheckEnds` context construction
  with inclusive range-difference unions. The result is still point-for-point identical, while
  the repeated marking work changes from `O(EL)` to `O(E + L)` per affected role/context.
- Limited routine progress serialization/posting/rendering to once every 500 ms. Worker scan
  budgets adapt toward a 40 ms slice and use `scheduler.yield()` where available, retaining prompt
  Stop handling without paying for a fixed JavaScript/WASM crossing every 512 triplets.
- Enabled release `-O3`, link-time optimization, and WebAssembly SIMD without enabling unsafe
  floating-point reassociation.

## Regression boundary

The focused ten-sequence/two-mosaic regression completes two cyclic events and reports 162 clean
triplets pruned, 168 method calls skipped, 18 same-origin schedule combinations bypassed, and two
event-free fragments removed with swap/reindex compaction. Its selected event roles, coordinates,
rounds, support order, and p-values retain the Session 22 digest
`5ad90dbeeecd3ea531d52455dd3ded89498c8d0aeefc5d73c2885e451648e6fa`.

The linked BootScan/public-ABI, SISCAN, supplied event-tree, and PHYLPRO/brute-force host gates also
pass. The Emscripten build and instantiated production-WASM smoke tests remain authoritative in
GitHub Pages Actions because this local toolchain does not provide Emscripten.

## Schema and fidelity status

The project schema stays `v1alpha19`: the new counters are transient progress telemetry, and the
analysis/result shape is unchanged. These shortcuts require exact working-row identity and a fixed
initial correction factor; any erased/new row runs fresh kernels, and edit/rejection rebuilds clear
the shortlist. Native golden comparison remains required before calling the broader port parity
validated.

The recombinant-identification families Darren listed (`ParsimonyO`, `ParsimonyI`, `dmax(VISRD)`,
`SetDistT`, `SetDistP`, O:E, `distrank`, `conflict`, `oucheck`, and `TrpScore`) remain a separate
functionality phase; this checkpoint does not claim they were completed merely by making the scan
faster.

See [`native-cyclic-shortlist-trace.md`](native-cyclic-shortlist-trace.md) for the supplied-source
mapping and exact cache invalidation boundary.
