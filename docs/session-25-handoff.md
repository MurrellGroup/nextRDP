# Session 25 handoff — cyclic fragment evidence reconstruction

Version: `0.25.0-session-25`
Project schema: `org.rdp-web.project/v1alpha19`

This checkpoint starts from the reattached Session 24 source archive (SHA-256
`510a6782c1dd818a91ebed4e3dd7bb28fdbfb6e763b6d2ce1bf6464404ce2b3b`). It
reconstructs the source-backed cyclic fixes that remained recoverable after the prior scratch
workspace was lost. It contains no round number, sequence name, breakpoint coordinate, or test-file
special case. The supplied desktop sources remained read-only and were not compiled. No alternate RDP implementation was consulted.

## Source semantics restored

- Event reconciliation now keeps the supplied internal `ISeqs` order: daughter, minor parent, major parent.
  Public JSON, UI, and exports still use recombinant, major parent, minor parent.
  The distinction matters in asymmetric `FinalTrim` branches.
- `TraceSub`-equivalent evidence aliases are recovered from the exact working triplets stored in
  cached XOverList-style signal summaries. The strongest support supplies a general alias for each
  original sequence and the anchor signal wins. Old project records without exact provenance keep
  the fragment-event fallback.
- When fragment aliases are needed, grouping kernels receive an analytical state view with the
  retained fragment rows overlaid into their original identity slots. This applies to all event
  candidates, not merely the three representatives, without changing public identities.
- The selected-role outlier rank uses the full fragment-expanded working row count (`NextNo`) and
  denominator, matching the supplied late-grouping path instead of ranking only original rows.
- Ordinary candidate identities are compared with the reported representative identities during
  final cleanup, rather than with transient working-fragment indices.
- BURT triplet preparation uses the same daughter/minor/major source ordering. BURT remains a later
  breakpoint-polishing step and does not recalculate discovery p-values.

## Performance boundary

The Session 24 XOverList/BestXOList-style clean-triplet and signal-summary reuse remains unchanged.
The reconstructed evidence path runs only when a retained fragment alias exists, so primary scanning
and pre-fragment events do not allocate an analytical copy. That copy contains only row-count,
length, and state bytes; it deliberately omits sequence strings, summaries, and the `O(N^2)` pair
similarity table.

The expanded `NextNo` rank is more work than the incorrect original-only rank. Its symmetric pair
matrix now requests each unordered sequence pair once, adds the result to both rows in the same
per-row order, and reuses those accumulated rows for the three representatives. This roughly halves
pair-cache calls and removes six representative re-sums without changing the intended floating-point
summation order. The expensive JC cache miss was already symmetric, so this is not a twofold
whole-stage or whole-run claim. Preventing a false event choice may also avoid downstream cyclic rounds and rescans,
but no whole-dataset speedup is claimed without matched native/Web benchmarks.

## GitHub Pages correction

The artifact verifier no longer assumes that Emscripten's private pthread helper basename is exactly
`rdp-core-threads.worker.js` or `.mjs`. It discovers emitted `.worker.js`/`.worker.mjs` helpers and
requires the threaded loader to reference at least one real emitted file. This retains the deployment
safety check while allowing CMake output-name and Emscripten-version naming changes.

## Validation and confidence

- Source and delimiter contracts pass.
- BootScan, SISCAN, event-tree, PHYLPRO, cyclic-pruning, and deterministic multicore host-linked
  gates compile and pass.
- The unrelated cyclic fixture retains selected-result digest
  `5ad90dbeeecd3ea531d52455dd3ded89498c8d0aeefc5d73c2885e451648e6fa` while still pruning 162 clean
  triplets and skipping 168 method scans.
- The local environment has neither Emscripten nor installed npm dependencies, so the dual-WASM,
  TypeScript, Vite, and instantiated Pages artifact gates remain GitHub Actions checks.
- The earlier opaque Darren comparison files were not reattached in this work block, so the native
  cycle-eight comparison was not rerun here.

Confidence is high that the ordering, fragment-alias, and expanded-rank changes implement general
supplied-source behavior rather than fixture repair. Confidence is moderate that they improve later
cyclic selection broadly. Full cyclic parity is not claimed: the last available targeted run still
continued beyond the desktop event count and diverged again after the repaired cycle. Multiple
independent native goldens remain required.
