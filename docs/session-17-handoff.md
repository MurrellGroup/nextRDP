# Session 17 handoff

## Checkpoint boundary

Session 17 advances the engine to `0.17.0-session-17` while retaining project schema
`org.rdp-web.project/v1alpha16`; no durable project fields changed. It resolves the reported
GENECONV probability-scope question, corrects the lifetime of the source multiple-comparison
factor, and activates a state-safe cyclic shortlist derived only from the supplied RDP5 source.
No alternate RDP implementation was consulted.

## GENECONV probability conclusion

- `GCCalcPValP2` calculates the raw Karlin–Altschul probability from the selected scored fragment.
- Ordinary `GCXoverD(0)` stores `PVals × MCCorrection` in `XOverList.Probability`. That corresponds
  to WebRDP's **Project corrected** value. WebRDP additionally exposes the pre-factor value as
  **Raw KA P**.
- Consequently, an RDP5 value near `1e-14` and WebRDP raw value near `1e-18` are expected if the
  initial scan plan contains about 10,000 opportunities. The event card now prints the exact
  multiplier beside the corrected value.
- BURT/`PolishBP` is later in the workflow. It may refine the reconciled event boundaries, but does
  not recalculate fragment score, raw KA P, or corrected discovery P.
- The WASM implementation preserves supplied float narrowing and KA tail branches. Its bounded
  bisection route is used only if the supplied Newton iteration becomes unstable; the exported
  `geneconvNumericalFallbackTracks` counter identifies that case.
- The DLL's `NDiff[3] == 1` comparison typo is treated as the intended assignment, matching the
  guards on outer tracks 4 and 5 and avoiding undefined division-by-zero behavior on a degenerate
  outer-track-3 profile.

Without Darren's exact saved project/dataset, the raw-versus-corrected multiplier is the strongest
explanation rather than a fixture-level proof. Comparing RDP5's stored value to WebRDP's
**Project corrected** value and checking `geneconvNumericalFallbackTracks` will distinguish the
ordinary case from the two explicit numerical safety boundaries.

## Cyclic shortcut implemented

The supplied source retains detected summaries in `XOverList`/`XOverDefine`, copies actionable
records into `BestXOList`, uses `Worthwhilescan` method bits and `StoreLPV` thresholds for later
analysis lists, and can build two-anchor substitution lists in `FindBetterRecSignal`.

The browser implementation now retains signal-bearing summaries by exact mutable working triplet.
After erasure it marks every changed row and newly created fragment dirty. Triplets touching a
dirty row run all enabled kernels; unchanged triplets replay cached signals and skip stable method
work, including known-empty summaries without storing an entry for every empty triplet. This keeps
memory proportional to threshold-passing summaries rather than `choose(N,3)`.

First-round 3SEQ summaries are deliberately excluded because its next pass activates supplied
post-erasure split handling. Unchanged triplets refresh only 3SEQ on that transition. Later 3SEQ
summaries can be reused. Manual correction/rejection rebuilds clear all shortlist state and start
with a full scan.

`MakeMCCorrection` is now frozen from the initial scan plan, as in `BuildFirstXOList`. Current-round
`totalTriplets` remains dynamic for honest progress while `correctionTests` remains stable for all
probability calculations and cache validation.

## Verification completed

- TypeScript/Vite production build passed.
- Emscripten 5.0.1 single-worker WASM build passed.
- Source-contract and delimiter-balance checks passed.
- GitHub Pages artifact structure, hidden `.nojekyll`, relative assets, WASM header, and real FASTA
  upload passed.
- The production verifier's deterministic two-event scan passed: 476 scheduled triplets, 280
  triplet-kernel evaluations, 196 unchanged summaries reused, 196 RDP scans skipped, 12 cached
  signals replayed, and the initial correction fixed at 120 while later working schedules changed.
- A temporary full-rescan build evaluated all 476 kernels on the same fixture. Its event roles,
  breakpoints, rounds, support-signal order, individual signal triplets, and p-values matched the
  optimized 280-kernel run exactly; the final deliverable was rebuilt with the shortlist enabled.
- A separate all-five-family optimized run completed 22 cyclic events: 131,488 scheduled triplets,
  32,911 triplet-kernel evaluations, 98,661 unchanged summary reuses, 1,199 cached signals, and
  493,221 individual method scans skipped. Its selected-result SHA-256 matched the temporary
  full-rescan build across all 22 events and 422 retained signals.

## Remaining boundary

The shortcut avoids expensive alignment/profile/probability kernels but still walks the lightweight
eligible-triplet schedule to preserve deterministic order, query/reference constraints, same-origin
fragment exclusion, progress, and cancellation. Porting the narrower native
`FindBetterRecSignal` two-anchor list can remove more scheduler overhead, but should wait for
authorized golden fixtures covering exact tie order, `BanTriplet`, masks, grouped references, and
fragment aliases.

The GENECONV and multi-method algorithms remain source-shaped active rather than desktop-parity
validated. Darren's exact comparison dataset/project would be the most useful next fixture if the
**Project corrected** values still differ materially after this checkpoint.
