# Session 19 handoff — primary BootScan discovery and pair cache

Version: `0.19.0-session-19`
Project schema: `org.rdp-web.project/v1alpha17`

## Completed in this checkpoint

- Added source-shaped automated distance-mode BootScan as a primary discovery family. The common
  cyclic order is now RDP → GENECONV → BootScan → MaxChi → CHIMAERA → 3SEQ.
- Traced and implemented the supplied `BSXoverR`, `SEQBOOT2`, `FastBootDist`, `GetPltVal`,
  `ScanBSPlots`, `FindBeginBS`, `FindEndBS`, `MakeBSEvent`, `BSSubSeq`, `MakeScoresBS`, and
  `ProbCalc` path without consulting alternate reference implementations.
- Preserved the seeded Microsoft-rand bootstrap weights, replicate-zero original window,
  source-generated-but-discarded tail draw per site, Jukes–Cantor saturation behavior, strict
  closest-pair votes, circular window traversal, and source-shaped support-region continuation.
- Added a bounded 64 MiB FIFO pair-profile cache. Shared two-sequence window/bootstrap distances are
  reused across triplets in the same round; cache telemetry is exported and shown during scanning.
  The cache is invalidated after erasure. The separate XOverList/XOverDefine and BestXOList-style
  triplet shortlist continues to replay unchanged whole-triplet summaries across rounds.
- Preserved `MakeScoresBS`'s `XPosDiff` invariant-boundary convention, its raw binomial tail,
  informative-length multiplier, 169-site scaling/exponent, and project correction as distinct
  values. BURT/BenHMM remains a later coordinate-polishing step and does not recalculate the
  BootScan discovery p-value.
- Added method-labelled signals, BootScan discovery evidence, closest-pair support plots, scan/cache
  counters, JSON/CSV export, hardened `v1alpha17` restore, and Windows 95-styled settings/progress/
  review controls. Imports remain backward compatible through `v1alpha1`.

## Verification contract

The source contract checks C ABI arity/export parity, v17 import/export coverage, method order,
cache bounds/invalidation, BootScan evidence fields, worker restore, UI controls, trace disclosures,
and Pages smoke-test presence. A native host check links every ported C++ translation unit before
running the deterministic BootScan/cache fixture. That host-linked check, source contracts,
TypeScript, script syntax, and the production Vite build pass in this checkpoint. The production
artifact verifier instantiates the real WASM module, loads FASTA through `HEAPU8`, preserves the
exact cyclic-shortlist/graceful-stop regressions, and adds a primary-BootScan/cache regression. It
runs in GitHub Actions because this checkpoint environment does not include Emscripten.

## Honest boundary

Primary distance-mode BootScan is source-shaped active and unvalidated. It is not claimed to match
every desktop event until authorized native saved-output fixtures compare candidate endpoints,
provisional roles, ordering, raw probability, and corrected probability. Tree/similarity/
permutation/manual BootScan modes, literal full edge-warning/catalogue behavior, and full late event
reconstruction remain pending. No alternate RDP implementation was consulted.

The supplied native projects were read only and were not compiled.
