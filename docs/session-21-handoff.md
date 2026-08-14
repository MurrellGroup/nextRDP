# Session 21 handoff — supplied SISCAN workflow

Version: `0.21.0-session-21`
Project schema: `org.rdp-web.project/v1alpha19`

Session 21 adds the ordinary supplied SISCAN path without consulting alternate RDP implementations.
Primary discovery is optional and confirmation is enabled by default, matching the RDP5 workflow.
No alternate RDP implementation was consulted.

## Delivered

- `wasm/src/siscan.cpp` and `.hpp` implement the supplied `SSXoverC` → `GetSSOL` → `Get3Score` →
  `GetPScores2` → `DoPerms3` → `MakeZValue2` → `DoSums` → `FindMaxZ` → `ShrinkRegionC` path.
- The nearest fourth sequence uses a round-wide source-shaped WPGMA/cophenetic context. Original
  identities prevent a cyclic fragment from becoming its own outlier, and disabled rows are
  removed before pair scoring/tree merging rather than consuming context work.
- MakeVRand uses the Microsoft CRT stream and an extensible flat byte prefix shared by triplets and
  later rounds. Only the state-dependent WPGMA context is invalidated after erasure.
- The existing XOverList/BestXOList-style triplet shortlist bypasses SISCAN for replayable unchanged
  triplets; affected triplets still run the source-shaped kernel.
- The default fast scan preserves the supplied `QuickCheckB` missing-braces control-flow quirk.
- Primary SISCAN candidates enter the same strongest-first cyclic order:
  RDP → GENECONV → BootScan → MaxChi → CHIMAERA → SISCAN → 3SEQ.
- Fixed-region SISCAN confirmation runs on the event representative and finalized distance-list
  rows without changing reconciled roles or coordinates.
- Discovery and confirmation retain normal-tail, region-adjusted, window-adjusted, and
  project-corrected P-values separately, plus category family/score, outlier, Z, and permutation
  workload.
- The review screen exposes discovery, representative, finalized-list, and a signed three-curve
  SISCAN plot that keeps the greatest-absolute eligible category for each sister pair/window.
  Project JSON and expanded CSV retain the same audit trail.
- Project schema `v1alpha19` restores method code 6, settings, discovery traces, scan work, WPGMA
  context, and random-prefix generation counters. Older schemas conservatively restore SISCAN
  disabled.
- `scripts/verify-siscan-core.cpp` checks the planted tract, exact seed-3 Microsoft-CRT prefix and
  replay, WPGMA cache, recheck, plot, cyclic invalidation, and same-browser-context restart. Cache
  telemetry is reset for every new analysis. GitHub Actions runs this gate before the production
  Emscripten/Pages build.

## Native defaults represented

- primary discovery off;
- confirmation on;
- 200-site window, 20-site step;
- 100 scan permutations and 1,000 final permutations;
- seed 3;
- gap stripping, all one/two/three-variable categories, nearest outlier, fast scan.

## Deliberate boundaries

- SISCAN is source-shaped and host-regression tested but not native-golden validated.
- Random-outlier, most-distant-outlier, alternative gap/category modes, and manual plot controls are
  not exposed as claimed-parity modes.
- LARD remains unavailable because the supplied VB project delegates it to an external executable
  whose source was not supplied.
- The manual's complete 15 category plus 9 summed-category plot remains open; the compact review
  plot is not claimed as full manual-plot parity.
- PHYLPRO and the remaining manual/full late-catalogue paths remain open.
- No supplied native source, C++, WebAssembly, or TypeScript compiler was executed while tracing the
  references. Only the independently ported code is compiled by project checks.

See [the detailed source trace](native-siscan-discovery-trace.md) before native comparisons.
