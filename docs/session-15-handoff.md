# Session 15 handoff

## Checkpoint boundary

Session 15 advances the source checkpoint to `0.15.0-session-15` and project schema
`org.rdp-web.project/v1alpha15`. It completes the supplied discovery-time post-erasure 3SEQ route
through `FindSubSeqTS2`, `CheckSplit3Seq`, and `SubPVal`, then adds the non-coordinate-changing
representative/finalized-list `TSXOver(1)` Findall recheck. The ordinary automated 3SEQ path from
Session 14 remains active. No alternate RDP or 3SEQ implementation was consulted.

Per the project instruction, no supplied/native source, C++, WebAssembly, TypeScript compiler,
Vite build, preview server, dependency install, or project/runtime test was executed. “Implemented”
means source and cross-layer contracts exist; it does not mean compiled or native-parity validated.

## Implemented this session

- Resolved `FindSubSeqTS2`'s position alias directly from the supplied DLL: after an informative
  site, `XPosDiff[x]` is the number of retained sites at or before `x`. The compact browser route
  uses `upper_bound` over retained coordinates instead of allocating another alignment-sized map.
- Added source-shaped beginning- and ending-side missing/erased-state searches. Each candidate
  sub-tract obtains the supplied `SubPVal` maximum-minus-minimum excursion and passes it through the
  same bounded exact hypergeometric random-walk/`SiegmundDiscrete` probability route as ordinary
  discovery.
- Preserved `TSXOver` ordering: the original corrected candidate must pass first; the selected
  orientation is split; the reverse orientation is attempted only when its unsplit P is strictly
  lower than the trimmed selected P; and the method-specific Dunn–Šidák correction plus the source
  call gate run again.
- Reused `MaxChiWorkspace::triplet_missing_data` and the compact three-target walk. This adds no
  alignment-byte pass, bounds missing-run searches to the alignment, and reuses the exact
  probability cache.
- Kept the compact exact DP state and every accumulation in supplied `Single` precision, reducing
  its WASM storage relative to extended precision while retaining the desktop table's numeric type.
- Retained split state, revised bounds, sub-tract probability excursion, raw/corrected probability,
  exact/fallback route, roles, and direction through project JSON, restore, CSV, and review UI.
- Corrected the bounded scaled `GetTSPVal` adaptation to use the supplied post-truncation exponent
  ratio and positive-underflow `1e-300` floor.
- Corrected visible and serialized trace order to `FindSubSeqTS` → `Seq3PVals`/`GetTSPVal` →
  `CheckwrapC` → `TSXOver` and retained the distinct probability-before-wrap behavior.
- Advanced project import through v15. Pre-v14 projects still keep 3SEQ disabled; v14 supports the
  ordinary path but rejects a stray split claim that only v15 can author.
- Extended static source contracts across the kernel, fused scheduler, restore gate, schema,
  worker, UI, trace, and Pages package.
- Added `threeseq_recheck_prepared` for the supplied `FindallFlag = 1` route. It performs the normal
  initial selection/low-information/corrected gate, forces `CheckSplit3Seq` on both orientations,
  applies the source's distinct correction behavior after `SwapRound`, and accounts for the
  inverse-parent/inverse-interval list copy created for every qualifying orientation.
- Runs that Findall evidence on the event representatives and every finalized nonrepresentative
  distance-list triplet. JSON, CSV, review cards, the late matrix, and project-status metadata retain
  workloads, list-entry counts, best target/direction, bounds, excursions, exact/fallback route,
  split state, raw/corrected P, and the source cutoff result without changing reconciled coordinates.

## Fidelity boundary

The ordinary and post-erasure-split 3SEQ discovery paths and late Findall recheck are
**source-shaped active unvalidated**.
Authorized desktop saved-output fixtures still need to cover coordinate-zero/one aliases, wrapped
missing runs, unchanged versus moved split bounds, strict orientation retry, exact-table edges,
float narrowing, correction underflow, both Findall orientations, inverse list copies, and native
list-cap replacement behavior. Manual 100-permutation display envelopes and full late event-
catalogue reconstruction remain separate open work.

MaxChi, CHIMAERA, GENECONV, BURT/BenHMM, query/reference, and late reconciliation likewise remain
source-shaped rather than parity validated. No result should be called desktop-equivalent until the
authorized native golden corpus exists.

## Static checks allowed for this checkpoint

The source checker covers C ABI/header/CMake/worker names and arities, version/schema gates,
3SEQ pair-slot maps and four-site floor, exact/fallback bounds and scaling, pre-/post-wrap excursions,
split coordinate mapping, missing-run bounds, orientation retry, scheduler/options/counters,
restore evidence, UI/plot/CSV surfaces, old-schema behavior, and source-trace documentation.
Delimiter balance, JavaScript syntax, shell syntax, YAML parsing, stale-current-version scans,
Pages workflow contracts, and archive integrity are also allowed without executing project code.

GitHub Pages remains configured for **Settings → Pages → Source: GitHub Actions**. The workflow
builds the compatibility WASM module and Vite site remotely, validates source/types/artifacts,
uploads `dist` including `.nojekyll`, and deploys only from the repository default branch or a
manual dispatch.

## Known open work

1. Validate ordinary and post-erasure 3SEQ against authorized native fixtures, especially split
   endpoint aliases, wrapped missing runs, reverse-orientation selection, scaled probability, and
   second-stage correction gates.
2. Validate the representative/finalized-list `TSXOver(1)` pass across target rotations, both
   orientations, split/no-split cases, correction underflow, inverse interval copies, and list caps.
3. Trace and port the manual 100-permutation 3SEQ envelope as a separately labelled review mode.
4. Continue the supplied BOOTSCAN or SISCAN family and its late corroboration, without consulting
   alternate reference implementations.
5. Run the first tiny resource-limited compile/runtime checkpoint only after explicit permission.

## Resume map

- 3SEQ kernel/contracts: `wasm/src/threeseq.hpp`, `wasm/src/threeseq.cpp`
- Fused scheduler/results/plot/CSV: `wasm/src/rdp_method.hpp`, `wasm/src/rdp_method.cpp`
- C ABI/schema/restore: `wasm/include/rdp_api.h`, `wasm/src/rdp_api.cpp`
- Worker scan/restore: `src/workers/analysis.worker.ts`
- Types/defaults: `src/lib/types.ts`, `src/App.tsx`
- Settings/progress/review/plot: `src/components/SettingsStep.tsx`,
  `src/components/ScanStep.tsx`, `src/components/ReviewStep.tsx`,
  `src/components/SignalPlot.tsx`
- Supplied-source trace: `docs/native-threeseq-discovery-trace.md`
- Static contracts: `scripts/check-source-contract.mjs`, `scripts/check-source-balance.mjs`

This checkpoint remains source-only and must be treated as unvalidated until the authorized native
golden/runtime phase.
