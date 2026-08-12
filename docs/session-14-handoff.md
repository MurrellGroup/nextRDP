# Session 14 handoff

## Checkpoint boundary

Session 14 advances the source checkpoint to `0.14.0-session-14` and project schema
`org.rdp-web.project/v1alpha14`. It adds the supplied ordinary automated 3SEQ path as a fifth
primary discovery stream beside RDP, GENECONV, MaxChi, and CHIMAERA. Any enabled stream can author
the strongest corrected signal entering cyclic event selection, tract erasure/fragment re-entry,
three-set reconciliation, BURT polishing, ordered review/repair, and export.

The implementation was derived only from the supplied RDP5 manual, VB source, DNA5 DLL source, and
companion DNA DLL source. Its active trace follows `FindSubSeqTS`, `FindSubSeqTS2`, `Seq3PVals`,
`Get3SeqPvalC`, `GetTSPVal`, `CheckwrapC`, `SiegmundDiscrete`, `SwapRound`, and `TSXOver`. No alternate
RDP or 3SEQ implementation was consulted.

Per the project instruction, no supplied/native source, C++, WebAssembly, TypeScript compiler,
Vite build, preview server, dependency install, or project/runtime test was executed. “Implemented”
means source and cross-layer contracts exist; it does not mean compiled or native-parity validated.

## Implemented this session

- Three source-order target rotations using the shared one-pass equality profile. The pinned
  target-to-parent pair-slot maps are `{0,2,1}` and `{1,0,2}`; all-different, parent-equal, and
  missing sites are excluded, and the source's four-site minimum is retained.
- Source-shaped maximum descent/ascent discovery with the supplied probability-before-`CheckwrapC`
  call order, prefix extension, beginning-site advance, and linear conversion. Pre-wrap probability
  and post-wrap boundary excursions are retained separately. A strictly smaller reverse-tail
  probability swaps parents, counts, excursions, and coordinates; ties retain the descent.
- An exact hypergeometric random-walk tail evaluator that replaces the desktop four-dimensional
  table with a compact dynamic program, narrows results to the supplied `Single` precision, caches
  `(m,n,k)` states, and caps both transitions and cache growth.
- The supplied large-profile `SiegmundDiscrete` path with `ApproxNu`, normal CDF/PDF helpers, and a
  bounded scaled-exact fallback when the approximation leaves `(0,1)`. Exact and fallback counts
  remain separate in progress, JSON, restore, and review.
- Literal post-orientation low-information exits, zero-P rejection, and the method-specific
  Dunn–Šidák correction above `p = 10^-15` with the supplied smaller-tail product route, using the
  active exploratory or query/reference opportunity count.
- Combined RDP → GENECONV → MaxChi → CHIMAERA → 3SEQ method-major tie ordering in the existing
  strongest-first cyclic scheduler. 3SEQ candidates retain provisional roles and enter the same
  event grouping, erasure, fragment, consensus, confidence, and export workflow.
- A settings switch enabled by default, method-specific correction explanation, live target-walk /
  exact / approximation / candidate counters, and a review warning when the bounded source
  approximation was used.
- Per-signal evidence for target rotation, walk direction, information-rich sites, both parent-
  match counts, probability and boundary excursions, raw/corrected probability, probability route,
  and split state.
- On-demand three-target random-walk plotting with a signed y-axis, zero-aware scaling, selected-
  tract highlighting, candidate-target legends, and disclosure for later original-alignment
  reconstructions.
- ABI and restore coverage: method code `4`, `rdp_restore_threeseq_discovery`, four authoritative
  workload counters, scan/restore enable state, and the Emscripten export. Restore rejects invalid
  counts, excursions, role/pair rotations, low-information calls, and unsupported split claims.
- Backward restoration: `v1alpha14` restores 3SEQ; all older schemas keep it disabled. Existing
  `v1alpha1`–`v1alpha13` method gates remain intact. Unknown saved method labels are rejected rather
  than silently restored as RDP.
- Static GitHub Pages workflow coverage remains in place for default-branch or manual Actions
  deployment, relative project paths, a single-worker compatibility WASM build, `.nojekyll`, and
  artifact validation.

## Fidelity boundary

3SEQ is **source-shaped active unvalidated**. Ordinary discovery and its cyclic scheduling are
represented end to end. The exact dynamic program is intentionally more memory efficient than the
desktop table but needs authorized native saved-output comparison for float rounding and cutoff
edges. The supplied `FindSubSeqTS2` coordinate alias, `CheckSplit3Seq`/`SubPVal` post-erasure split
and re-probability path, manual 100-permutation envelopes, and representative/final-list
corroboration are not claimed. `missingDataSplitApplied` therefore remains false.

MaxChi, CHIMAERA, GENECONV, BURT/BenHMM, query/reference, and late reconciliation also remain
source-shaped rather than parity validated. No result should be called desktop-equivalent until the
authorized native golden corpus exists.

## Static checks allowed for this checkpoint

The source checker covers C ABI/header/CMake/worker names and arities, version/schema gates, the
3SEQ pair-slot maps and four-site floor, exact/fallback bounds, scheduler/options/counters, restore
evidence, UI/plot/CSV surfaces, old-schema behavior, and source-trace documentation. Delimiter
balance, JavaScript syntax, shell syntax, YAML parsing, stale-current-version scans, Pages workflow
contracts, and archive integrity are also checked without executing project code.

GitHub Pages remains configured for **Settings → Pages → Source: GitHub Actions**. The workflow
builds the compatibility WASM module and Vite site remotely, validates source/types/artifacts,
uploads `dist` including `.nojekyll`, and deploys only from the repository default branch or a
manual dispatch.

## Known open work

1. Create authorized native 3SEQ fixtures spanning both orientations, target rotations, four-site
   and low-information exits, exact-table edges, Siegmund/scaled fallbacks, circular/linear bounds,
   all-different/missing states, Dunn–Šidák underflow, and combined event ordering.
2. Resolve the `FindSubSeqTS2` position alias and port `CheckSplit3Seq`/`SubPVal` without silently
   replacing source behavior with an inferred missing-data policy.
3. Trace and port 3SEQ permutation/manual display and late-list corroboration as separately
   labelled modes.
4. Continue with BOOTSCAN or SISCAN from the supplied sources only, while extending the authorized
   golden corpus for all active methods.
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
