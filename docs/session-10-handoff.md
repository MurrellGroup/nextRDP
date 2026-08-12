# Session 10 handoff

## Checkpoint boundary

Session 10 advances the source checkpoint to `0.10.0-session-10` and project schema
`org.rdp-web.project/v1alpha10`. It adds ordinary-triplet MaxChi exploratory discovery beside the
existing RDP discovery stream and carries either method through the complete browser workflow:

dataset → curation/settings → combined cyclic scan → event reconciliation → ordered review/repair →
project/CSV/alignment exports.

Per the project instruction, no native source, C++, WebAssembly, TypeScript compiler, Vite build,
preview server, dependency install, or project/runtime test was executed. “Implemented” below means
source and cross-layer contracts exist; it does not mean compiled or native-parity validated.

## Implemented this session

- Source-shaped `MCXoverF` discovery in `wasm/src/maxchi.cpp`:
  - all three variable-site pair profiles and native critical-difference screen;
  - fixed/fallback window selection;
  - current triplet `MissingData`, prior-erasure, and linear-end window bans;
  - raw boundary/pair peak priority;
  - symmetric grow, `FindSide`, left/right breakpoint optimization, and tract mapping;
  - literal `SmoothChiValsP` terminal-padding/off-by-one behavior;
  - completed and rejected peak destruction;
  - accepted-hit waste reset, three-wasted stop, and 100-attempt cap; and
  - raw, within-triplet, and current-round corrected probabilities.
- A performance-preserving peak scheduler: three rolling χ² profiles, linear-time heap construction,
  and lazy invalidation of destroyed cells instead of a complete maximum rescan after every peak.
  When both methods are enabled, one alignment-byte pass prepares the RDP categories, MaxChi
  matches/variable prefixes, and MaxChi missing/erasure map. The supplied global raw-maximum gate
  runs before smoothing/retry accounting, avoiding that work for an insignificant triplet.
- Combined RDP/MaxChi strongest-first cyclic scheduling. Signals retain their method; either can
  anchor a reconciled event; support can span methods; the selected group tract is erased before a
  fresh complete pass.
- Method-aware result and project contracts, progress counters, signal plots, event badges, anchor
  discovery diagnostics, and CSV columns.
- C ABI, CMake export list, worker calls, and restore logic for MaxChi settings, signal methods, and
  per-signal discovery traces. Completed-project reload preserves cumulative triplets, scan rounds,
  final-round processed/total triplets, the terminal reason, and all MaxChi
  profile/peak/candidate/limit counters. Replaying its events no longer inflates saved cumulative
  work.
- Backward project loading through `v1alpha1`. Projects through `v1alpha9` intentionally restore
  MaxChi discovery disabled, preserving the detection semantics under which they were saved.
- Settings for enabling ordinary MaxChi discovery and selecting its requested window, defaulting to
  enabled and 70 sites.
- A method-aware review plot: RDP displays pair identity; MaxChi displays the three χ² profiles.
  Bounded thinning always retains both breakpoints, the selected MaxChi peak, and every pair
  maximum. A later-round or fragment-assisted plot is explicitly labelled as an original-alignment
  reconstruction because compact checkpoints retain its statistics, not every historical working
  profile point.
- A dedicated MaxChi anchor card with attempt, tract side, peak/pair, initial/grown windows, critical
  difference, flank χ², raw/within/project probabilities, and filter status.
- Visible scan/review warnings when the supplied 100-peak bound leaves positive raw peaks unexplored.
- A source trace in `docs/native-maxchi-discovery-trace.md` and expanded validation fixtures for the
  newly active lifecycle.

The Session 9 `FastRecCheckMC2` MaxChi representative/finalized-list recheck remains active and
separate. It corroborates an already reconciled event and does not reposition detected or manually
edited coordinates.

## Important fidelity decisions

- The active `Module5.bas` `MCXoverF` path is authoritative for this phase; older/sibling VB copies
  with disabled optimizer branches were not substituted.
- `FindMChiP`'s strict raw tie order is retained. A heap changes how the next maximum is found, not
  the intended chi/boundary/pair priority.
- The smoothing routine's twelve-term/eleven-divisor behavior and terminal zero padding are
  intentional supplied-source quirks.
- Destruction traversal is bounded. Undefined native memory reads are not reproduced.
- MaxChi's raw triplet role is provisional. It is visibly retained and then arbitrated by the shared
  three-role evidence/consensus path; native role fixtures remain mandatory.
- The ordinary triplet path is not the manual pair/doublet or permutation path. Those modes remain
  absent rather than being approximated under the same label.

## Static checks allowed and performed

The source-only checker verifies:

- every public `rdp_api.h` entry has a keepalive definition, CMake export, and worker contract/call;
- package, native engine, result, and project schema versions agree;
- `v1alpha1`–`v1alpha10` import coverage and the older-project MaxChi-disabled boundary exist;
- MaxChi kernel/options/progress/result/restore/UI contract markers exist; and
- C++/TypeScript delimiter balance across the checked source set.

Shell/YAML/archive syntax and integrity may also be inspected without executing project code. The
GitHub Pages workflow remains configured for **Settings → Pages → Source: GitHub Actions**. Its first
actual run will install/build remotely only when the user pushes or manually dispatches it.

## Known open work

1. Establish supplied-desktop golden fixtures for MaxChi profile values, native lookup rounding,
   `MaxX = 0`, smoothing padding, exact ties, missing/linear bans, grown windows, left/right
   optimization, accepted/rejected destruction basins, three/100 retry behavior, and combined
   cross-method event order.
2. Validate and, where golden evidence requires it, refine preliminary MaxChi role selection.
3. Trace and port supplied MaxChi permutation and manual pair/doublet modes as distinct workflows.
4. Add the remaining detection-method families and their post-group method-stack rechecks without
   hiding per-method evidence or overstating full native weighting parity.
5. Run the first tiny, resource-limited compile/runtime checkpoint only after explicit permission.

## Resume map

- MaxChi discovery kernel: `wasm/src/maxchi.cpp`, `wasm/src/maxchi.hpp`
- Combined scanner/reconciliation/results/CSV/plots: `wasm/src/rdp_method.cpp`,
  `wasm/src/rdp_method.hpp`
- C ABI and project serialization: `wasm/include/rdp_api.h`, `wasm/src/rdp_api.cpp`
- Worker/import boundary: `src/workers/analysis.worker.ts`
- Shared types: `src/lib/types.ts`
- Settings/progress/review/plot/export UI: `src/components/SettingsStep.tsx`,
  `ScanStep.tsx`, `ReviewStep.tsx`, `SignalPlot.tsx`, `ExportStep.tsx`
- Static source contracts: `scripts/check-source-contract.mjs`,
  `scripts/check-source-balance.mjs`
- Supplied-source mapping: `docs/native-maxchi-discovery-trace.md`
- Golden corpus: `docs/validation-plan.md`

The checkpoint remains source-only and should continue to be treated as unvalidated until the
authorized golden/runtime phase.
