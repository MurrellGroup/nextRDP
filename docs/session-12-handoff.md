# Session 12 handoff

## Checkpoint boundary

Session 12 advances the source checkpoint to `0.12.0-session-12` and project schema
`org.rdp-web.project/v1alpha12`. It adds source-shaped ordinary-triplet CHIMAERA discovery to the
existing RDP, MaxChi, exploratory/query-reference, cyclic reconciliation, ordered review/repair,
and export workflow.

The implementation was derived only from the supplied RDP5 manual, VB source, DNA5 DLL source, and
companion DNA DLL source. The principal trace is manual §8.5.1 plus
`MathFuncsDll.cpp::FindSubSeqDP3/6`, `FastRecCheckChim`, `AlistChi`, and the active VB `CXoverA`
dispatch. No alternate reference implementation was consulted.

Per the project instruction, no supplied/native source, C++, WebAssembly, TypeScript compiler,
Vite build, preview server, dependency install, or project/runtime test was executed. “Implemented”
means source and cross-layer contracts exist; it does not mean compiled or native-parity validated.

## Implemented this session

- A CHIMAERA kernel that rotates all three possible recombinant targets, builds the manual's
  information-rich binary parent-match strings, scans raw χ² peaks, grows windows, selects a tract
  side, optimizes boundaries, destroys peak basins, and respects three-wasted/100-attempt bounds.
- A fused preparation path: the existing one-pass triplet profile records RDP categories, all three
  MaxChi equality tracks, the alignment-coordinate prefix, and missing/erasure state. CHIMAERA
  derives each target profile from those cached bytes instead of rereading the alignment three times.
- Combined RDP/MaxChi/CHIMAERA strongest-first cyclic scheduling. Every signal retains its method;
  any enabled method can anchor an event before the shared support grouping, role evidence,
  breakpoint confidence, tract erasure, fragment re-entry, and full fresh pass.
- Independent CHIMAERA options with supplied default window 60, input validation, settings toggle,
  scan-plan copy, and progress counters for target profiles, peak attempts, emitted candidates, and
  target profiles reaching the retry cap.
- Full signal evidence in JSON and reloadable projects: target rotation, information-rich site
  count, peak/order/coordinate, initial and grown windows, critical difference, χ², raw/within-
  triplet/project probability, flank χ² values, inside/outside parent-one match rates, and filter flags.
- An on-demand CHIMAERA review plot containing exactly the retained target/parent-one trace, forced
  breakpoint and peak samples, and the same historical-profile disclosure used for later cyclic
  RDP/MaxChi signals.
- A dedicated responsive review evidence card and method-aware null-result, retry-cap, discovery,
  review, export, project, and CSV wording.
- A persisted event-level and CSV caution when MaxChi and CHIMAERA are the only supporting methods;
  the review does not present those closely related methods as independent confirmation.
- Source-shaped secondary `FastRecCheckChim` corroboration for each representative and finalized-
  list triplet. It rotates three targets, preserves raw/within/project probability scopes, reuses
  the MaxChi-prepared bytes, survives reload, and never changes event coordinates.
- ABI and restore coverage:
  - CHIMAERA enable/window arguments in scan and restore begin;
  - method code `2` for CHIMAERA signals;
  - a complete `rdp_restore_chimaera_discovery` record;
  - four authoritative workload counters in completed-project restore; and
  - the new keepalive symbol in the Emscripten export list.
- Backward restoration:
  - `v1alpha12` restores CHIMAERA settings, traces, counters, query/reference groups, and events;
  - `v1alpha10`–`v1alpha11` preserve MaxChi discovery but keep CHIMAERA disabled; and
  - `v1alpha1`–`v1alpha9` retain their earlier RDP-only discovery semantics.
- Dataset workflow repairs retained from the opening phase of this session: selected/filter-wide
  reference-group assignment, group compaction, exact scan eligibility feedback, corrected
  query/reference signal/event metadata, terminal zero-work progress, and live cyclic role counts.

## Fidelity and performance boundary

CHIMAERA is **source-shaped active unvalidated**. The binary profile and full event-producing
`CXoverA`-shaped lifecycle are present, but native lookup rounding, destruction edges, parent order,
and combined method ordering need authorized golden fixtures. The shared late consensus can change
the provisional target/parents, as the manual's event-identification stage requires.

The browser does not claim CHIMAERA permutation modes or full late-list method-stack event
reconstruction. Those remain explicit pending stages. Bounded `FastRecCheckMC2` and
`FastRecCheckChim` strongest-peak rechecks are active and must not be interpreted as two independent
confirmations; the manual cautions that MaxChi and CHIMAERA are closely related.

## Static checks allowed for this checkpoint

The source-only checker verifies API/header/CMake/worker alignment, package/native/result versions,
`v1alpha1`–`v1alpha12` import coverage, CHIMAERA source markers, fused preparation, options,
progress/results/types, full restore evidence, review/plot/export surfaces, and historical schema
gates. Delimiter, shell syntax, JavaScript-module syntax, YAML parse, Pages workflow markers, stale
version searches, and archive integrity are also checked without executing project code.

GitHub Pages remains configured for **Settings → Pages → Source: GitHub Actions**. The action builds
the single-worker WASM module and Vite site remotely, validates source/types and the final artifact,
uploads `dist` including `.nojekyll`, and deploys from the repository default branch or a manual
dispatch.

## Known open work

1. Establish authorized native CHIMAERA fixtures for all rotations, information-rich maps, window
   fallback, missing/linear bans, peak ties, growth, flank side, boundary optimization, destruction,
   retry counters, probabilities, roles, and cross-method cyclic order.
2. Establish the existing MaxChi, BURT/BenHMM, query/reference, and reconciliation golden corpora.
3. Trace and port CHIMAERA permutation and full late-list event-reconstruction modes as distinct
   source workflows beyond the active strongest-target statistic.
4. Port the next primary method family without reusing nonfaithful third-party reference code.
5. Run the first tiny resource-limited compile/runtime checkpoint only after explicit permission.

## Resume map

- CHIMAERA kernel and shared MaxChi primitives: `wasm/src/chimaera.hpp`, `wasm/src/maxchi.cpp`
- Scheduler/options/results/plots/CSV: `wasm/src/rdp_method.hpp`, `wasm/src/rdp_method.cpp`
- C ABI/project schema: `wasm/include/rdp_api.h`, `wasm/src/rdp_api.cpp`
- Worker scan/restore bridge: `src/workers/analysis.worker.ts`
- Web contracts/defaults: `src/lib/types.ts`, `src/App.tsx`
- Settings/progress/review/export: `src/components/SettingsStep.tsx`,
  `src/components/ScanStep.tsx`, `src/components/ReviewStep.tsx`,
  `src/components/SignalPlot.tsx`, `src/components/ExportStep.tsx`
- Source trace: `docs/native-chimaera-discovery-trace.md`
- Next-family supplied-source trace: `docs/native-geneconv-discovery-trace.md`
- Static contracts: `scripts/check-source-contract.mjs`, `scripts/check-source-balance.mjs`

The checkpoint remains source-only and must be treated as unvalidated until the authorized native
golden/runtime phase.
