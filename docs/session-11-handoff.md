# Session 11 handoff

## Checkpoint boundary

Session 11 advances the source checkpoint to `0.11.0-session-11` and project schema
`org.rdp-web.project/v1alpha11`. It adds the manual's automated query-vs-reference primary-analysis
scheme to the existing end-to-end RDP/MaxChi workflow:

dataset and role assignment → settings → constrained or exploratory cyclic scan → event
reconciliation → ordered review/repair → project, CSV, and alignment exports.

The implementation was derived only from the supplied manual and source, principally manual §4.2
and `Module3.bas::MakeAnalysisListQvR`. No alternate reference implementation was consulted.

Per the project instruction, no native source, C++, WebAssembly, TypeScript compiler, Vite build,
preview server, dependency install, or project/runtime test was executed. “Implemented” below means
source and cross-layer contracts exist; it does not mean compiled or native-parity validated.

## Implemented this session

- A first-class `AnalysisMode` in C++, the C ABI, worker, TypeScript state, results, persistence,
  and exports:
  - `exploratory` retains the existing all-distinct-origin triplet walk;
  - `query-reference` requires one query and two references from different positive group IDs; and
  - reference/query status constrains triplet composition, not recombinant/parent inference.
- Dataset role assignment:
  - blank or group `0` means query;
  - any positive 32-bit integer means reference group;
  - each row has an editable group control;
  - selected rows, including all matches beyond the 500-row render cap, can be assigned together;
  - group IDs can be compacted by first input appearance without changing membership;
  - “Detect REF names” maps the manual's `REF-A<name>`-style prefixes; and
  - “All queries” clears every assignment.
- Exact eligibility and plan feedback using the same enabled/unmasked and
  `max(5, RDP window)` usable-site threshold as the core. Settings block a constrained scan unless
  at least one query, two references, two groups, and one cross-group record triplet exist.
- A memory-bounded constrained scheduler:
  - active query and reference record vectors are `O(W)`;
  - a three-cursor reference-pair/query walk replaces the supplied growing `AnalysisList(2, T)`;
  - reference pairs remain the outer loops and queries the inner loop;
  - the exact scheduled count is computed without enumerating triplets; and
  - worker batches and cancellation remain bounded at 512 numerical triplets.
- Cyclic identity handling:
  - fragments inherit the reference group of their original sequence;
  - two working copies of one original cannot share a triplet;
  - active roles, groups, counts, correction, and cursors refresh after every tract erasure;
  - event-zero edits and all-rejected restore prefixes also force this refresh; and
  - a constrained round has its own `no-eligible-query-reference-triplets` terminal reason.
- Statistical accounting:
  - progress reports the exact cross-group reference-record/query workload;
  - the supplied `choose(reference groups, 2) × query` correction remains independently labelled;
  - the browser cyclic adaptation counts unique active query origins so fragments do not multiply
    the correction opportunity; and
  - both workload and correction products use saturating 64-bit arithmetic before the supplied cap.
- Live/review role context:
  - progress reports the current cyclic round's actual correction factor;
  - progress also reports current active working rows, query/reference records, and groups after
    tract erasure or fragment re-entry, including an honest zero-work terminal round;
  - all three current event-role cards and manual role choices expose query/reference input status;
  - a current reference-as-recombinant call is flagged in amber as in the manual; and
  - event JSON, CSV, and export summaries retain that distinction.
- Result/project/CSV coverage:
  - analysis scheme and query/reference summary;
  - full original-sequence reference-group vector;
  - triplet constraint and source correction-rule identifiers;
  - scheme-aware progress, null-result, cyclic-rescan, review, and export wording; and
  - query/reference counts and group count in every event CSV row.
- Backward restoration:
  - `v1alpha11` restores scheme and group IDs through a copied `uint32_t` worker buffer;
  - `v1alpha1`–`v1alpha10` remain fully exploratory; and
  - the pre-v10 MaxChi compatibility rule remains unchanged;
  - saved v11 signals are rejected if they repeat an origin or violate the constrained composition;
    and
  - absent legacy factors are rebuilt from the selected scheme, while C ABI role buffers must match
    the restored alignment exactly.
- MaxChi plot fidelity fixes:
  - reconstructed plots now pass the selected original triplet's supplied ten-character input-
    missing map to `make_banned_windows`; and
  - forced breakpoints, selected peak, pair maxima, and later-round reconstruction labels remain.

## Important fidelity and performance decisions

- The active supplied `MakeAnalysisListQvR` branch schedules every cross-group reference-record pair
  with every query. That ordering and constraint are authoritative for this phase.
- The source materializes `AnalysisList` in 10,000-column increments. The port generates the same
  logical initial catalogue lazily, preventing triplet-count-sized browser memory use.
- The source sets `MCCorrection = choose(RefNum, 2) × QNum`, which is not necessarily its materialized
  record-triplet count when groups contain multiple references. The port deliberately preserves
  those as separate concepts and exposes both.
- During cyclic fragment re-entry, the port uses unique query origins for that correction while the
  supplied routine refers to its current desktop record count. This avoids same-origin aliases
  manufacturing statistical opportunities but is explicitly marked for native golden comparison.
- The name parser is an editable convenience for the manual's documented prefix convention. It is
  not presented as a complete port of every `CheckQueryReference` auto-grouping heuristic.
- The desktop auto-assignment routine clears its mask array. The browser deliberately keeps
  query/reference assignment separate from visible curation and never silently unmasks rows;
  “Enable all” remains an explicit action.
- References are not assumed nonrecombinant by the detection or late-role path. This preserves the
  manual's warning that a reference can itself be identified as recombinant. Such current calls are
  highlighted in amber in event review and identified explicitly in CSV output.

## Static checks allowed and performed

The source-only checker verifies:

- every public `rdp_api.h` entry has a keepalive definition, CMake export, and worker contract/call;
- package, native engine, result, and `v1alpha11` schema versions agree;
- `v1alpha1`–`v1alpha11` import coverage and the pre-v11 exploratory boundary exist;
- query/reference options, scheduler markers, correction semantics, worker pointer lifetime,
  results/types, and all workflow UI surfaces are connected;
- the Session 10 MaxChi contracts remain present; and
- C++/TypeScript delimiter balance holds across the checked source set.

Shell, YAML, JavaScript-module syntax, and archive integrity are also checked without executing
project code. The Pages workflow remains configured for **Settings → Pages → Source: GitHub
Actions**; its first run will install and build remotely only when pushed to the default branch or
manually dispatched.

## Known open work

1. Establish supplied-desktop query/reference golden fixtures for pair/query order, grouped and
   ungrouped records, masks, empty groups, reference-as-recombinant events, correction factors,
   fragment re-entry, repairs, and project replay.
2. Compare the editable prefix mapping with native `CheckQueryReference` behavior and port only the
   source-supported heuristics that prove necessary; never silently alter explicit user groups.
3. Establish the MaxChi discovery/confirmation and BURT/BenHMM golden fixtures already specified in
   `docs/validation-plan.md`.
4. Trace and port supplied MaxChi permutation/manual-doublet modes as distinct workflows.
5. Add remaining detection-method families and post-group method-stack rechecks without hiding
   per-method evidence or overstating native weighting parity.
6. Run the first tiny, resource-limited compile/runtime checkpoint only after explicit permission.

## Resume map

- Analysis mode/options/scheduler: `wasm/src/rdp_method.hpp`, `wasm/src/rdp_method.cpp`
- C ABI and project schema: `wasm/include/rdp_api.h`, `wasm/src/rdp_api.cpp`
- Worker scan/restore bridge: `src/workers/analysis.worker.ts`
- Shared project/result types: `src/lib/types.ts`
- Dataset role editor and plan derivation: `src/App.tsx`, `src/components/DatasetStep.tsx`
- Scheme settings/progress/review/export: `src/components/SettingsStep.tsx`,
  `src/components/ScanStep.tsx`, `src/components/ReviewStep.tsx`,
  `src/components/ExportStep.tsx`
- Query/reference role context in lazy review payloads: `src/components/EventAlignmentInspector.tsx`,
  `src/components/EventTreeInspector.tsx`
- Styling: `src/styles.css`
- Static contracts: `scripts/check-source-contract.mjs`, `scripts/check-source-balance.mjs`
- Supplied workflow/fidelity mapping: `docs/workflow-mapping.md`, `docs/fidelity-notes.md`
- Golden corpus: `docs/validation-plan.md`

The checkpoint remains source-only and must be treated as unvalidated until the authorized native
golden/runtime phase.
