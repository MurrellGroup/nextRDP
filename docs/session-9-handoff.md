# Session 9 handoff

This checkpoint is source-only. No C++/Emscripten compilation, TypeScript checking, dependency
installation, Vite bundling, preview, server, browser runtime, or project execution was invoked.
The executable project check is limited to Node source-text contract/delimiter scanners; shell and
archive integrity checks do not execute the application.

## Reference provenance

Only the user-supplied reference set was used. No alternate RDP implementation or reference-code
site was consulted.

| Attachment | SHA-256 |
| --- | --- |
| `dnaDLLSource(3).zip` | `8a79a2adf7e733c57fb91cb0c50c5fd0e1183d528e7db45175efbcddc524b165` |
| `dna5DLLSource(3).zip` | `1fb01d49ed4513de9bd262aa943f82b16fcac3396cb449f2deb7751a10bf9102` |
| `VB Source(3).zip` | `6b1400378ac15c276a2030807c922e1826ddea557606e4fb718812acdf3e6933` |
| `RDP5Manual(3).pdf` | `ce9b4b0c88f3a93076e5d90502f855e247b0df6cc1bec82401e82dbfd3dc416c` |

The manual's load → inspect/curate → preliminary scan → hypothesis refinement → ordered review and
repair → save/export workflow remains the product contract. MaxChi is added inside hypothesis
refinement and cannot bypass or rewrite primary RDP discovery.

## Completed in this phase

1. Traced the supplied MaxChi call graph from VB `GetCriticalDiff` and list orchestration through
   DNA5 `FindSubSeqCP`, `MakeWindowSizeP`, `WinScoreCalcP`, `MakeBanWinP`, `CalcChiVals*`,
   `ChiPVal2P`, `GrowMChiWin*`, `FastRecCheckMC2`, and `AlistMC3`. The audit mapping is in
   `docs/native-maxchi-recheck-trace.md`.
2. Added `maxchi.cpp`/`maxchi.hpp` as an isolated C++20 confirmation kernel and included it in the
   WASM target. The kernel deliberately has no ABI endpoint and runs inside event reconciliation,
   keeping the existing worker boundary compact.
3. Ported the fixed-window profile gate: at least seven usable variable sites; 70-site default;
   source short-profile half-window arithmetic; `LowestProb / 6` critical χ² derivation; and strict
   absolute match-difference screening.
4. Ported the three MaxChi variable-site tracks (`s0=s1`, `s0=s2`, `s1=s2`) with original alignment
   coordinates. Compact-state zero remains the browser's existing gap/ambiguity exclusion.
5. Ported `MakeBanWinP`-shaped boundary handling, including distinct zero/length cells, both
   half-window ban directions, the trailing origin ban, native input missing runs, accumulated prior
   erasures, fragment gaps, and the additional linear-end exclusion.
6. Implemented the supplied 2×2 χ² statistic and source `NormalZ`/`ChiPVal2P` polynomial including
   its extreme-tail fallback. Lookup-table rounding remains an explicit golden-test boundary.
7. Replaced the first draft's repeated half-window sums with six rolling totals across the three
   pair tracks. After initial totals, each boundary costs constant work, so the strongest-peak scan
   is linear in variable-site count rather than variable-site count times window width.
8. Added a source-shaped strongest-window growth pass: `H/4` start with six-site floor, source
   `MaxX = 0 → 1` entry behavior, paired expansion, non-decreasing χ² retention, failure reset/count,
   and evaluation-before-`MDMap` stopping. The browser conservatively stops before the two
   half-windows overlap; this is documented rather than hidden.
9. Kept raw `ChiPVal2P`, within-triplet `V / min(initial,grown) × 3`, and project-corrected
   probabilities as separate values. Each evidence record also retains whether Bonferroni was
   applied and the exact event-round triplet opportunity count used, avoiding a later-round count
   being displayed for an earlier event.
10. Added event-level corroboration for the actual working representatives. Fragment-assisted
    events select their retained working fragment rows before MaxChi runs, while result identities
    stay mapped to the original recombinant and parents.
11. Added finalized-list corroboration beside the primary-RDP post-group recheck. Each surviving
    nonrepresentative candidate is compared with its role's two working representatives; the own
    representative follows the supplied loop's skip and is covered by the event-level card.
12. Added capacity-retaining MaxChi workspace buffers for alignment-coordinate prefixes, all three
    match tracks, ban/`MDMap` arrays, and the triplet missing mask. Repeated finalized-list rechecks
    reuse allocations rather than constructing alignment-sized vectors per candidate.
13. Extended result/project JSON with `maxChiTripletRecheck` and `postGroupMaxChiRecheck`, plus
    machine-readable late-consensus fields that distinguish strongest-peak confirmation from full
    native MaxChi discovery. Project schema advanced to `v1alpha9`; imports accept `v1alpha1` through
    `v1alpha9`.
14. Extended event CSV with MaxChi status, variable sites, windows, pair, peak coordinate, χ², raw
    and within-triplet probabilities, event-round correction count, corrected probability, and hit.
    Per-list evidence remains in the late-consensus diagnostic cell and project JSON.
15. Added a modern responsive event card and finalized-list column. The UI shows availability,
    missing/edge filtering, strongest pair/peak, all probability scopes, and the explicit statement
    that MaxChi does not discover or reposition an event in this checkpoint.
16. Updated the settings method panel to say “RDP full · MaxChi recheck”: RDP is included in event
    discovery, MaxChi is included in confirmation, and every other method stays visibly queued.
17. Advanced package/native/result labelling to `0.9.0-session-9`, updated source-contract coverage
    for MaxChi/CMake/schema/UI fields, and refreshed README, status, architecture, workflow,
    fidelity, validation, changelog, and late-consensus documents.
18. Preserved the existing GitHub Pages Actions deployment. With the project contents at repository
    root and **Settings → Pages → Source: GitHub Actions**, a default-branch push or manual dispatch
    performs the locked Node 20/Emscripten 5.0.1 build, verifies the static artifact, and deploys it.

## Deliberate fidelity boundary

Session 9 does **not** claim a complete MaxChi port. Primary RDP remains the sole discovery and
coordinate authority. The following supplied paths are still pending:

- `SmoothChiValsP` and ordered smoothed-peak selection;
- `DestroyPeakP`, the up-to-100 retry loop, and `WasteOfTime` termination;
- MaxChi breakpoint-pair/event construction and persistent event-catalogue writes;
- exploratory `AlistMC3` scheduling across all eligible triplets;
- full native method-origin dispatch in the late recheck stack;
- native lookup-table/float rounding and desktop-vs-browser golden validation.

Result metadata keeps `eventDiscoveryApplied: false`,
`maxChiKernelStatus: "source-shaped-strongest-peak-unvalidated"`, and
`nativeMaxChiFullRecheckComplete: false`. The UI and export caveat repeat that boundary.

## Static inspection performed

- Rechecked supplied attachment hashes and inspected only the supplied DLL/VB/manual sources.
- Ran `npm run check:source`: 36 public ABI functions align across the C header, C++ keepalive
  definitions, CMake exports, and worker calls; package/native/result version and `v1alpha9` schema
  align; MaxChi evidence/CMake/UI contracts are present.
- Scanned delimiters across 29 C++/header/TypeScript/TSX source files.
- Checked `scripts/build-wasm.sh` with `bash -n` and reviewed the GitHub Pages workflow/relative Vite
  base without invoking either build path.
- Rechecked the workflow's action majors against their official repositories/Marketplace on
  2026-08-12: checkout 7, setup-node 7, configure-pages 6, upload-pages-artifact 5,
  deploy-pages 5, and setup-emsdk 16 are published lines.
- Verified the delivered ZIP with `unzip -t` after packaging.

## Recommended next fidelity phase

1. Port MaxChi's smoothing, ordered peak selection/destruction, retry bounds, and event construction
   from the supplied `FastRecCheckMC2` path while retaining the current rolling scan and explicit
   method stream.
2. Add source-derived discovery scheduling without allowing a partial MaxChi call to become the
   RDP coordinate authority or bypass strongest-first cyclic event reconciliation.
3. Build the focused MaxChi golden corpus in `docs/validation-plan.md`: zero/length bans,
   `MaxX = 0`, strict critical differences/ties, linear ends, missing/erasure union, growth, extreme
   `ChiPVal2P`, round correction counts, and multi-peak retry order.
4. Continue with GENECONV and the remaining method families only from the supplied sources, keeping
   per-method evidence and native full/half role contributions explicit.
5. If execution is later authorized, start with a disposable resource-limited build and the
   smallest three-sequence fixtures; do not begin with a production alignment.
