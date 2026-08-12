# RDP Web

RDP Web is a browser-native port of the Recombination Detection Program workflow. It is designed
as a static site: alignments are parsed and analysed locally in a Web Worker, while the numerical
core runs in WebAssembly.

> **Session 9 source checkpoint — intentionally uncompiled.** Per the project instruction, this
> archive contains source only. No native code, WebAssembly, npm build, preview server, or test
> executable was run while preparing it.

## What this checkpoint contains

- The RDP5 dataset → settings → primary scan → event reconciliation → ordered review → export workflow.
- FASTA, GDE, CLUSTAL/MUSCLE, sequential/interleaved PHYLIP, NEXUS, and MEGA alignment readers.
- Alignment diagnostics, pairwise identity checks, and the supplied RDP5 auto-mask workflow.
- The manual's distinct enabled/masked/disabled row states: masked sequences skip primary triplets
  but remain eligible for secondary evidence and tree placement; disabled sequences skip primary
  and secondary evidence while remaining available as phylogenetic context. Modern bulk controls
  restore the supplied auto-mask, enable all, mask all, or disable all without rendering every row.
- A C++20/WASM primary RDP scanner derived only from the supplied RDP sources:
  information-rich triplet sites, rolling pair counts, candidate tract boundaries, the binomial
  tail calculation, 169-site scaling, and RDP5's multiple-testing cap.
- Bounded 512-triplet worker batches, cancellation, reusable scan buffers, and plot downsampling.
- Strongest-first cyclic detection: scan all eligible triplets, reconcile the best event, infer its
  co-recombinant group, erase the event tract, then start a fresh full pass until no signal remains.
- Source-shaped fragment re-entry below the supplied 100,000-site cutoff. Erased tracts become
  gap-padded working sequences with original/event provenance; same-origin copies cannot share a
  triplet, short/duplicate fragments are omitted, and a visible 256-fragment cap bounds browser cost.
- Event hypotheses use the supplied detectable-signal rule: two shared original sequence identities
  and greater than 30% symmetric tract overlap, including fragment-assisted signals.
- Follow-up RDP profile checks of every masked (but not disabled) sequence, retaining trace-like
  profiles separately from statistically significant signals.
- Three role hypotheses per event, with each anchor sequence treated in turn as the presumed
  recombinant as required by the RDP5 reconciliation workflow.
- The supplied distance-pattern correlation layer: source-defined breakpoint flanks, five region
  profiles, three paired six-value Pearson tests, five inverse-category relabellings, dominant-pattern
  and distance-triangle warnings, `MakeGoodC` overlap eligibility, the active `MakeACOR` affinity
  gate, exact `MakeRList` dual-correlation override, positive aggregate scoring, `StripDupInv`, and
  the active first `FinalTrim` duplicate-correlation cleanup.
- The source `FinalTrim` `OKSeq 6` nearest-nonrecombinant fixed-point pass, including post-
  `StripDupInv` swap-last list order, exact `0.83`/`0.95`/`0.99` gates, direct-event availability,
  outside/inside tree bounds, and paired breakpoint-distance veto. Its retained lists feed both
  ascending final expansions and selected-role pruning.
- The complete active `FinalTrim` matrix-score family for `OKSeq` 7–14: collapsed/raw tree
  position, whole-tract and breakpoint JC distances, explicit source-zero slots 10/11, and
  `FindActualEvents`/`MakeMatchMatX2P` detected-region distance for 14. Warning gates, asymmetric
  ties, saturation, penalties, closest-pair modifiers, the active bare-`CompMat` index quirk, and
  the raw subtotal feed the active completed consensus score and remain auditable.
- The supplied `CalcMatchY` evidence path for `OKSeq` 17 and 18: four bounded 40-variable-site
  flank walks, VB half-to-even window rounding, signed match states, circular rolling smoothing,
  regional-product score, six breakpoint checkpoints, standard threshold class, and the opening
  `ConsensusOK` raw-tree topology-consistency filter. Raw/filtered classes and bounded fallback
  state remain auditable; available rows now drive active grouping.
- The active RFF=0 `FinalTrim` completion: nearest-nonrecombinant thresholds are retained for the
  correlation-gated expansion, raw whole-tree parent bounds drive the second expansion, and the
  selected-role branch preserves swap-last deletion plus the source's inherited third-list index.
- Completed `ConsensusOK` scoring and membership: `OKSeq` 0–6 and 15, `CheckPatternX`, `RCorrX`,
  the source's declared-`Long` `NS` narrowing, topology-filtered `OKSeq` 18, primary thresholds,
  exact-distance equivalence widening, straggler collection, and all-role empty fallback now rebuild
  each distance list before the manual's two-of-three group is formed.
- The shared selected-role conservative cleanup after the RFF guard: direct-distance outlier ranks,
  raw/direct topology constraints, swap-last removal, bounded-cluster admission, and the source's
  always-true `x = x` strict-inlier branch now finalize each distance list.
- The primary-RDP post-group recheck: every finalized nonrepresentative list candidate is rerun
  against its role's two representatives with the supplied `LowP * 100000` lift. Emitted signals,
  candidate-recombinant signals, event-overlapping traces, ordinary corrected significance, and the
  best tract/probabilities remain separately auditable.
- A source-shaped MaxChi confirmation kernel from the supplied `FastRecCheckMC2` path. It builds
  the three variable-site pair profiles, applies the DLL match-difference screen and `ChiPVal2P`
  approximation, excludes native `MissingData`/prior-erasure and linear-edge windows, grows the
  strongest peak, and retains raw, within-triplet, and project-corrected probabilities separately.
  The event representative triplet and every finalized nonrepresentative distance-list row are
  rechecked. The scan uses rolling match totals, so its strongest-peak pass is linear in variable
  sites rather than linear in both sites and window width.
- Six Jukes–Cantor neighbour-joining trees per event, deterministic ten-replicate column bootstrap,
  50% support collapse, and paired-tree phylogenetic-correlation membership.
- Iterative detectable-set closure and the manual's complete “present in at least two of three”
  co-recombinant group for every presumed-recombinant role.
- An auditable source-decision-tree role recommendation: nine displayed metrics cover direct, raw-
  tree and collapsed-tree PhylPro families, leave-one-role-out and displacement scores, weighted
  triplet ordering changes, and three-set context. Eight voting methods use the supplied full/half
  contribution weights; the unported native method families remain explicitly absent.
- Ordered event review with accept/reject decisions, editable roles, breakpoints, and complete
  co-recombinant membership. The current and automatic groups remain separate and auditable.
- An on-demand graphical breakpoint inspector derived from the original alignment. It prioritizes
  the recombinant, both parents, current/automatic co-groups, masked traces, and supporting
  evidence; shows bounded windows at both breakpoints; handles circular origin wrapping; and
  colour-codes parent matches without transferring the full alignment to the main thread. Persistent
  event context applies the manual's RDP-specific deleted-tract rule: a boundary within one RDP
  window, counted in information-rich triplet positions, is marked uncertain, while immediate
  tract contact remains distinguishable in JSON, CSV, and review UI. Its nearest-informative-state
  bracket is explicitly a review aid. The separate supplied BURT/BenHMM path now performs the
  three-state seeded HMM training and reports signed source-labelled 99%/95% ranges, HMM positions,
  breakpoint movement, missing/gap adjustments, and revert state without conflating either output.
  The supplied default-enabled “polish breakpoints” option is available in scan settings, retained
  in project checkpoints, and can be disabled to preserve the primary RDP coordinates.
- An on-demand graphical tree inspector for the same six event regions used by reconciliation.
  It compares whole-tract and both breakpoint pairs, labels bootstrap support and retained-fragment
  leaves, optionally collapses branches below 50%, and transfers compact saved edge lists rather
  than rebuilding trees or moving distance matrices to the interface.
- A batched erase/fragment/re-scan cycle that re-identifies later events after a correction or
  rejection. Core/API guards enforce review order; mid-repair project reload drops the stale tail,
  remaps retained signal anchors, and resumes at the changed event.
- Reloadable project JSON (`v1alpha9`, accepting `v1alpha1`–`v1alpha9`), expanded event-level CSV,
  full, enabled-only, and masked/disabled-only curation FASTA directly from the loaded dataset,
  accepted-group sequence removal,
  accepted-tract column removal, tract-masking FASTA, and event-ordered mosaic-fragment FASTA.
  Event-derived alignment readiness is enforced in both the interface and WASM boundary, and
  exports that would contain no sequences or no alignment columns fail clearly.
- Unsaved-analysis protection around the manual's save-often loop. Finished scans and subsequent
  review/repair changes visibly require a fresh project checkpoint; tab exit and destructive
  dataset/settings replacement warn until that checkpoint is downloaded.
- A GitHub Pages Actions workflow that installs locked JavaScript dependencies, provisions a pinned
  Emscripten toolchain, checks the C ABI/worker/version/schema contract and TypeScript, builds the
  compatibility WASM target and Vite site, validates the deployment artifact, and publishes it
  without a separate hosting service.
- A responsive, accessible React interface that requires no application server.

This is **not yet a parity-validated RDP5 replacement**. The primary RDP workflow now reaches an
end-to-end reviewed-event result and alignment variants, including the active late list-build and
selected-role cleanup paths.
MaxChi evidence is currently confirmation-only: primary RDP still discovers and coordinates every
event, and the native MaxChi multi-peak smoothing/destroy/retry scheduler is not yet active.
Unported role-method families, their remaining post-group method-stack rechecks, broader statistical
breakpoint probability/confidence diagnostics, and golden native comparison are still open. The
manual's RDP deleted-tract one-window uncertainty rule and the primary-RDP post-group
signal/probability recheck are active, as is the distinct BURT/BenHMM statistical 99%/95%
confidence and repositioning path; other detection and breakpoint-probability families remain later
milestones. Exports keep these fidelity boundaries explicit and do not label the implemented
native-weight subset as full desktop method parity.

See [STATUS.md](STATUS.md), [docs/fidelity-notes.md](docs/fidelity-notes.md), the
[late-consensus source trace](docs/native-late-consensus-trace.md), and the
[breakpoint-uncertainty source trace](docs/native-breakpoint-uncertainty-trace.md), the
[BURT/confidence source trace](docs/native-breakpoint-confidence-trace.md), and the
[MaxChi recheck source trace](docs/native-maxchi-recheck-trace.md), and the
[session-9 handoff](docs/session-9-handoff.md) before interpreting results or starting the next phase.

## Deploy with GitHub Pages Actions

The repository includes `.github/workflows/deploy-pages.yml`. To publish it:

1. Put the **contents of this `rdp-wasm` directory** at the root of a GitHub repository.
2. In GitHub, open **Settings → Pages** and set **Source** to **GitHub Actions**.
3. Push the workflow to the repository's default branch. The same workflow can also be started
   manually from the Actions tab.

Only a default-branch push deploys automatically; pushes to other branches create a skipped run.
The workflow publishes `dist/` only after the expected HTML, Emscripten loader, and valid WASM
binary are present and every built URL is safe for a Pages project subdirectory. The existing
relative Vite base also supports a user/organization root site or custom domain.

GitHub Pages does not provide the COOP/COEP headers required for `SharedArrayBuffer`, so this
workflow deliberately builds the single-worker compatibility module. Numerical work still stays
off the UI thread in the dedicated analysis worker.

## Build later (not run for this checkpoint)

Prerequisites:

- Node.js 20 or newer and npm
- Emscripten with `emcmake` on `PATH`
- CMake 3.20 or newer

```bash
npm ci
npm run build
```

The deployable static site will be written to `dist/`. Upload the **contents** of that directory to
the static host. Vite uses relative asset paths, so a project subdirectory is supported. The host
must serve `.wasm` files with the `application/wasm` media type.

For local UI development after a WASM build:

```bash
npm run dev
```

The default module uses one dedicated analysis worker and needs no special response headers. An
optional pthread-capable module can be produced with `npm run build:wasm:threads`; browsers will
only select it when the host supplies cross-origin isolation headers. Otherwise the single-worker
module remains the automatic fallback.

The generated files belong in `public/wasm/`:

- `rdp-core.mjs`
- `rdp-core.wasm`
- optional `rdp-core-threads.*` files

## Project layout

| Path | Purpose |
| --- | --- |
| `.github/workflows/` | Locked GitHub Pages build, artifact validation, and deployment |
| `src/` | React workflow, worker client, review plots, and exports |
| `src/workers/` | Isolated WASM bridge and bounded scan scheduler |
| `wasm/src/` | Alignment readers, primary RDP method, MaxChi recheck, BURT/BenHMM confidence, phylogenetics, reconciliation, trace checks, and exports |
| `wasm/include/` | Stable C ABI consumed by the worker |
| `scripts/` | Explicit WASM build entry point and Pages artifact verifier |
| `docs/` | Workflow, fidelity, architecture, and validation handoff |
| `package-lock.json` | Reproducible dependency graph used by `npm ci` and Actions |

## Data handling

No sequence is uploaded by the application. The selected file, encoded alignment, scan state, and
review decisions remain in the browser tab until the user explicitly downloads an export. The
project has no analytics, cookies, remote API, or persistence layer.

## Source basis

This checkpoint was produced from the supplied `DNA5.dll` C++ source, VB6 application source, and
RDP5 instruction manual. No alternate RDP implementation was consulted.
