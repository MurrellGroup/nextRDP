# RDP Web

RDP Web is a browser-native port of the Recombination Detection Program workflow. It is designed
as a static site: alignments are parsed and analysed locally in a Web Worker, while the numerical
core runs in WebAssembly.

> **Session 5 source checkpoint — intentionally uncompiled.** Per the project instruction, this
> archive contains source only. No native code, WebAssembly, npm build, preview server, or test
> executable was run while preparing it.

## What this checkpoint contains

- The RDP5 dataset → settings → primary scan → event reconciliation → ordered review → export workflow.
- FASTA, GDE, CLUSTAL/MUSCLE, sequential/interleaved PHYLIP, NEXUS, and MEGA alignment readers.
- Alignment diagnostics, pairwise identity checks, and the supplied RDP5 auto-mask workflow.
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
- Follow-up RDP profile checks of every masked sequence, retaining trace-like profiles separately
  from statistically significant signals.
- Three role hypotheses per event, with each anchor sequence treated in turn as the presumed
  recombinant as required by the RDP5 reconciliation workflow.
- The supplied distance-pattern correlation layer: source-defined breakpoint flanks, five region
  profiles, three paired six-value Pearson tests, five inverse-category relabellings, dominant-pattern
  and distance-triangle warnings, `MakeGoodC` overlap eligibility, the active `MakeACOR` affinity
  gate, exact `MakeRList` dual-correlation override, positive aggregate scoring, `StripDupInv`, and
  the active first `FinalTrim` duplicate-correlation cleanup as non-pruning diagnostics.
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
- A batched erase/fragment/re-scan cycle that re-identifies later events after a correction or
  rejection. Core/API guards enforce review order; mid-repair project reload drops the stale tail,
  remaps retained signal anchors, and resumes at the changed event.
- Reloadable project JSON (`v1alpha6`, accepting `v1alpha1`–`v1alpha6`), expanded event-level CSV,
  accepted-tract masking FASTA, and event-ordered mosaic-fragment FASTA. Final alignment readiness
  is enforced in both the interface and WASM boundary.
- A responsive, accessible React interface that requires no application server.

This is **not yet a parity-validated RDP5 replacement**. The primary RDP workflow now reaches an
end-to-end reviewed-event result and alignment variants, but unported role-method families, the
remaining late `ConsensusOK`/`FinalTrim` score-and-prune stack, breakpoint uncertainty beside
previously deleted tracts, and golden native comparison are still open. The other detection methods
and graphical alignment/tree review panels remain later milestones. Exports keep these fidelity
boundaries explicit and do not label the implemented native-weight subset as the full consensus.

See [STATUS.md](STATUS.md), [docs/fidelity-notes.md](docs/fidelity-notes.md), and the
[session-5 handoff](docs/session-5-handoff.md) before interpreting results or starting the next phase.

## Build later (not run for this checkpoint)

Prerequisites:

- Node.js 20 or newer and npm
- Emscripten with `emcmake` on `PATH`
- CMake 3.20 or newer

```bash
npm install
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
| `src/` | React workflow, worker client, review plots, and exports |
| `src/workers/` | Isolated WASM bridge and bounded scan scheduler |
| `wasm/src/` | Alignment readers, primary RDP method, phylogenetics, reconciliation, trace checks, and exports |
| `wasm/include/` | Stable C ABI consumed by the worker |
| `scripts/` | Explicit WASM build entry point |
| `docs/` | Workflow, fidelity, architecture, and validation handoff |

## Data handling

No sequence is uploaded by the application. The selected file, encoded alignment, scan state, and
review decisions remain in the browser tab until the user explicitly downloads an export. The
project has no analytics, cookies, remote API, or persistence layer.

## Source basis

This checkpoint was produced from the supplied `DNA5.dll` C++ source, VB6 application source, and
RDP5 instruction manual. No alternate RDP implementation was consulted.
