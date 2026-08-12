# RDP Web

RDP Web is a browser-native port of the Recombination Detection Program workflow. It is designed
as a static site: alignments are parsed and analysed locally in a Web Worker, while the numerical
core runs in WebAssembly.

> **Session 17 source checkpoint — cyclic shortlist and probability-scope build verified.** This
> archive contains source only. Its Node 20 / Emscripten 5.0.1 production build,
> TypeScript/source contracts, Pages artifact, real WASM FASTA-upload test, and deterministic
> two-event cyclic-shortlist test passed before packaging.

## What this checkpoint contains

- The RDP5 dataset → settings → primary scan → event reconciliation → ordered review → export workflow.
- FASTA, GDE, CLUSTAL/MUSCLE, sequential/interleaved PHYLIP, NEXUS, and MEGA alignment readers.
- Alignment diagnostics, pairwise identity checks, and the supplied RDP5 auto-mask workflow.
- The manual's distinct enabled/masked/disabled row states: masked sequences skip primary triplets
  but remain eligible for secondary evidence and tree placement; disabled sequences skip primary
  and secondary evidence while remaining available as phylogenetic context. Modern bulk controls
  restore the supplied auto-mask, enable all, mask all, or disable all without rendering every row.
- The manual's automated query-vs-reference workflow. A dataset-level role editor accepts explicit
  numeric reference groups or detects documented `REF-A<name>`-style prefixes. Filter-aware
  multi-row selection can assign or clear a group without rendering every selected row, and group
  IDs can be compacted by first input appearance without changing membership. The primary scheduler
  then lazily emits exactly one query plus two references from different groups. Reference records
  remain eligible to be inferred as recombinants, as the manual requires. The exact record-triplet
  workload and the supplied `reference-group pairs × queries` correction factor remain separate and
  visible, avoiding an `O(T)` materialized analysis list on large inputs. Reference-as-recombinant
  calls receive the manual's distinct amber review treatment and explicit JSON/CSV input-role data.
- A C++20/WASM primary RDP scanner derived only from the supplied RDP sources:
  information-rich triplet sites, rolling pair counts, candidate tract boundaries, the binomial
  tail calculation, 169-site scaling, and RDP5's multiple-testing cap.
- A second, independently auditable C++20/WASM discovery stream for the supplied MaxChi
  `MCXoverF` workflow: all three variable-site pair profiles, native critical-difference screening,
  raw strongest-peak order, source-shaped window growth, `FindSide`, `OptLeftBPMC`/
  `OptRightBPMC`, tract construction, literal twelve-term/eleven-divisor `SmoothChiValsP`, completed/rejected peak
  destruction, and the supplied three-wasted/100-attempt retry bounds. A linear-time heap build plus
  lazy destroyed-peak rejection preserves the source's chi/boundary/pair order without repeatedly
  rescanning every surviving profile cell.
- A third independently labelled discovery stream for the supplied CHIMAERA workflow. Every triplet
  member rotates through the candidate-recombinant role; monomorphic/all-different sites are removed;
  the target is encoded as a binary match to either parent; and raw χ² peaks enter the same
  source-shaped growth, tract-side, boundary-optimization, destruction, and retry lifecycle. The
  supplied default is 60 information-rich sites. RDP categories, MaxChi equality tracks, CHIMAERA
  target inputs, and the triplet missing/erasure map are prepared in one alignment-byte pass.
- A fourth independently labelled discovery stream for the supplied ordinary automated GENECONV
  workflow. The shared non-monomorphic profile becomes three inner pair-match and three outer
  discordant-sequence tracks; source mismatch penalties, lambda/K calculation, strict critical
  scores, Karlin–Altschul tails, six-track provisional roles, stable lowest-P selection, and the
  configured overlap rule are retained. The active automated defaults are ignored indels, `G=1`,
  one overlap, and inactive minimum-fragment controls. Prefix/next-lower/rightmost-maximum queries
  preserve the start-extension stop/tie rules in `O(R log R)` per track, and lazy range coverage
  keeps fragment selection bounded without another alignment scan.
- A fifth independently labelled discovery stream for the supplied automated 3SEQ workflow. Each
  triplet member rotates through the candidate-recombinant role in source order; only sites where
  the parents differ and the target matches exactly one are retained; and the source maximum
  descent/ascent plus `CheckwrapC` tract construction is preserved. A compact exact hypergeometric
  random-walk dynamic program replaces the desktop four-dimensional lookup table within a bounded
  transition budget, with the supplied `SiegmundDiscrete` and scaled-exact routes explicitly
  labelled beyond it. The pre-`CheckwrapC` probability excursion is retained separately from any
  larger origin-extended boundary excursion, matching the supplied call order. 3SEQ applies the
  source's `p > 10^-15` Dunn–Šidák / smaller-tail product branch and adds no alignment-byte pass
  because it reuses the MaxChi/CHIMAERA equality profile.
- Bounded 512-triplet worker batches, cancellation, reusable scan buffers, and method-aware plot
  downsampling that forces both breakpoints, the selected method peak, and applicable profile maxima
  into the browser payload. CHIMAERA displays only its selected target/parent-one trace; GENECONV
  displays a three-colour `-log10(raw KA P)` inner/outer fragment envelope; 3SEQ displays all three
  signed target-specific random walks on a zero-aware axis. Later-round
  reconstructions are labelled when the exact historical erased/fragment profile was not serialized.
- Combined strongest-first cyclic detection: scan all eligible triplets in the supplied
  RDP/GENECONV/MaxChi/CHIMAERA/3SEQ method-major order, reconcile the best event, infer its co-recombinant group, erase the event
  tract, then retain XOverList-style summaries for unchanged exact working triplets. Triplets that
  touch erased rows or new fragments run fresh kernels; unchanged triplets replay signals and skip
  stable method work. 3SEQ refreshes once when its post-erasure split mode first activates.
- The source `MakeMCCorrection` factor is fixed from the initial scan plan while each later round's
  actual fragment-expanded workload remains visible separately. GENECONV review distinguishes raw
  `GCCalcPValP2` probability from the RDP5 `XOverList`-equivalent project-corrected value and makes
  clear that later BURT polishing does not recalculate either value.
- Event review and CSV export explicitly flag MaxChi+CHIMAERA-only support because the manual treats
  those methods as closely related rather than independent confirmation.
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
- A separate source-shaped MaxChi confirmation kernel from the supplied `FastRecCheckMC2` path. It builds
  the three variable-site pair profiles, applies the DLL match-difference screen and `ChiPVal2P`
  approximation, excludes native `MissingData`/prior-erasure and linear-edge windows, grows the
  strongest peak, and retains raw, within-triplet, and project-corrected probabilities separately.
  The event representative triplet and every finalized nonrepresentative distance-list row are
  rechecked. The scan uses rolling match totals, so its strongest-peak pass is linear in variable
  sites rather than linear in both sites and window width.
- Source-shaped CHIMAERA, GENECONV, and 3SEQ late corroboration over those same representative and
  finalized-list triplets. `FastRecCheckChim` rotates all three targets, ordinary `GCXoverD` retains
  its best six-track KA fragment, and `TSXOver(1)` evaluates both split walk orientations plus the
  supplied inverse-parent/inverse-interval list copy. These records keep their own probability
  scopes and never move a reconciled event.
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
- Reloadable project JSON (`v1alpha16`, accepting `v1alpha1`–`v1alpha16`), expanded event-level CSV,
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

This is **not yet a parity-validated RDP5 replacement**. The fully exploratory and automated
query-vs-reference RDP/MaxChi/CHIMAERA/GENECONV/3SEQ workflows now reach
an end-to-end reviewed-event result and alignment variants, including the active late list-build and
selected-role cleanup paths. MaxChi and CHIMAERA exploratory discovery are active in source, but
their indexing, smoothing/destruction basins, preliminary role assignment, and cross-method event
order still require native golden comparison before either can be called parity validated. Native
MaxChi manual-doublet/permutation, CHIMAERA permutation/full late-event-reconstruction, and
GENECONV permutation/manual-pair/alternative-indel/full late-event-reconstruction modes are not yet
represented. The ordinary ignored-indel KA GENECONV path is source-shaped active but unvalidated.
The ordinary 3SEQ exact/random-walk path is likewise source-shaped active but unvalidated. Later
cyclic rounds now retain the supplied `FindSubSeqTS2` inclusive position map and
`CheckSplit3Seq`/`SubPVal` missing-run trim, reverse-orientation retry, and corrected-P re-gate.
The supplied two-orientation `TSXOver(1)` representative/finalized-list recheck is also active;
manual permutation envelopes and full late event-catalogue reconstruction remain explicit
boundaries. Source-shaped FastRecCheckChim strongest-target, ordinary-kernel GENECONV, and 3SEQ
Findall rechecks remain non-coordinate-changing and require native golden comparison.
Unported role-method families and their post-group method-stack rechecks, broader statistical
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
[MaxChi discovery trace](docs/native-maxchi-discovery-trace.md), the
[MaxChi recheck source trace](docs/native-maxchi-recheck-trace.md), and the
[CHIMAERA discovery trace](docs/native-chimaera-discovery-trace.md), the
[GENECONV discovery trace](docs/native-geneconv-discovery-trace.md), the
[3SEQ discovery trace](docs/native-threeseq-discovery-trace.md), the
[cyclic-shortlist trace](docs/native-cyclic-shortlist-trace.md), and the
[Session 17 handoff](docs/session-17-handoff.md) before interpreting results or starting the next phase.

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

## Automated query vs reference workflow

1. Load the aligned dataset. Names using the manual's `REF-A<sequence name>` convention are
   detected into editable reference groups; “Detect REF names” can reapply that mapping.
2. In the dataset table, leave a role blank for a query or enter a positive integer for its
   reference group. Select individual rows—or every row matching the current name filter—to assign
   a group or make them queries in one action. “Compact groups” renumbers IDs by first appearance
   without changing group membership. Records in the same reference group are not paired with one
   another. Role assignment does not silently change enabled, masked, or disabled curation.
3. On the settings page, choose **Query vs reference**. The plan must contain at least one eligible
   query and eligible references from two groups. The screen shows the exact scheduled record
   triplets and the separately capped group-pair × query correction factor.
4. Run and review in the same strongest-first cyclic order as an exploratory analysis. The role
   inference remains free; a reference called as recombinant is highlighted in amber and keeps its
   input group on the role cards, alignment view, tree view, JSON, and CSV.
5. Download the `.rdpweb.json` checkpoint to preserve all assignments and review decisions. The
   CSV is a readable event summary, not a substitute for that complete project state.

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
| `wasm/src/` | Alignment readers, RDP/MaxChi/CHIMAERA/GENECONV/3SEQ discovery and active rechecks, BURT/BenHMM confidence, phylogenetics, reconciliation, trace checks, and exports |
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
