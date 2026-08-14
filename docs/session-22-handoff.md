# Session 22 handoff — supplied PHYLPRO event review

Version: `0.22.0-session-22`
Project schema: `org.rdp-web.project/v1alpha19`

Session 22 adds the supplied PHYLPRO left/right phylogenetic-profile calculation as a lazy,
review-only diagnostic. No alternate RDP implementation was consulted, and no supplied native
source was compiled.

## Added in this checkpoint

- Traced `FindSubSeqPP` → `PXoverD` → `MakePDstMat`/`UpdatePDstMat` → `PPRegression` in the supplied
  VB/DLL sources and mapped the manual's diagnostic-only workflow boundary.
- Added a C++20 PHYLPRO kernel with circular source windows, complete-window linear handling,
  pairwise-ignore or whole-column-strip missing policies, optional self inclusion, VB half-to-even
  window rounding, and the supplied Pearson/zero-variance behavior.
- Reduced ordinary plot work from an all-pairs matrix to the three target-to-context rows that the
  selected-event plot consumes. This preserves their coefficients while changing rolling work from
  `O(LN²)` to `O(LN)` and working distance memory from `O(N²)` to `O(N)`.
- Added a bounded on-demand C/WASM ABI, worker request, typed client, and Windows 95 review panel.
  The full profile is evaluated before retaining at most about 2,048 display samples plus ends,
  minima, and breakpoint-nearest points.
- Kept PHYLPRO outside discovery and reconciliation. The supplied RDP5 route has no implemented
  PHYLPRO significance test, so the browser reports no p-value and cannot create or reorder events.
- Confirmed that the supplied `Match(66/68/72/85)` assignments correctly reset RDP's encoded
  A/C/G/T counters each column; the browser's fresh state set is equivalent. Explicitly repaired
  the separate compact/original `RevSeq` direction defect in the supplied `FindSubSeqPP` snapshot;
  masked originals remain context, disabled rows are excluded, and cyclic fragments are not used.
- Added a deterministic host regression that checks every optimized point against brute-force
  recomputation, including a planted mosaic trough, both topologies, both gap modes, self policy,
  disabled context, and window caps. The linked public-API regression and production Pages WASM
  smoke test also call the lazy endpoint.
- Added PHYLPRO to source/ABI contracts and the GitHub Pages gate. The reloadable project schema
  stays `v1alpha19` because generated plot points are transient and reproducible from already-saved
  state.

## Fidelity status

The kernel is source-shaped active, host-regression-tested, and native-unvalidated. Circular
distance windows and the selected three curves follow the supplied intended route. Linear
complete-window behavior, explicit context indexing, and per-column state reset are documented
correctness adaptations. Explicit compact context mapping repairs the supplied `RevSeq` direction;
the manual's prose says PHYLPRO conceptually uses all sites, while the
active source maps polymorphic columns; the API and UI report the active-source mapping.
Browser-normalized ambiguity remains missing (and U remains T), while the literal DLL compares
non-gap raw character codes; ambiguity-heavy curves are consequently a golden-test boundary.

See [`native-phylpro-review-trace.md`](native-phylpro-review-trace.md) for the calculation,
optimization, compact-index repair, and validation boundary.

## Validation performed for this source checkpoint

- strict TypeScript checking;
- strict host C++20 PHYLPRO/brute-force regression;
- linked BootScan/public-ABI regression including the PHYLPRO endpoint;
- focused SISCAN and supplied event-tree regressions;
- source/ABI/version/schema/workflow contracts and shell syntax;
- production Vite build and source-only archive audit.

The Emscripten build and instantiated production-WASM Pages verifier remain authoritative in
GitHub Actions because this local toolchain does not provide Emscripten. Authorized desktop
PHYLPRO curves remain required for native golden comparison.

## Next high-value work

1. Capture authorized RDP5 PHYLPRO curves for gap/mask subsets, odd/even windows, circular edges,
   role ordering, self inclusion, flat rows, the encoded counter reset, and the documented
   `FindSubSeqPP` compact-index defect.
2. Capture authorized SISCAN, event-tree, BootScan, 3SEQ, GENECONV, MaxChi/CHIMAERA, and BURT golden
   fixtures before changing source-shaped labels to parity-validated.
3. Trace the remaining recombinant-identification/manual method families and complete only routes
   present in the supplied sources; keep external-program boundaries explicit.
