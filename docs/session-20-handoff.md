# Session 20 handoff — supplied event-tree calculation path

Version: `0.20.0-session-20`
Project schema: `org.rdp-web.project/v1alpha18`

Session 20 replaces the generic event-tree implementation with the active calculation path traced
through the supplied RDP5 sources and manual. No alternate RDP implementation was consulted.

## Implemented in this checkpoint

- Traced `TestMoveInTree` through `MakeNJTreesP`, bundled Clearcut, `SEQBOOT2`,
  `FastBootDistIP6`, `TreeRepsP`, `TreeGroupsXP`, `CollapseNodes`, `Tree2ArrayP`, and
  `MakeTreeArrayXP2`.
- Replaced double-precision generic NJ with a bounded source-shaped single-precision Clearcut
  kernel, including ordered strict tie selection, the first-round terminal-`r2` initialization
  behavior, packed collapse order, final-clade distance behavior, and five-decimal branch
  serialization.
- Reproduced the Microsoft CRT `rand()` stream used by `SEQBOOT2`, including its two discarded
  values, retained unresampled replicate zero, site-major draw order, and shared configured seed.
- Reproduced the `FastBootDistIP6` missing-state and saturation rules without the old ten-site
  minimum.
- Reproduced the base-tree support pseudocount, VB6 `CLng` percentage rounding, and strict
  below-50% collapse decision.
- Replaced raw patristic analytical values with the source's `TreeMidP`/`UltraTreeDistP`
  midpoint-ultrametric, four-decimal, ascending rank-coded `rank / 1000` matrices. Weak nodes promote to
  existing parent ranks as in `CollapseNodesXP3`; the collapsed scale is not recompressed.
- Routed those matrices into the already-active TreePhPr, CollapsedTreePhPr, TreeSubPhPr,
  TreeSubDist, TrpScore, phylogenetic-correlation, role-consensus, final-trim, and recheck paths.
- Kept visual edge lengths separate from analytical ranks and exposed kernel, seed, rank-level,
  normalized-negative-branch, bootstrap, and support provenance in project/result and lazy tree
  JSON.
- Added a linked event-tree core regression and made it a required GitHub Actions/Pages gate.

The supplied manual and late VB source diverge at the breakpoint flanks: the manual specifies four
20-variable-site bootstrapped trees, while the late code builds 60-site direct-distance `DMatS`
flanks and reserves the bootstrapped `FAMat`/`SAMat` kernel for the whole outside/inside pair. This
checkpoint deliberately preserves the manual's advertised six-tree workflow and exposes the flank
target in JSON/UI. That choice is explicit rather than being labelled desktop-binary parity.

The browser intentionally repairs one unsafe writer/parser defect: the supplied writer can omit an
internal root-child length and the parser then consumes the next unrelated decimal token. The port
serializes every actual edge before applying the supplied four-digit clamp/midpoint/rank pipeline.
That deliberate numerical fix is covered in the source trace and remains a native-golden comparison
item.

## Performance behavior

Each of the six bounded regional families is computed once per reconciled event. A single compact
site-weight matrix supplies all ten bootstrap matrices, and the resulting evidence objects are
reused by every role metric. The existing 100-sequence tree panel cap and lazy edge-list inspector
remain in place, so no distance matrix or alignment slice crosses the worker boundary merely to
draw a tree.

## Validation completed for this checkpoint

- Strict host-linked compilation of all C++ translation units through the existing BootScan/public
  API regression.
- Dedicated event-tree core regression covering the seed-3 bootstrap stream, retained replicate,
  source support values, five-decimal serialization, rank matrices, determinism, and short-region
  valid-site rule.
- Source/ABI/version/schema checks, TypeScript checking, script syntax checks, and production web
  build.

The real Emscripten build and instantiated-WASM browser smoke suite remain GitHub Actions gates
because this scratch environment does not provide the configured Emscripten 5.0.1 toolchain.

## Fidelity boundary and next work

This is source-faithful implementation work, not a claim of native-parity validation. The next
useful validation phase is a small authorized native-vs-WASM corpus that records the six base
distance matrices, bootstrap node percentages, collapsed ranked matrices, role scores, and final
role for identical event inputs. Tree display/export breadth (alternative user-selected tree
methods and full general-purpose tree-editor workflows) remains distinct from the primary automatic
event-analysis path completed here.
