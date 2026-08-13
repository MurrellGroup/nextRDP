# Supplied-source trace — event-tree calculation kernel

This note records the Session 20 trace used to replace the generic event-tree helper with the
calculation path actually reached by RDP5's event reconciliation. Only the supplied VB, `DNA.dll`,
`DNA5.dll`, and manual sources were inspected. No alternate RDP implementation was consulted.

## Active call path

The supplied desktop tree kernel enters `Module3.bas::TestMoveInTree`. For the whole inside/outside
event partition it:

1. prepares the selected event-tree sequences and the inside/outside or breakpoint-flank
   subalignment;
2. calls `MakeNJTreesP`, which builds the unbootstrapped distance/tree matrices and invokes the
   bundled Clearcut implementation in traditional neighbour-joining mode;
3. calls `SEQBOOT2` to create a site-major weight matrix with the original data in replicate zero
   and resampled replicates in columns 1 through `Reps`;
4. calls `FastBootDistIP6(1, Reps, ...)` to create single-precision Jukes–Cantor matrices for the
   resampled columns;
5. calls `TreeRepsP`, which runs Clearcut for each bootstrap matrix and matches its groups against
   the base-tree groups;
6. converts every node count with
   `CLng((DLen + 1) / (Reps + 1) * 100)`, thereby retaining the base tree as a pseudocount;
7. calls `CollapseNodes` with a cutoff of 50; and
8. uses the raw and collapsed topology-distance matrices in the later phylogenetic-correlation and
   recombinant-role decisions.

The manual describes six event subalignments: four short flanks around the two breakpoints and the
whole inside/outside regions. It describes Jukes–Cantor distance, bootstrapped NJ, and collapse of
branches below 50%, which agrees with this active source path.

There is a manual/source divergence worth making explicit. The later supplied VB path builds the
whole inside/outside `FAMat`/`SAMat` through the bootstrapped tree kernel, but builds its four
breakpoint `DMatS` inputs from `MakeBPosLR(VSN=60)` plus direct distances. The manual instead says
that all six subalignments receive bootstrapped NJ trees and that each breakpoint flank spans 20
variable positions. RDP Web follows that documented six-tree workflow and applies the traced kernel
to every family, using a separate 20-informative-site boundary pass; it does not silently reuse the
60-site distance-correlation flanks. Authorized golden comparison should therefore report both the
current desktop matrices and the manual-intended six-tree result.

## Observable numerical details retained

### `SEQBOOT2`

`SEQBOOT2` calls `srand(seed)`, discards two `rand()` values before the draw loop, and then iterates
sites outside replicates. The supplied Windows build therefore exposes the Microsoft CRT 15-bit
`rand()` stream. Session 20 implements that LCG directly. It keeps replicate zero at weight one for
every site and uses the same configured `BSRndNumSeed` for each regional call, as the desktop does.

### `FastBootDistIP6`

Only coordinates where both states are non-missing contribute. A pair with at least one comparable
site uses single-precision Jukes–Cantor correction; identity at or below 0.25, or no comparable
sites, produces the source saturation distance 10. The earlier browser helper incorrectly imposed
a ten-comparable-site floor on this event-tree path.

### Bundled Clearcut

The port retains the observable behavior of the supplied Clearcut copy rather than substituting a
generic double-precision NJ routine:

- matrices and row transforms are single precision;
- the minimum transformed-distance scan is ordered and uses strict `<`, retaining the first tie;
- the final initial `r2` slot remains zero-initialized on the first join, matching `NJ_init_r`;
- the packed-matrix collapse order follows `NJ_decompose`/`NJ_collapse`;
- the final two clades each receive the remaining distance from the `NJ_LAST` path; and
- the Newick path suppresses a negative branch sign and truncates the magnitude to five decimals
  before `Tree2ArrayP` reads it back.

The graph endpoint exposes those serialized branch lengths. They are useful for a faithful visual
tree, but they are not the matrix values used by RDP's analytical role decisions.

### Midpoint-ultrametric ranked analytical distances

`Tree2ArrayP` parses branches at four fractional digits, clamps them to `[0, 1]`, and first builds
leaf-to-leaf path lengths. `TreeMidP` selects the first maximum-distance pair, locates its midpoint,
and `UltraTreeDistP` lengthens each terminal path so all leaves are equally distant from that root.
`MakeTreeArrayXP2` rounds the resulting ultrametric distances to four decimals, sorts distinct
levels ascending, and replaces each distance with `rank / 1000`. `CollapseNodesXP3` does not rerank
the matrix: it promotes a weak node's pairs to an existing parent-node rank, retaining any gaps in
the original rank scale. Session 20 now supplies those raw and collapsed rank matrices to
`TreePhPr`, `CollapsedTreePhPr`,
`TreeSubPhPr`, `TreeSubDist`, triplet-tree scoring, phylogenetic-correlation set construction, and
the existing downstream consensus/recheck stages.

The supplied Newick writer omits a length after one internal root child, while `TreeToArrayP` scans
forward for the next decimal without checking that it belongs to that node. Reproducing that
out-of-token read would make an internal edge depend on the first numeric token in the opposite
subtree. Session 20 deliberately repairs this handoff: every actual Clearcut edge is serialized,
parsed at the source's four-digit precision, clamped, and then passed through the same mathematical
midpoint/ultrametric/rank transformation. This can differ from a desktop result that triggers the
parser defect; it is an explicit numerical-correctness fix, not hidden parity drift.

## Deliberate browser bounds

The existing event-tree panel cap remains 100 working sequences and the event workflow uses ten
bootstrap replicates. The six matrices are built once per reconciled event and reused across all
downstream role metrics. The lazy graphical endpoint transfers only the already-saved bounded edge
lists; it does not rebuild a tree or copy alignment rows.

## Validation boundary

The linked host regression fixes the seed-3 `SEQBOOT2` weight stream, matrix dimensions, rooted
Clearcut-shaped node/edge counts, five-decimal branches, the allowed pseudocount support
percentages, 50% collapse metadata, rank/1000 analytical matrices, determinism, and the
one-comparable-site rule. It does not claim native-vs-WASM golden parity: exact supplied-binary
fixtures remain an authorized validation phase.
