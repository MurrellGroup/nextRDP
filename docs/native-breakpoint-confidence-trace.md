# Native BURT / breakpoint-confidence trace

This is the source-to-port map for the statistical breakpoint-polishing path in the user-supplied
RDP5 VB and DNA5 DLL sources. No alternate implementation was consulted. The active Session 8 port
lives in `wasm/src/burt_confidence.cpp`; current parent-state brackets remain review aids, while
BURT ranges are emitted separately as signed statistical evidence.

## Product role

The manual's BURT section describes a windowless hidden Markov model over the information-rich
three-sequence alphabet. The supplied active source invokes it after the preliminary triplet and
role refinement, before storing the event's ten breakpoint-confidence values. It may move the two
reported breakpoints as well as produce 99% and 95% ranges, so it cannot be safely approximated by
decorating the existing breakpoints with a generic interval.

The supplied application initializes `PolishBPFlag = 1` (`MainForm22.frm` 29988), so the browser
setting is also enabled by default. The settings screen exposes the flag because disabling it is a
legitimate desktop workflow choice; the selected value crosses the C ABI, is retained in project
JSON, and is restored with old projects defaulting to enabled.

The manual describes the conceptual method (including a step-up state count and ten initial
conditions). The supplied active source is the parity authority for this port: its `BenHMM` call
uses three hidden states, passes `HMMCycles = 20`, and the DLL loops from zero through that value
inclusive (21 starts).

## Active control flow

| Supplied source | Lines | Active behavior to preserve |
| --- | ---: | --- |
| `Module3.bas` event processing | 43874–43894 | Optionally improves the triplet, then calls `PolishBP(20, 0, ...)` when breakpoint polishing is enabled and the alignment is not compressed |
| `Module3.bas` event storage | 44127–44133 | Copies all ten `CIOut` values into `BPCIs(:, SEventNumber)` |
| `Module2.bas` role replacement | 17549–17558 | Re-runs `PolishBP(20, 0, ...)` if the selected representative triplet changes; the comments deliberately avoid updating the already-built distance matrices |
| `Module4.bas` `PolishBP` entry | 9615–9707 | Saves input roles/breakpoints, handles reassortment-only segment boundaries, computes event span, and calls the single/double `BenHMM` path |
| `Module4.bas` interval selection | 9711–9744 | Calls `MatchBPtoCI` independently for beginning and ending and writes the 99%, HMM position, and 95% values into the ten-element output |
| `Module4.bas` repositioning | 9750–10970 | Applies sign-based match status, circular distance/event-span gates, missing-data scans, segment-boundary logic, a minimum three-informative-site guard, and final gap relocation |
| `Module4.bas` `MatchBPtoCI` | 10971–11074 | Normalizes expanded circular coordinates, tests interval containment, measures distance in `XPosDiff` space, and chooses the nearest containing interval before the nearest non-containing interval |
| `Module4.bas` `BenHMM` | 11076–11866 | Sorts triplet identities, recodes information-rich sites, expands circular data, trains the HMM, computes posteriors, and constructs all candidate breakpoint intervals |
| `MathFuncsDll.cpp` `DoHMMCyclesSerial` | 25350–25578 | Runs randomized Viterbi training and retains the highest-likelihood transition/emission matrices and lattice path |
| `MathFuncsDll.cpp` helpers | 28047–28643 | Supplies Viterbi, forward, reverse, count-update, and lattice backtracking primitives used by the active path |

## Information-rich alphabet and circular expansion

`BenHMM` first sorts the three input sequence indices. For every alignment coordinate at which all
three supplied `SeqNum` values are not ASCII 46, it emits exactly one symbol when two sequences
agree and the third differs:

| Pair that agrees | Source symbol |
| --- | ---: |
| sorted sequence 1 + 2 | `0` |
| sorted sequence 2 + 3 | `1` |
| sorted sequence 1 + 3 | `2` |

`XDiffPos` maps each emitted symbol back to the alignment; `XPosDiff` carries the cumulative emitted
position at every alignment coordinate. This is not identical to treating every polymorphic column
as informative.

For circular data the source rotates by half of the information-rich length and adds overlap on
both sides. After posterior calculation it selects the centered half again. The unusual array
offsets and endpoint assignments at `Module4.bas` 11151–11188 and 11527–11536 are parity-sensitive;
they must be translated from VB's inclusive arrays explicitly rather than normalized into a
generic modulo view.

## Training model

The active call supplies three symbols and three hidden states. `DoHMMCyclesSerial` preserves these
source choices:

- initial self-transition probability `1 - 5 / alignmentLength`; the remaining mass is divided
  between the other states (the DLL stores and subtracts that mass as `float` before taking logs);
- uniform initial hidden-state probabilities;
- randomized, imbalanced emission seeds generated from the supplied global seed;
- cycles `0 ... HMMCycles` inclusive;
- at most 100 Viterbi-training iterations per start;
- `0.01` pseudocounts for every transition and emission cell;
- convergence when the lattice-path likelihood exactly equals the preceding value;
- retention of only the highest-likelihood start's matrices and lattice path.

The browser port must preserve the DLL's flattened matrix order and random-number behavior if exact
desktop reproduction is required. Replacing this with a third-party HMM, EM implementation, or a
deterministic heuristic would not be a faithful port.

## Posterior ranges

After the best trained matrices are restored, `BenHMM` runs log-space forward and reverse passes,
adds their values by state, subtracts the per-site maximum, exponentiates, and normalizes. A lattice
state transition between positions `x` and `x + 1` creates a candidate breakpoint centered between
the corresponding `XDiffPos` coordinates. VB narrows the complete midpoint expression to `Long`;
rounding only the half-distance can differ by one at an exact half when the left coordinate is odd,
so the port applies half-to-even conversion after adding the left coordinate.

For each candidate, the active source searches left and right until **any** state posterior exceeds
the following strict thresholds:

| Stored range | Strict source threshold | `CIs` rows |
| --- | ---: | --- |
| labelled 95% | `> 0.995` | 3, 4 |
| labelled 99% | `> 0.999` | 0, 1 |
| HMM breakpoint | lattice transition midpoint, repaired into the 95% range | 2 |

The threshold values differ from the rounded conceptual wording in the manual and are therefore
recorded verbatim in result capability metadata. Linear searches clamp to alignment endpoints.
Circular searches wrap, replace sentinel endpoint values using the first/last information-rich
coordinates, and may produce an interval whose left coordinate is numerically greater than its
right coordinate.

## Ten-value output contract

`PolishBP` maps the two intervals into `CIOut(0 ... 9)` as follows:

| Index | Meaning in supplied source |
| ---: | --- |
| 0 | beginning, 99% left bound |
| 1 | beginning, 99% right bound |
| 2 | beginning, selected HMM breakpoint |
| 3 | ending, 99% left bound |
| 4 | ending, 99% right bound |
| 5 | ending, selected HMM breakpoint |
| 6 | beginning, 95% left bound |
| 7 | beginning, 95% right bound |
| 8 | ending, 95% left bound |
| 9 | ending, 95% right bound |

`MatchBPtoCI` negates all five values when the reported breakpoint is outside the selected interval.
`PolishBP` may later accept that nearby result and restore positive signs. Negative values and `-1`
therefore carry workflow meaning and must not be represented merely as unsigned coordinates.

## Repositioning dependencies

The HMM output is not the final result. The remaining `PolishBP` body includes several dependent
decisions that a faithful port must implement together:

1. Match beginning and ending independently, then reject implausibly distant alternatives relative
   to half the original event span (using circular distance when enabled).
2. If only one distinct HMM transition is available, assign it to the closer reported endpoint and
   invalidate the other interval.
3. Scan the selected ranges against the current triplet's accumulated `MissingData`, moving
   breakpoints to the first usable boundary and altering HMM midpoint values when missing runs
   intervene. At this point in the supplied cycle that mask contains input missing runs plus every
   prior `ModSeqNumY` erasure affecting one of the three representatives, but not the current event.
4. Preserve the large reassortment-only segment-boundary branch separately; it is not part of the
   ordinary RDP web scope unless segmented analyses are added.
5. Revert both polished breakpoints if either inside or outside region has fewer than three usable,
   differing triplet sites (`Module4.bas` 10817–10884).
6. Move a final breakpoint off ASCII-46 sites, retaining the source's repeated `Seq2` comparison in
   both final searches (`Module4.bas` 10901–10951).

The current `CheckEnds` uncertainty mask is related input state but is not a substitute for these
steps. `CheckEnds` answers whether a reported boundary should be warned after cyclic erasure;
BURT/`PolishBP` estimates and may reposition the boundary itself.

## Browser implementation seam

The implementation is isolated from the already-large RDP scanner and provides:

- a reusable, capacity-retaining three-symbol profile (`symbols`, information-coordinate map, and
  cumulative information-position map);
- fixed `3 x 3` transition/emission storage and contiguous lattice, backpointer, forward, reverse,
  posterior, and path buffers sized only to the expanded information-rich profile;
- a local 15-bit Microsoft C-runtime `srand`/`rand` adapter seeded with the supplied default `3`, so
  the Windows DNA5 DLL's emission starts do not depend on Emscripten libc;
- signed `BreakpointConfidence` output retaining both 95% and 99% bounds, selected HMM point,
  containment and wrap state, input versus polished coordinate, missing/gap movement, and revert
  state;
- result/project JSON, on-demand alignment JSON, CSV, and UI fields, with the nearest-informative-
  state parent bracket kept as a separate review diagnostic.
- a default-enabled scan setting that bypasses HMM work completely when disabled, reports the
  reason as `disabled`, and preserves detected or manually edited coordinates.
- a capacity-retaining one-dimensional triplet missing mask that unions source-shaped input missing
  runs with prior inclusive erased tracts, avoiding a fresh native-style `N x L` allocation.

The work is linear in the number of information-rich positions for each Viterbi/forward/reverse
pass. With three fixed states, the dominant cost is the 21 starts times up to 100 training passes.
Scratch reuse and contiguous fixed-state matrices should make the WASM implementation at least as
cache-friendly as the supplied DLL without changing the statistical path. In particular, lattice,
backpointer, and path buffers are cleared once per event, matching the DLL's `calloc` lifetime,
rather than once per Viterbi iteration.

## Validation boundary

No Session 8 source was compiled or executed. Before this implementation can be called native-
validated, parity fixtures need to cover at least:

- no information-rich sites and fewer than four information-rich sites;
- a single lattice switch, multiple switches, and an HMM interval not containing the RDP endpoint;
- linear endpoint clamps and a circular 99%/95% range that wraps coordinate 1;
- a long input-missing run and a prior representative erasure inside each side of a confidence range;
- rejection by the three-usable-sites inside/outside guard;
- role replacement followed by confidence recomputation;
- repeatability under the supplied random seed and exact-equality convergence behavior.
