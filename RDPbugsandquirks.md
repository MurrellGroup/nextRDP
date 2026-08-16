# RDP bugs and quirks

nextRDP currently aims to reproduce supplied native RDP behavior, including behavior that would
not be chosen for a new implementation. This ledger records candidates for a later, deliberately
incompatible cleanup. Entries describe observed source/runtime behavior, not changes that should be
made during the compatibility phase.

## Cyclic erasure preserves both breakpoint sites

The installed RDP 5.93 executable's `ModSeqNumY` behavior leaves both reported breakpoint sites in
the background sequence and omits both from the retained recombinant fragment. For a non-wrapping
event reported as `B..E`, the fragment therefore contains only `B+1..E-1`. The routine's own
comment says this deliberately leaves two sites behind to penalize later scans that read across the
erased region.

The supplied older `dna/threshold.CPP` source contradicts both that comment and the executable: its
active `ModSeqNumY` loops use inclusive `B..E` bounds. Dataset7 exposes the binary behavior exactly:
using the strict interior reproduces every captured fragment-to-background float distance, while
inclusive endpoints change the neighbour-joining topology and the second event's selected role.

A future cleanup should represent erased intervals explicitly instead of encoding the read-through
penalty as two residual nucleotides. Compatibility code follows the observed executable rather than
the mismatched source body.

Observed binary: RDP 5.93 `ModSeqNumY` followed by `MakeNJTreesP2`. Contradictory supplied source:
`dna/threshold.CPP`, `ModSeqNumY`, around lines 1834-1870.

## Primary RDP candidate boundaries prefer the preceding run

`FindNextP` locates a rolling-window peak, but the peak's information-rich site need not itself
support the candidate pair. `DefineEventP2` decrements before inspecting the site category, so a
peak immediately after a candidate-supporting run is attached to that preceding run. Only when the
immediately preceding site is not supportive does it search forward for the next supporting run.

This asymmetric rule is easy to misread and can move a breakpoint across several information-rich
sites. In long tracts it can also change the scaled mismatch count and therefore the event order.
A future cleaned-up detector should define explicitly whether a window peak belongs to the nearest,
preceding, following, or statistically strongest supporting run.

Supplied source: `dna5/MathFuncsDll.cpp`, `DefineEventP2`, around lines 17379-17409.

## An origin candidate can rewrite the terminal rolling count

When `FindFirstCOP` accepts a candidate at circular index one, it conditionally copies that
candidate pair's terminal rolling count over the high pair's terminal count. The mutation affects
the later tract walk and lives in the shared rolling-count array, so it can also be observed by a
later role rotation. This is stateful control flow hidden inside what otherwise looks like a peak
finder.

A future cleanup should return an explicit boundary condition instead of modifying another pair's
profile. Compatibility code must reproduce the mutation's observable effect.

Supplied source: `dna5/MathFuncsDll.cpp`, `FindFirstCOP`, its `X = 1` branch.

## The final information-rich site uses an N-minus-one convention

After `FindSubSeqPB3` returns `N` information-rich sites and `XOHomologyP2` builds the circular
rolling profiles, `FastRecCheckP` sets `LenXoverSeq = abs(LenXoverSeq) - 1`. `FindNextP` therefore
does not use the last information-rich site as an ordinary peak-search starting position, and
`ProbCalcP2` receives `N - 1` as its opportunity-length factor. During circular tract traversal,
index zero aliases the final site, so the nucleotide information is not simply deleted; search,
traversal, and probability calculations use different effective bounds.

This is a fragile indexing convention rather than a clean zero- or one-based representation. A
future implementation should use all `N` sites consistently and recalibrate/validate the resulting
probabilities as an intentional compatibility break.

Supplied source: `dna5/MathFuncsDll.cpp`, `FastRecCheckP`, the assignment immediately after the
`AvHomol` calculation and the subsequent `FindNextP`/`ProbCalcP2` calls.

## Long tracts are quantized to 169 sites

For detected regions of at least 170 information-rich sites, RDP rescales the mismatch count to a
169-site tract, rounds it to an integer, calculates the binomial probability at that quantized
count, and raises the result to `original_length / 169`. A one-site breakpoint change can cross a
rounding boundary and alter a very small p-value by orders of magnitude. In supplied Dataset6, a
boundary shift changed the scaled mismatch count from 29 to 28 and changed the B/G/O probability by
about 2,200-fold.

A future cleanup should replace this discontinuous approximation with a numerically stable direct
tail calculation over the actual tract length, with an explicit migration note because event order
can change.

Supplied source: `dna5/MathFuncsDll.cpp`, `FastRecCheckP`, the `XOverLength >= 170` branch.

## Candidate-cycle duplicate suppression leaks across role rotations

`FastRecCheckP` initializes `OldX` once before rotating through the alternative high/medium/low pair
assignments. It is not reset for each rotation. If the first peak in a new role rotation equals the
last peak examined in the previous rotation, the `x != OldX` guard ends that scan instead of
examining the candidate under its new role assignment.

A future cleanup should scope duplicate suppression to one candidate-pair scan or use an explicit
`(pair, position)` key. Compatibility code retains the supplied cross-rotation state because it
changes Dataset6's early event order.

Supplied source: `dna5/MathFuncsDll.cpp`, `FastRecCheckP`, `OldX = -1` before the `FindCycle` loop.

## Candidate scanning resumes two sites after the first drop

`DefineEventP2` returns the first information-rich position after the candidate-supporting run, and
the `FastRecCheckP` caller then advances once more before asking `FindNextP` for another peak. The
first two positions after a detected run therefore cannot begin the next candidate in that role
rotation.

A future cleanup should make consumed ranges explicit and test whether skipping both sites is
statistically intended. Compatibility code preserves the supplied resume position because it
changes Dataset6's selected candidates.

Supplied source: `dna5/MathFuncsDll.cpp`, `DefineEventP2`'s return value and the subsequent
`FindNextP(EndXO + 1, ...)` call in `FastRecCheckP`.

## Sub-threshold hits keep synthetic fragments alive

During cyclic scanning, `DropSeqs` removes a synthetic fragment only when its `NumRecsI` count is
zero (or the row is too short). `NumRecsI` is accumulated from the complete internal `XOverList`,
not merely from statistically significant events that can appear in the final result table. A
fragment can therefore remain in the mutable alignment solely because it participated in a
sub-threshold candidate. That hidden row can later alter `CheckMatrixP` pruning, neighbour-joining
trees, role consensus, and ultimately event order.

Dataset8 exposes this at its third event: the retained fragment of B has no reportable signal in
the preceding round, but native RDP keeps it and uses it in the role-scoring panel. Treating the
public signal list as the fragment-retention ledger deletes B prematurely and reverses the inferred
recombinant.

A future cleanup should give fragment lifetime an explicit policy, rather than coupling it to an
otherwise invisible candidate cache. Compatibility code must retain the distinction between
internal candidates and reported events.

Supplied source: `Module3b.bas`, the `NumRecsI` accumulation before cyclic rescanning and the later
`CheckDrop`/`DropSeqs` calls; `Module2.bas`, `DropSeqs`.

## Cyclic split detection ignores gaps from the loaded alignment

`ModSeqNum` freshly zeroes `MissingData` before replaying prior accepted events. It does not rebuild
that matrix from literal gaps already present in the input alignment. Consequently, `CheckSplit`
can split a candidate at an event-erased coordinate but will read straight through an arbitrarily
long input-gap run. Dataset9 exposes this distinction: treating its long Q input gap as
`MissingData` incorrectly suppresses the native A/Q/U event.

A future cleanup should use one explicit missing-data definition throughout the detector. If input
gaps are intentionally excluded from split handling, that should be a named policy rather than an
incidental consequence of clearing and replaying a side matrix.

Supplied source: `Module2.bas`, `ModSeqNum`, the `ReDim MissingData(...)` calls before step replay;
`Module30.bas`, `XOver`, the cyclic `CheckSplit2` call.

## Split testing happens only after the unsplit event passes

In a cyclic pass, RDP first calculates and corrects the probability of the complete candidate that
spans an erased tract. Only if that unsplit probability passes `LowestProb` does it sample
`MissingData` every ten alignment positions and consider splitting. Scoring each gap-bounded piece
directly can therefore create events that native RDP never examines, while rejecting every spanning
candidate loses native events whose split pieces remain significant.

When a split is taken, the left piece also omits the informative site immediately before the gap:
`ETarget = XDiffPos(XPosDiff(Z) - 1)`. Breakpoint centering can subsequently place the reported
coordinate at the retained edge of the erased tract, hiding that excluded site in the public table.

A future cleanup should detect missing-data partitions first and define the statistical test over
those observed regions directly. Compatibility code preserves the preliminary-unsplit gate and the
off-by-one site exclusion because both change cyclic event order.

Supplied source: `Module30.bas`, `XOver`, from the first `ProbabilityXOver < LowestProb` test through
`CheckSplit2`, `FindMissing`, and the two `SplitEvent` calls.

## Long split pieces omit the ordinary length exponent

Ordinary RDP events with at least 170 informative sites are reduced to a quantized 169-site count
and then have their probability raised to `original_length / 169`. The `SplitEvent` branch repeats
the 169-site count reduction but never applies that exponent. Long split pieces are therefore much
less significant than an otherwise equivalent ordinary event. Applying the shared ordinary-event
helper unchanged manufactured a Dataset9 E/I/Y residual that native RDP does not report.

A future cleanup should use the same statistically justified long-tract calculation for ordinary
and split events. Compatibility code deliberately preserves the missing exponent.

Supplied source: `Module30.bas`, `XOver`; compare the ordinary `AFact` block before the cutoff with
the probability calculation inside the `For A = 0 To 1` split loop.

## Sequence erasure and split missing-data bounds disagree

The installed executable leaves the two reported breakpoint bases readable in the background
sequence, but its cyclic `MissingData` replay treats the accepted event interval inclusively for
`CheckSplit` and `CentreBP`. Thus the compressed sequence and its parallel missing-data mask disagree
at both bounds. Dataset3 and Dataset4 expose the effect as exact one-base split endpoints.

A future cleanup should derive both views from one interval object. Compatibility code keeps the
two representations separate because merging them moves reported breakpoints.

Observed binary: RDP 5.93 cyclic split results in Datasets3 and 4. Supplied source paths:
`dna/threshold.CPP`, `ModSeqNumY`; `Module3.bas`, `CentreBP`.

## Raw RDP hits are filed by list occupancy, not detected roles

In long-winded mode, the native `XOver` routine does not file each raw hit under the recombinant
role inferred by the detector. It instead chooses the uniquely least-populated `CurrentXOver`
list among the triplet, falling back to the greatest `StoreLPV` value on a tie. `StoreLPV` starts
at one, so ordinary initial ties select `Seq1`. The chosen list determines the stored
daughter/minor/major ordering even though the choice is bookkeeping rather than biological
evidence, and that ordering can leak into later event processing.

This state belongs to the shared initial `XOverList` construction and its subsequent rewrite
lifecycle. Reinitializing and replaying it independently in each cyclic rescan does not reproduce
native behavior and can reverse otherwise correct parent assignments. A future cleanup should
store detected roles separately from list placement and make the latter irrelevant to event
selection. Compatibility work must reproduce the global lifecycle before enabling this rule.

Supplied source: `Module30.bas`, the long-winded `CurrentXOver`/`StoreLPV` selection around labels
55556-55588. Observed binary: the `AlistRDP4` initial scan and later `FinalTrim` inputs for
Dataset6.

## NOPINI is not the final public parent order

`NOPINI` records the daughter/minor/major column ordering used by the decision-score tables, but
native RDP does not copy those three columns directly into every public event record. `MakeBreaks`
revisits each retained `BestXOList` record for every member and enabled program, finds its nearest
and farthest parents from the full `FAMat`/`SAMat` matrices, then switches back to
`FMatSmall`/`SMatSmall` for part of the major/minor decision. A later `SwapFlag` consistency check
can separately exchange the alternative major/minor lists without rewriting the already processed
main list.

Consequently, two events can have the same `NOPINI` permutation and different public parent order.
The third supplied events in Dataset6 and Dataset8 both store `[2, 0, 1]`, while their public
major/minor fields require opposite interpretations. Reconstructing roles from only the final
triplet and `NOPINI`, or applying the list rewrite once per cyclic round, causes regressions in
otherwise correct early events.

A future cleanup should maintain one explicit biological-role record and derive presentation and
alternative-list order from it. Compatibility work must retain the complete per-list,
per-program rewrite lifecycle before replacing the current verified pre-`MakeBreaks` mapping.

Supplied source: `Module3.bas`, the `NOPINI`/`SwapFlag` block and calls to `MakeBreaks` and
`MakeAlternative`; `Module2.bas`, `MakeBreaks`. Observed binary: Dataset6 and Dataset8 r20 state
tables and captured `MakeNJTreesP2` inputs.
