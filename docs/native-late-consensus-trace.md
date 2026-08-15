# Supplied-source trace: active late RDP consensus

This trace began with the Session 8 primary-RDP translation and now records the MaxChi, CHIMAERA,
GENECONV, and 3SEQ corroboration addenda, using only the supplied RDP4/RDP5 DLL and VB sources. It is
an audit aid, not a claim of runtime parity. No alternate implementation was used.

## Active call shape

The active VB orchestration in `Module3.bas` has two mutually exclusive RFF=0 shapes:

- When all three role-compatibility results are available, `RetrimFlag` is set near lines
  45340–45343 and `FinalTrim` is called near line 45390 with `RFF = 0` and `RWinPP = 4`. That pass
  computes the late evidence without selected-role pruning. A later call near line 45755 uses
  `RFF = 1` after role selection.
- Otherwise the first call is skipped. The call near line 45755 receives `RFF = 0` and the selected
  `RWinPP < 3`, so it performs the complete RFF=0 path including selected-role pruning.

The browser collapses those two orchestration shapes around its explicit role reorder.
`refresh_role_hypotheses` calculates the role vote, the event is reordered if that vote changes the
recombinant, and the routine is rerun with role zero as the current selected role. The complete RFF=0
list build runs there, followed once by the selected-role block shared below the RFF guard. The
primary-RDP part of the later signal/probability revalidation then runs over the finalized lists.
The browser also exposes bounded MaxChi, CHIMAERA, GENECONV, and 3SEQ corroboration over the same
finalized rows; full method-origin dispatch and global event-catalogue integration remain a
documented boundary. Together with primary RDP, these are five active method-family rechecks.

## Session 25 cyclic identity corrections

The supplied late path does not use the public display order for its three representatives.
`ISeqs` is daughter, minor parent, major parent. Session 25 restores that internal ordering in
`FinalTrim`/`ConsensusOK` and BURT inputs while leaving JSON, UI, CSV, and project roles in the
public recombinant, major parent, minor parent order.

After cyclic erasure, `TraceSub` resolves original identities to the retained working fragments
that produced cached XOverList evidence. The port now retains exact working-triplet provenance in
the signal shortlist, assigns strongest-support aliases for all original candidates, and gives the
anchor triplet precedence. This is an identity mapping rule, not a cycle-number rule. Restored old
projects without exact working triplets keep the earlier fragment-event lookup.

The selected-role movement/rank block iterates through source `NextNo`, including working fragments.
Session 25 therefore uses the expanded population and denominator rather than only the original
alignment rows. Candidate membership remains keyed by original identity; comparisons that skip the
role representative use the reported identity rather than its transient working-fragment index.

## FinalTrim mapping

`Module2.bas` defines `FinalTrim` at line 23342 and gates this family with `If RFF = 0` at line
23390. The browser mapping is:

| Supplied lines | Native state transition | Port treatment |
| --- | --- | --- |
| 23398–23477 | Copy correlations, clear warned/inverted values, suppress duplicate candidate/pair values, set `OKSeq 5` | `raw_correlations`, `cleaned_correlations`, `duplicate_filtered`, `duplicate_cleaned_support` |
| 23479–23832 | Snapshot `OKSeq 6`; repeatedly rebuild nearest-nonrecombinant lists and remove rows until the last index stabilizes | `nearest_lists` fixed-point block, preserving `MakeRList` ascending order and `StripDupInv` swap-last order |
| 23835–23868 | Mark current rows and all representatives; append ascending candidates inside both retained nearest-neighbour limits with any cleaned `r > 0.83` | `first_expansion_added` and ordered `finaltrim_lists` append |
| 23873–24300 | Populate `OKSeq 7–14`, including explicit-zero 10/11 and event-region 14 | `FinalTrimMatrixEvidence`, `SourceDetectedTractGrid`, and match-matrix prefix calculations |
| 24303–24323 | Append remaining ascending candidates strictly closer to the anchor than both parents in both raw whole-region trees | `second_expansion_added` |
| 24331–24420 | If closest pairs differ and a role is selected, prune the three mapped lists with collapsed/raw thresholds and an all-three-`r < 0.99` removal ticket | active selected-role block and `selected_role_pruned` |
| 24425–24429 | Mark surviving rows as `OKSeq 15` | `finaltrim_membership` |
| 24430–24442 | Build `RCorrX` as the maximum un-warned, non-inverted correlation | `maximum_direct_correlation` |
| 24444–24459 | Run `CalcMatchY`, then clear/rebuild lists through `ConsensusOK` | `source_calc_match`, topology class grid, completed score and three rebuild passes |

Two order-sensitive details are intentionally retained:

1. List removal copies the final row into the removed position. The port does not replace this with
   stable filtering before the next native stage.
2. After pruning `INList(1)`, the source does not reset `Y` before the `INList(2)` loop at line
   24393. The third loop therefore begins at the final size of the second list. The port preserves
   this inherited index.

The collapsed first-list inside comparator at line 24340 is spelled
`SCMatSmall(INList(1), ISeqs(INList(0)))`; the fallback spelling uses the opposite row. The matrix is
symmetric, but the port retains the active spelling in code and documentation.

## ConsensusOK mapping

`ConsensusOK` begins at line 22079.

| Supplied lines | Native behavior | Port treatment |
| --- | --- | --- |
| 22120–22177 | Score nonrepresentatives from `OKSeq 0–6`, multiply by 2 for 15, apply `RCorrX`, then the `NS As Long` matrix sum; own representative is 1000 and the other two are zero | `ConsensusScoreEvidence.final_score`; `source_vb_long` narrows after every native assignment |
| 22178–22187 | Save all three lists and clear them | `finaltrim_lists` retained while `consensus_lists` start empty |
| 22207–22297 | Reject inconsistent `OKSeq 18` rows using both raw whole-region trees | `consensus_match_class` topology pass over every ordinary sequence |
| 22313–22483 | Primary ascending rebuild from CScore, regional match, class, raw dominance, collapsed topology, exact raw position, and zero-distance fallbacks | `primary_membership` |
| 22488–22534 | Snapshot membership; add exact six-distance equivalents or low-score direct neighbours | `equivalent_membership` with a fixed pass snapshot |
| 22545–22655 | Snapshot again; collect raw/direct distance stragglers | `straggler_membership` with a fixed widened snapshot |
| 22705–22715 | If any role is empty, restore all three saved lists | `consensus_fallback_restored` |

The role-zero exact-position fallback at line 22447 uses literal full-matrix row `0`, while roles one
and two use `ISeqs(role)`. The port deliberately retains that literal sequence-zero behavior.
`ConservativeGroup` is treated as its active default zero; no user option in the supplied workflow
enables the conservative branch.

## Shared selected-role cleanup

After the RFF=0 guard closes at line 24473, both native call shapes converge when `RWinPP < 3` at
line 24510. The active `ConservativeGroup = 0` block is mapped as follows:

| Supplied lines | Native behavior | Port treatment |
| --- | --- | --- |
| 24540–24571 | Sum all direct outside/inside distances into Single arrays and rank each representative by total movement | `movement_distance` and `representative_rank`, narrowing every accumulation back to `float` |
| 24576–24829 | Mark rebuilt-list rows that fail raw/direct outlier, clade equality, sixfold-distance, region-direction, or extreme-rank constraints | `selected_tree_cleanup_pruned`; all comparisons retain strict/equal source operators |
| 24830–24845 | Remove marked rows by copying the final entry into the current position | swap-last mutation of each `consensus_lists` row |
| 24866–24918 | Admit candidates within the remaining raw/direct cluster bounds and closer in both raw trees | first selected-tree admission branch |
| 24919–24944 | Admit strict direct+raw four-matrix inliers | unconditional fallback branch; the source condition ends with `Or x = x` and is therefore always true |

Line 24895 compares outside `FAMatSmall` with `HDF`, the maximum direct outside distance. The port
retains that mixed-matrix comparison rather than substituting `FMatSmall`.

The distance matrix is symmetric. Session 25 evaluates each unordered pair once, adds it to both
movement rows, and reuses the completed rows for representative totals. With ascending lower/upper
indices, every row still receives contributions in ascending sequence order, so Single narrowing
and source rank thresholds are unchanged. Pair-cache calls are roughly halved and six representative
re-sums disappear; symmetric caching already prevented duplicate expensive JC calculations.

## Post-group method rechecks

The native recheck block begins at line 25000. The primary-RDP branch can emit widened-threshold
records, while the six other discovery methods currently represented in the browser retain
bounded, non-coordinate-changing corroboration from their supplied ordinary statistics:

| Supplied lines | Native behavior | Port treatment |
| --- | --- | --- |
| 25008–25015 | Temporarily widen `LowestProb` to at least `LowP * 100000` and the corrected project threshold | `post_group_local_cutoff = min(1, max(project cutoff, best local P * 100000))` |
| 25087–25214 | For every role list, skip its own representative and rerun the event's originating program against the other two representatives | finalized `consensus_lists`; the RDP path calls `triplet_signals(..., enforce_cutoff=false)`, while BootScan, MaxChi, CHIMAERA, GENECONV, SISCAN, and 3SEQ retain their separate bounded evidence records |
| 25137–25139 | Dispatch `XOver` when `ProgramFlag = 0` | primary RDP profile/candidate path; no substitute implementation is used |
| 25220–25353 | On the selected list, run enabled programs other than the event's originating method | distance-mode BootScan, `FastRecCheckMC2`-shaped MaxChi, three-target `FastRecCheckChim`, ordinary six-track `GCXoverD`, fixed-region SISCAN, and three-rotation/two-orientation `TSXOver(1)` corroboration are retained without inserting new events; full native dispatcher/event reconstruction and unported families remain pending |
| 26092–26120 | Copy emitted records back into the persistent event catalogue | compact per-sequence recheck evidence records emitted/candidate/overlap counts and the best tract/probabilities; the ordinary cyclic scan remains the authoritative global signal catalogue |

The browser distinguishes four things that the native global `XOverList` representation combines:
all RDP records inside the widened threshold, records naming the list candidate as recombinant,
records overlapping the reconciled event by more than 30%, and records also passing the ordinary
corrected project cutoff. This keeps loosened-threshold traces from being mislabeled significant.
An unavailable information-rich profile is retained explicitly rather than treated as a negative
signal.

The supplied disk-only `RDP5ExcludeList` substitution path after line 25500 has no browser analogue:
the browser does not offload alignment records to an external exclude file. Every loaded sequence is
already considered in memory, while the existing masked-sequence trace path remains separately
auditable.

## BootScan, MaxChi, CHIMAERA, GENECONV, SISCAN, and 3SEQ corroboration addendum

The browser now applies the supplied `FastRecCheckMC2` statistic to the representative event triplet
and to finalized nonrepresentative distance-list rows. This is intentionally narrower than calling
the complete native program dispatcher in lines 25087–25353: it reports the strongest admissible
MaxChi peak and all probability scopes but does not insert a MaxChi-discovered event into the global
catalogue, reposition the RDP event, or run the native smoothing/`DestroyPeakP` retry loop.

The line-level numerical mapping and its explicit deviations are recorded in
[`native-maxchi-recheck-trace.md`](native-maxchi-recheck-trace.md). Result metadata therefore uses
`source-shaped-strongest-peak-unvalidated` and keeps `nativeMaxChiFullRecheckComplete` false.

Session 12 added the supplied three-target `FastRecCheckChim` statistic over the same prepared
equality/missing-data workspace. Session 13 added the supplied ordinary six-track GENECONV
fragment/KA kernel over that workspace. Both are exposed beside MaxChi for the representative
triplet and every finalized nonrepresentative distance-list row. None of these corroboration passes
moves a reconciled breakpoint or silently writes a new discovery event.

Session 14 added ordinary and cyclic 3SEQ discovery to the shared event scheduler. Session 15 adds
the supplied post-erasure `CheckSplit3Seq`/`SubPVal` discovery re-probability and the distinct late
`TSXOver(1)` representative/finalized-list route. That route preserves the initial selected-
orientation gate, forces split evaluation of both Findall orientations, applies the post-swap
correction branch, and accounts for the inverse-parent/inverse-interval list copy per accepted
orientation. It remains non-coordinate-changing audit evidence rather than global event-catalogue
reconstruction.

Session 21 adds the supplied fixed-region SISCAN route. It chooses the nearest eligible fourth
sequence from the round-wide source WPGMA context, scores the already reconciled region with the
final permutation count and shared MakeVRand prefix, and retains normal-tail, region-, window-, and
project-adjusted probabilities. It cannot move event coordinates. See
[`native-siscan-discovery-trace.md`](native-siscan-discovery-trace.md).

Session 22 adds PHYLPRO only as the supplied graphical event-review diagnostic. It is not inserted
into this post-group dispatcher or the global event catalogue because the supplied RDP5 route has
no implemented PHYLPRO significance test. Its three current-role curves and explicit source repairs
are traced separately in
[`native-phylpro-review-trace.md`](native-phylpro-review-trace.md).

## Matrix correspondence and bounded performance

| Native matrix | Browser value |
| --- | --- |
| `FMat` / `SMat` | direct whole-outside / whole-inside Jukes–Cantor distance |
| `FAMat` / `SAMat` | raw whole-outside / whole-inside midpoint-ultrametric topology rank |
| `FCMat` / `SCMat` | weak-node parent-promoted whole-outside / whole-inside topology rank |
| `*MatSmall(role, y)` | distance from the working representative for `role` to ordinary candidate `y` |

Candidate-to-candidate checks in the equivalence and straggler passes can revisit the same pair many
times. Session 8 memoized symmetric `(sequence A, sequence B)` values separately for the two regions
and three distance families. Saved tree-panel values remain constant-time; a candidate outside the
bounded panel falls back to a cached direct distance. This changes repeated-work cost, not branch
ordering or thresholds.

## Remaining boundary

- MaxChi strongest-peak, CHIMAERA strongest-target, ordinary-kernel GENECONV, and 3SEQ Findall
  corroboration are active. Their full native method dispatch and late event-catalogue
  reconstruction are not.
- Post-group rechecks for other method families not represented by this checkpoint also remain
  inactive.
- The 100-sequence event-tree panel and labelled direct-distance fallback are explicit browser
  bounds rather than native full-matrix storage.
- Compilation, runtime execution, and native golden comparison were intentionally not performed.
