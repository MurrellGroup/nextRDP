# Supplied-source trace: CHIMAERA exploratory discovery

This trace documents the Session 12 ordinary-triplet CHIMAERA discovery port. It was derived only
from the supplied RDP5 manual, VB6 application source, DNA5 DLL source, and companion `dna.dll`
source. No alternate RDP implementation was consulted. The implementation is source-shaped and
active, but remains unvalidated until it is compared with authorized native saved-output fixtures.

## Manual workflow represented

Section 8.5.1 of the supplied manual describes a three-rotation triplet method:

1. choose each triplet member in turn as the possible recombinant;
2. discard monomorphic sites and sites where neither possible parent matches that target;
3. encode a target match to parent one as `1` and a target match to parent two as `0`;
4. slide a centrally divided window over that information-rich binary string and calculate χ²;
5. infer a tract from each qualifying peak with the same boundary machinery used by MaxChi; and
6. pass significant evidence into the ordinary event-identification workflow.

The browser preserves that separation. CHIMAERA produces method-labelled discovery signals;
shared event reconciliation then groups overlapping support, evaluates all three possible roles,
polishes breakpoints when enabled, erases fixed tracts, and starts a complete new scan round.

## Active supplied call path

| Supplied location | Native responsibility | Port location |
| --- | --- | --- |
| `DNA5/MathFuncsDll.cpp::AlistChi`, around 23231–23404 | Walk analysis-list triplets and call `FastRecCheckChim` for three target rotations | `RdpScanner::scan_triplet`; `chimaera_discover_prepared` |
| `MathFuncsDll.cpp::FindSubSeqDP3` / `FindSubSeqDP6`, around 4835–4984 | Build target-specific information-rich position lists for all three rotations | `build_chimaera_target_profile` |
| `MathFuncsDll.cpp::FastRecCheckChim`, around 19823–20290 | Fixed/proportional window choice, missing/edge bans, initial χ² screen, and bounded peak loop | `discover_chimaera_target` |
| `MathFuncsDll.cpp::WinScoreCalc4P2`, `CalcChiVals5P` / `CalcChiVals3P` | Construct one binary target/parent match track and its rolling χ² profile | cached pair-match transform; `calculate_chi_profiles` |
| `MathFuncsDll.cpp::FindMChi3P` and `SmoothChiVals3P` | Raw strongest-peak order and smoothed destruction support | lazy raw-χ² heap; `source_smooth_chi` |
| `VB Source/Module5.bas`, CHIMAERA branches around 6141–6192 | Run `AlistChi`, rotate stored roles, then invoke full `CXoverA` for qualifying target calls | candidate target role plus shared source-shaped tract construction |
| `VB Source` `CXoverA` and the supplied MaxChi helpers | Grow the window, choose tract side, optimize the second boundary, destroy the peak basin, and retry | `grow_discovery_peak`, `find_side_chi`, breakpoint optimizers, and destroy helpers |

The supplied VB default `CWinSize` is 60 information-rich sites. It remains independently
configurable from the MaxChi default of 70 variable sites.

## Target rotation and binary profile

For a canonical triplet `(s0, s1, s2)`, the source rotations represented by `YP` become:

| Target | Parent one | Parent two | Cached equality used for `1` | Other accepted equality |
| --- | --- | --- | --- | --- |
| `s0` | `s1` | `s2` | `s0 = s1` | `s0 = s2` |
| `s1` | `s2` | `s0` | `s1 = s2` | `s0 = s1` |
| `s2` | `s0` | `s1` | `s0 = s2` | `s1 = s2` |

The first shared alignment-byte pass already records the three equality tracks needed by MaxChi.
CHIMAERA filters and remaps those cached bytes for each target instead of reading all three aligned
rows another three times. All-different sites have no target-parent equality and are removed;
monomorphic sites never enter the cached variable profile. The result is the same information-rich
binary domain described by the manual and `FindSubSeqDP3/6`.

## Peak, probability, and tract lifecycle

Each usable target profile retains:

- source-derived fixed-window fallback and critical match-difference gate;
- missing/earlier-erasure window bans and linear-edge bans;
- strict raw χ² peak order, with lower coordinate winning equal values;
- source-shaped symmetric window growth;
- raw χ² tail, within-triplet `raw × positions / window × 3`, and optional current-round correction;
- far-side left/right χ² tract selection and optimized second breakpoint;
- literal twelve-term/eleven-divisor smoothing used only to locate destruction basins;
- three consecutive wasted attempts as an early stop; and
- the supplied 100-peak retry bound, counted independently for each target rotation.

The implementation uses the same bounded heap adaptation as MaxChi: construct raw peak priority
once, then lazily discard entries cleared by a destruction basin. This replaces repeated complete
profile scans while preserving intended raw ordering. Workspaces retain capacity across triplets.

## Role and review boundary

The selected target is the provisional recombinant for the discovery signal. Parent order is
inferred from the target's parent-one match rate inside versus outside the proposed tract, then all
three roles enter the shared late decision-tree subset. This is intentionally auditable in project
JSON and the review card: target rotation, provisional parents, parent-one contrast, final consensus
recommendation, and any manual correction remain distinct.

The review plot reconstructs only the chosen target/parent-one χ² trace. First-round, nonfragment
signals are labelled exact; later signals are labelled original-alignment reconstructions because a
compact checkpoint does not retain every erased/fragment working-profile point.

## Secondary strongest-target recheck

The supplied secondary analysis lists invoke `AlistChi`, which calls `FastRecCheckChim` for all
three target rotations. The port now mirrors that bounded statistic for the current representative
triplet and each finalized nonrepresentative role-list row. MaxChi prepares the variable-site
equality and missing/erasure maps once; CHIMAERA filters its three target strings from those bytes,
applies the supplied window/critical screen, retains the best grown target peak, and reports raw,
within-triplet, and project-corrected probabilities. It is corroborative evidence only and never
rewrites the reconciled event coordinates.

## Deliberately unresolved modes

- CHIMAERA permutation significance modes;
- full native late-list event reconstruction and method-stack interactions beyond the retained
  strongest-target statistic;
- native lookup-table/floating-point rounding parity;
- authorized golden comparison of missing-boundary indices, destruction basins, equal raw peaks,
  preliminary parent order, and combined RDP/MaxChi/CHIMAERA event order; and
- manual plot and method-combination behavior outside the ordinary automated triplet workflow.

Those limitations are exposed in result metadata, UI copy, documentation, and the validation plan.
The ordinary enabled CHIMAERA path nevertheless reaches cyclic discovery, reconciliation, ordered
review/repair, reloadable project checkpoint, CSV summary, and accepted-event alignment exports in
source.
