# Supplied-source trace: MaxChi strongest-peak recheck

This document records the Session 9 MaxChi translation made only from the supplied DNA/DNA5 DLL
source, VB application source, and RDP5 manual. No alternate implementation or reference-code
project was consulted. The implementation is intentionally labelled **source-shaped and
unvalidated**. This document describes the bounded corroboration kernel; the subsequently added
ordinary exploratory discovery path is documented separately in
[`native-maxchi-discovery-trace.md`](native-maxchi-discovery-trace.md).

No source in this checkpoint was compiled or executed.

## Native call path

The supplied workflow establishes the following path:

| Supplied source | Native responsibility | Session 9 mapping |
| --- | --- | --- |
| `VB Source/Module30.bas`, `GetCriticalDiff`, around 9598–9690 | Derive a cheap absolute match-count threshold from `LowestProb / 6`, with a `0.0001` floor | `source_critical_difference`; strict screen remains separate from χ² evaluation |
| `DNA5/MathFuncsDll.cpp`, `FindSubSeqCP` / `FindSubSeqMCPB2` | Retain sites where all three sequences are usable and at least one state differs; map alignment ↔ variable-site coordinates | `build_variable_profile`, with one coordinate vector and three compact pair-match tracks |
| `DNA5/MathFuncsDll.cpp`, `MakeWindowSizeP`, around 18577–18623 | Start from the configured fixed window, shrink for short profiles, reject half-windows below six | `fixed_window_sites = 70` and the same fixed-window fallback arithmetic |
| `DNA5/MathFuncsDll.cpp`, `WinScoreCalcP`, around 4631–4705 | Build circular rolling match totals for all three sequence pairs | `strongest_peak`, maintaining six totals and updating them in constant time per variable site |
| `DNA5/MathFuncsDll.cpp`, `MakeBanWinP`, around 4554–4629 | Ban windows starting, ending, or traversing `MissingData`; retain distinct boundary indices zero and length | `make_banned_windows`, including zero/length aliases and the trailing origin ban |
| `DNA5/MathFuncsDll.cpp`, `CalcChiVals4P3` / `CalcChiValsP`, around 6358–6407 and 6675–6712 | Apply the strict critical-difference screen, calculate three 2×2 χ² tracks, and find the maximum | `chi_square` plus the rolling `strongest_peak` pass |
| `DNA5/MathFuncsDll.cpp`, `NormalZ`, around 4210–4259; `ChiPVal2P`, around 7494–7516 | Convert one-degree-of-freedom χ² to the supplied tail approximation, including its extreme fallback | `source_normal_z` and `source_chi_p_value` |
| `DNA5/MathFuncsDll.cpp`, `MakeTWinP`, around 18644–18665; `GetACP`; `GrowMChiWin2P2` / `GrowMChiWinP2` | Start near one quarter of the original half-window, expand both sides, retain non-decreasing χ² maxima, stop at missing boundaries or after failures | `grow_peak`, preserving the `MaxX == 0 → 1` entry quirk and missing-edge ordering |
| `DNA5/MathFuncsDll.cpp`, `FastRecCheckMC2`, around 19215–19510 | Gate the triplet, scan/smooth peaks, grow one candidate, correct its probability, construct/destroy candidates, and retry | strongest admissible peak, growth, and probabilities only |
| `DNA5/MathFuncsDll.cpp`, `AlistMC3`, around 23138–23220 | Schedule MaxChi over analysis lists and feed persistent event records | representative/final-list scheduling is mapped to bounded recheck evidence; ordinary exploratory discovery uses the separately traced `MCXoverF` lifecycle |

The RDP5 manual's distinction between discovery and later corroboration remains intact. The kernel
in this document cannot create or move an event. Separately, the active ordinary `MCXoverF` port can
author a method-labelled discovery signal that enters the shared cyclic event workflow.

## Variable-site and match profiles

For each triplet `(s0, s1, s2)`, the kernel walks the compact working alignment once. A site enters
the variable profile only when all three states are nonzero and the three states are not identical.
It records the original 1-based alignment coordinate and three bytes:

| Pair track | Stored value |
| --- | --- |
| 0 | `s0 == s1` |
| 1 | `s0 == s2` |
| 2 | `s1 == s2` |

The profile also retains a prefix map from every alignment coordinate to the current variable-site
index. That map is used for native-shaped `MissingData` bans without searching the coordinate list.

## Missing data, prior erasures, and linear ends

The caller constructs one compact triplet mask from two sources already present in the RDP cycle:

1. immutable input `MissingData` produced by the supplied ten-consecutive-low-state rule; and
2. coordinates that were present in the original row but are now zero in a working representative,
   covering prior tract erasure and the unavailable side of a re-entered fragment.

`MakeBanWinP` treats boundary zero and boundary `LenXoverSeq` as distinct array cells even though
they describe the same circular join. The port retains both cells. For each mapped missing site it
bans the two native half-window ranges, records both adjacent `MDMap` edges when available, and
applies the source trailing-window ban when either origin alias is marked. Linear sequences add the
same first/last boundary markers and trailing ban used by `FastRecCheckMC2`.

The result reports input/erasure filtering and linear-edge filtering as separate booleans. It does
not claim that a missing-data-filtered profile is an ordinary negative MaxChi result.

## Rolling χ² scan

The initial half-window is `int(70 / 2 + 0.51) = 35`, then follows `MakeWindowSizeP` for short
profiles. The `GetCriticalDiff` probability is `max(0.0001, cutoff / 6)`. The kernel finds the
smallest absolute extreme match difference whose χ² exceeds that threshold and subtracts one,
because the scan itself uses a strict `difference > CriticalDiff` comparison.

For each pair track, the first left and right half-window totals are calculated once. Advancing one
boundary removes one match byte and adds one match byte on each side. Consequently, scanning all
three tracks is `O(V)` for `V` variable sites and uses constant rolling state; it does not perform
`O(V × H)` repeated half-window summation and does not retain the complete χ² profile.

The statistic for left matches `A`, right matches `C`, and half-window `H` is:

\[
\chi^2 = \frac{2(AD-BC)^2}{H(A+C)(B+D)},\qquad B=H-A,\ D=H-C.
\]

Ties retain the earlier boundary/pair because the source maximum finder updates only on a strict
increase.

## Growth and probability scopes

If the initial within-triplet probability passes the project cutoff, growth starts from
`int(H / 4 + 0.51)`, clamped to `[6, H]` and half the variable profile. Both sides expand together;
a non-decreasing χ² resets the failure count. A new window is evaluated before its `MDMap` edges
stop further traversal, matching `GrowMChiWin2P2`'s ordering.

Three values remain distinct in JSON, CSV, and the review UI:

1. `localPValue = ChiPVal2P(maximumChiSquare)`;
2. `withinTripletPValue = min(1, local × V / min(initial H, grown H) × 3)`; and
3. `correctedPValue`, which additionally applies the current valid distinct-origin triplet count
   when Bonferroni correction is selected.

`sourceRecheckHit` compares the final corrected value with the ordinary project cutoff. A profile
with no peak above the critical-difference screen is still a completed profile and is represented
separately from an unavailable profile.

## Integration points

- `UniqueEvent.maxchi_triplet_recheck` stores a recheck of the current recombinant/major/minor
  representatives. It is recalculated when role hypotheses are refreshed.
- `DistanceCorrelationEvidence.post_group_maxchi_recheck` stores the same strongest-peak evidence
  for each finalized nonrepresentative distance-list row against its role's two representatives.
- The own-role representative follows the supplied final-list loop's representative skip; its
  event-level representative evidence remains available above the table.
- MaxChi evidence is serialized into ordinary results, reloadable projects, and event CSV. The UI
  shows profile availability, bans, pair, peak, χ², windows, and all three probability scopes.
- This strongest-peak recheck never writes the cyclic discovery catalogue; the separate
  exploratory MaxChi kernel can author ordinary method-labelled discovery signals.

## Deliberate boundary and deviations of this recheck

The following behavior is not claimed by this bounded recheck. Several items are represented by the
separate exploratory kernel, but they are deliberately not folded into corroboration results:

- `SmoothChiValsP`, ordered local-peak selection, `DestroyPeakP`, the up-to-100 retry loop, and its
  `WasteOfTime` behavior inside the recheck;
- MaxChi breakpoint-pair/event construction or writes back to the event catalogue from a recheck;
- full `AlistMC3` late dispatcher/event reconstruction and every enabled-method combination;
- native lookup-table rounding (`ChiTable2`); the port evaluates the same formula in `double`;
- growth beyond half the variable profile. The browser stops before the two half-windows overlap,
  while the supplied loop relies on failure/boundary termination and has no explicit overlap guard;
- desktop-vs-browser golden numeric validation.

Those boundaries are machine-readable as `eventDiscoveryApplied: false`,
`maxChiKernelStatus: "source-shaped-strongest-peak-unvalidated"`, and
`nativeMaxChiFullRecheckComplete: false`.

## Next native-parity fixtures

Before calling either MaxChi path validated, compare the supplied desktop/DLL and this kernel at full
precision for: `MaxX = 0`; exact critical-difference ties; zero/length `BanWin` aliases; missing sites before
the first variable site; linear end windows; a grown maximum smaller than, equal to, and larger than
the initial half-window; correction factors near the native cap; and multiple peaks requiring
destruction/retry. The full corpus is listed in [`validation-plan.md`](validation-plan.md).
