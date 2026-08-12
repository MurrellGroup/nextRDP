# Supplied-source trace: MaxChi exploratory discovery

This trace documents the Session 10 ordinary-triplet MaxChi discovery port. It uses only the
supplied RDP5 manual, VB6 application source, DNA5 DLL source, and the companion `dna.dll` source.
No alternate RDP implementation was consulted. The separate `FastRecCheckMC2` confirmation path is
covered in `native-maxchi-recheck-trace.md`.

## Manual workflow represented

Section 8.4.1 of the supplied manual describes the exploratory triplet workflow:

1. discard sites that are not variable across the selected sequences;
2. slide a centrally partitioned window and calculate all three pairwise χ² plots;
3. select the highest peak, grow its window one site at each end until the score declines;
4. compare the far-side left and right χ² values to choose the recombinant tract side; and
5. repeat with the next peak until no significant peak remains.

The browser keeps this discovery stage separate from later event reconciliation. Method-labelled RDP
and MaxChi signals enter the same strongest-first round, two-shared-identities/tract-overlap support
grouping, three-role evidence, erasure, fragment re-entry, and complete rescan cycle.

## Active supplied call path

| Supplied location | Native responsibility | Port location |
| --- | --- | --- |
| `VB Source/Module5.bas`, `MCXoverF`, around 49606–50910 | Active ordinary triplet scan, probability gates, peak loop, event construction, destruction, and retry exits | `wasm/src/maxchi.cpp::maxchi_discover`; `RdpScanner::scan_triplet` |
| `Module5.bas`, around 49684–49822 | Variable-site compression, window selection, missing/erased and linear-end bans, three raw χ² profiles | `build_variable_profile`, `choose_discovery_window`, `make_banned_windows`, `calculate_chi_profiles` |
| `DNA5/MathFuncsDll.cpp`, `FindMChiP`, around 4467–4492 | Strict global raw maximum over boundary outer/pair inner loops | one-time `HeapPeak` priority ordering with lazy destroyed-cell rejection |
| `DNA5/MathFuncsDll.cpp`, `SmoothChiValsP`, around 5008–5047 | Smoothed profiles used only to find destruction basins | `source_smooth_chi` |
| `DNA5/MathFuncsDll.cpp`, `GrowMChiWinP2` and `GrowMChiWin2P2`, around 4359 and 6508 | Symmetric window growth and missing-boundary stop | `grow_discovery_peak` |
| `dna.dll/threshold.CPP`, `FindSide`, around 12927 | Far-side left/right χ² comparison | `find_side_chi` |
| `threshold.CPP`, `OptLeftBPMC` / `OptRightBPMC`, around 27085 / 27177 | One-boundary optimization with scaled χ² comparison and failure limit | `optimize_left_breakpoint`, `optimize_right_breakpoint` |
| `threshold.CPP`, `DestroyPeaks`, around 12792 | Remove the smoothed basin around a completed tract | `destroy_completed_peak_region` |
| `DNA5/MathFuncsDll.cpp`, `DestroyPeakP`, around 18662 | Remove the smoothed basin around a rejected grown peak | `destroy_rejected_peak` |
| `Module5.bas`, around 50411–50433 and 50870–50910 | Ensure the source peak lies in a completed destruction interval; clear the exact raw maximum; stop after three wasted cycles | `include_source_peak_in_destroy_region`, exact-cell clear, `wasted_attempts` |

`Module31.bas` contains an older/sibling copy with disabled optimization branches. Session 10 follows
the active `Module5.bas` path, where both breakpoint optimizers are called.

## Index and ordering contract

The source mixes zero- and one-based positions. Session 10 retains these distinctions deliberately:

- raw `ChiVals` peaks are searched at indices `0 ... LenXoverSeq - 1`;
- sequence match and missing-boundary positions are treated as `1 ... LenXoverSeq`;
- `MaxX = 0` becomes `1` for growth, while the original raw index is retained for reporting and
  exact peak destruction;
- `FindMChiP` updates only on strict `>`; equal scores therefore prefer the lowest raw boundary, then
  pair 0, pair 1, pair 2;
- the accepted-interval adjustment uses ordinary numerical absolute differences, not shortest
  circular distance, exactly as the VB branch does; and
- alignment coordinates are emitted one based from the compressed variable-site map.

The smoothing code also preserves a supplied implementation quirk. `SmoothChiValsP` sums positions
`-5 ... +6`—twelve terms—but divides by eleven. Its terminal `ChiVals[LenXoverSeq]` access is a
zero-filled padding cell, not an alias for raw boundary zero. The port allocates that padding cell and
keeps the resulting asymmetry because it can change destruction basins.

## Candidate probability lifecycle

For each raw maximum, the port keeps the three native scopes separate:

- raw one-degree-of-freedom `ChiPVal2P`-shaped tail;
- within-triplet position/profile correction, `raw × LenXoverSeq / min(initial,grown window) × 3`;
- optional current-round triplet correction.

The raw maximum first passes the source pre-growth within-triplet gate. The peak is then grown using
the active DLL's highest-χ² rule. A corrected-significant result enters side selection and tract
construction. A nonsignificant grown result has its smoothed peak basin destroyed and counts as a
wasted attempt. A stored candidate resets the wasted counter. Three consecutive wasted attempts end
the triplet scan; a separate 100-attempt cap bounds the native `Redox` loop.

## Safe performance adaptation

The supplied VB loop calls `FindMChiP` after every destruction, rescanning all three raw profiles.
The port constructs the same raw priority order once with linear-time heap construction. Each
extraction is logarithmic; entries whose raw cells were cleared by a destruction basin are discarded
lazily. Raw χ² construction itself uses rolling left/right totals on all three pair tracks, avoiding
a half-window sum at every boundary. In the combined scan, RDP categories, MaxChi pair matches and
variable prefixes, and the triplet missing/erasure map are populated during one alignment-byte pass.
Capacity-retaining workspace vectors are reused across triplets.

The original destruction functions contain index expressions that can address outside their logical
profile range under malformed/edge state. The browser port uses bounded circular traversal and a
zero padding cell; it does not reproduce undefined memory reads. This is a safety adaptation and a
golden-comparison item, not a claim of bit-for-bit identity.

## Event-role boundary

The native scan initially stores a tract and later uses desktop event state, `FindDaughter`, existing
event counts, and `StoreLPV` to arbitrate recombinant and parents. Those global desktop structures do
not have a one-call equivalent at the raw discovery boundary. Session 10 derives a provisional role
from inside-versus-outside changes in the three pair similarities, then immediately subjects all
three roles to the shared supplied-source late evidence and weighted consensus. The raw peak pair,
provisional role, final recommendation, and any user edit all remain auditable. Native golden role
fixtures are required before this part is called parity validated.

## Deliberately unresolved modes

- MaxChi permutation significance modes;
- the manual pair/doublet workflow whose variable positions are derived from a larger alignment;
- native lookup-table and floating-point rounding parity;
- full desktop method-battery interactions with MaxChi role voting; and
- native golden comparison of smoothing aliases, destruction basins, missing-boundary edges, and
  combined RDP/MaxChi event order.

These limitations are exposed in result metadata, documentation, and the validation plan. They do
not prevent the ordinary enabled triplet workflow from reaching scan, reconciliation, ordered
review, project checkpoint, CSV, and accepted-alignment export in source.
