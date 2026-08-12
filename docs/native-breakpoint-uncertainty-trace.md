# Native RDP breakpoint-uncertainty trace

This trace maps the browser source to the user-supplied RDP5 VB and DLL sources. It uses no
alternate implementation. It records source semantics and remaining validation boundaries; it is
not a claim that the uncompiled port is numerically identical to the desktop build.

## Native control flow

| Supplied source | Lines | Behavior retained in the browser source |
| --- | ---: | --- |
| `Module3.bas` `ProcessEvent` | 15155–15245 | After at least one event has already been processed, RDP calls `CheckEndsVB` with the full RDP `XoverWindow`; `ChF=0` supplies the beginning warning and `ChF=1` the ending warning, then the two flags become `SBPFlag` 1/2/3 |
| `Module2.bas` `CheckEndsVB` | 29061–29095 | If a coordinate precedes the first mapped information-rich position, the beginning is moved forward for the check and the ending is moved backward; the original reported breakpoints are restored after checking |
| `threshold.CPP` `CheckEnds`, beginning | 25135–25190 | The beginning range starts one RDP window earlier in `XPosDiff`/`XDiffPos` space and ends at the adjusted beginning; any `MissingData` value in any triplet sequence raises the warning |
| `threshold.CPP` `CheckEnds`, ending | 25191–25259 | The ending branch preserves its source-asymmetric range construction: the lower bound is one RDP window before the ending and the target one window after it; any triplet `MissingData` raises the warning |
| `Module2.bas` `CheckEndsVB` linear guard | 29232–29257 | A linear beginning with fewer than `RL` preceding mapped positions, or ending without `RL` following positions, is warned even if no missing-data coordinate was encountered |
| `Module2.bas` `ModSeqNum` | 32099–32120 | Input `MissingData` is rebuilt after ten consecutive `SeqNum < 50` values, including the literal `Y-10` look-back when the counter first reaches ten |
| `Module2.bas` cyclic erasure | 19137–19166 | Erased event sequence is replaced by ASCII 46 and marked in `MissingData`; wrapped and nonwrapped tracts use separate loops |

The DLL's wrap branches contain a literal `CirF == 0` check around the two sentinel boundary
coordinates. The port retains that comparison rather than silently changing it to the apparently
more intuitive circular case. Strict `> 0`, `< LSS`, and linear-edge comparisons are also retained.

## Browser mapping

- `RdpScanner` precomputes the source-shaped long-run input `MissingData` mask once. It does not
  reinterpret every ambiguity as native `MissingData`.
- `refresh_breakpoint_context` layers prior cyclic erasures onto the three current representatives,
  using `ModSeqNumY`'s inclusive linear/wrapped coordinates, then reconstructs the current
  information-rich coordinate map after those erasures.
- `nativeCheckEndsApplied` remains false on the first event, matching the supplied
  `SEventNumber > 0` call guard; the range/reason check starts only after a cyclic erasure pass.
- The beginning and ending check-coordinate lists follow their distinct supplied `CheckEnds`
  branches. Source wrap omissions and linear edge gates remain explicit.
- A warning is decomposed into input `MissingData`, linear-edge, unavailable-profile, and erased-
  sequence reasons. Earlier event IDs are attributed only when that event erased one of the current
  representatives and its tract intersects the native check range.
- Immediate tract contact is retained separately from the broader native-range warning. The nearest
  erased tract is additionally reported in information-rich positions for review; that diagnostic
  is not presented as a statistical confidence interval.
- Event JSON, lazy alignment JSON, CSV, TypeScript contracts, and both review surfaces expose the
  same reason flags and native check range.

## Deliberate browser boundaries

- Working fragment records are mapped back to original identities for event roles. The desktop can
  carry extra sequence records directly through `MissingData`; the browser exposes fragment
  provenance separately and applies uncertainty to original representatives.
- The shared tract helper now follows the supplied active `ModSeqNumY` loops and includes both
  reported breakpoint coordinates for linear and wrapped tracts. Equal endpoints take the source's
  wrapped branch and cover the full alignment. A focused native fixture remains required to validate
  the active branch and its fragment/export consequences.
- Statistical breakpoint probability distributions and cross-method warning families remain
  separate desktop features outside this `CheckEnds` port. BURT 95%/99% intervals are now an
  independent active path mapped in `docs/native-breakpoint-confidence-trace.md`; they do not alter
  `CheckEnds` warning reasons.
- None of these paths was compiled or executed in Session 8. Native equality still requires focused
  linear/circular fixtures around the strict source boundaries.

## Cost

The immutable input missing mask is `O(NL)` once. Per-event reconstruction is conservatively
`O(EL)` over prior events and alignment length, followed by bounded reason storage. No full sequence
or missing-data mask crosses the worker boundary; the interface receives only flags, event IDs,
ranges, and counts.
