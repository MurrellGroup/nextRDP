# Supplied-source BootScan discovery trace

This note records the Session 19 trace used to add automated primary BootScan discovery. It is a
conservative implementation map, not a native-parity claim. No alternate RDP implementation was consulted.

## Supplied workflow followed

The active desktop route is spread across the supplied VB application and DLL sources:

1. `BSXoverR` prepares the automated distance-mode scan and reusable pair/window results.
2. `SEQBOOT2` constructs site-major bootstrap weights. Replicate zero is the unresampled window;
   the remaining replicates use the Microsoft C runtime random-number sequence after the two
   source warm-up calls. The supplied array also generates one unused tail replicate per site;
   its random draw is retained because it advances the stream used by the following site.
3. `FastBootDist` calculates pairwise window distances. The port follows the supplied
   Jukes–Cantor distance and its saturated/unavailable sentinel.
4. `GetPltVal` assigns a vote only to a strict closest pair. Equal minima do not vote.
5. `ScanBSPlots`, `FindBeginBS`, and `FindEndBS` identify supported regions, including the source's
   two padding windows, circular origin traversal, and weaker dominant-support continuation.
6. `MakeBSEvent` converts the supported pair and outside relationship into a provisional event,
   calls the supplied daughter/parent logic, and retains a separate bootstrap support probability.
7. `BSSubSeq`, `MakeScoresBS`, and `ProbCalc` calculate the ordinary BootScan event probability.
   This is the probability that enters `XOverList` and the normal project correction path.

The browser route is therefore:

`BSXoverR` → `SEQBOOT2` → `FastBootDist` → `GetPltVal` → `ScanBSPlots` → `MakeBSEvent` →
`BSSubSeq`/`MakeScoresBS`/`ProbCalc`.

## Probability scope and the BURT question

The supplied `MakeBSEvent` calculates the BootScan event probability before later BURT/BenHMM
breakpoint polishing. `MakeScoresBS` counts pair matches over the proposed tract, estimates the
independent match probability over the information-rich sequence, evaluates the upper binomial
tail, applies the supplied informative-length multiplier, and uses the 169-site rescale/exponent
route for long tracts. A tract covering the complete information-rich sequence takes the supplied
unit-background branch. The project correction is applied afterward.

The desktop code indexes that score through `XPosDiff(BE)` and `XPosDiff(EN)`, not by testing raw
alignment coordinates directly. Consequently, when a reported boundary lies on an invariant site,
the preceding information-rich score can be included. Session 19 preserves this observable
invariant-boundary convention. It also preserves VB6 `CLng` round-to-nearest-even behavior for the
169-site match rescale and window-overlap count. Earlier WebRDP BootScan corroboration used event coordinates and
could therefore disagree at very small probabilities. BURT can move the displayed event bounds
later, but it does not recalculate the saved BootScan p-value.

The primary implementation exposes three scopes separately: the mean-support-derived bootstrap
probability, the raw `MakeScoresBS`-shaped binomial probability, and the project-corrected
probability. Native saved-output comparison is still required before describing the last digits as
parity validated.

## Pair-profile reuse

The desktop `BSXoverR` path writes reusable pair/window distance summaries to temporary cache files.
Many triplets share two-sequence pairs, so recomputing every window and bootstrap replicate for all
three pairs would discard the intended shortcut.

The WASM port uses a bounded 64 MiB in-memory pair-profile cache instead:

- the key is the unordered pair of exact working-sequence IDs;
- geometry, replicate count, seed, and alignment length are workspace invariants;
- each profile is held by `shared_ptr` while a triplet combines its three curves, so FIFO eviction
  cannot invalidate an in-flight calculation;
- cache hits, misses, evictions, current bytes, and peak bytes are visible in progress/results;
- the cache is invalidated at every erasure/re-entry round boundary.

This cache is deliberately separate from the cyclic `XOverList/XOverDefine` and `BestXOList`
shortlist. The pair-profile cache avoids repeated distance/bootstrap work inside one mutable
alignment round. The shortlist replays already evaluated whole-triplet signal summaries across
rounds when all three exact working rows are unchanged. Changed rows and new fragments are scanned
fresh.

## Browser implementation boundary

Session 19 implements the automated distance-mode primary path, strict closest-pair voting,
source-shaped supported-region discovery, provisional roles, binomial scoring, correction,
strongest-first cyclic scheduling, plot reconstruction, project restore, CSV diagnostics, and the
bounded pair-profile cache. It is active but unvalidated.

The following remain explicit work:

- native golden comparison for region endpoints, role choices, tie behavior, and p-values;
- the complete literal `FindBeginBS`/`FindEndBS` edge-warning and overlap catalogue behavior;
- BootScan tree, similarity, permutation, manual, and other submodel modes;
- full late-list event reconstruction beyond the already active non-coordinate-changing recheck.
