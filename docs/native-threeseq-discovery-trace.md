# Supplied-source 3SEQ discovery trace

## Scope and provenance

This trace covers the ordinary automated triplet 3SEQ path used in the browser checkpoint. It was
derived only from the supplied RDP5 manual, VB source, DNA5 DLL source, and companion DNA DLL
source. No alternate RDP or 3SEQ implementation was consulted.

The layout-preserving manual reference is printed section 8.7. It describes the same core walk:
rotate each member of a triplet through the candidate-recombinant role, keep sites where the two
candidate parents differ and the target matches exactly one, encode the two match states as `+1`
and `-1`, and find the largest excursion. The ordinary automated source does not expose a 3SEQ
window setting.

## Active dispatch map

| Supplied source | Active responsibility | Port target |
| --- | --- | --- |
| `Module5.bas` and `Module3.bas` triplet loops | Call `TSXOver(0)` in target order 3, 1, 2 | Stable local target order `{2, 0, 1}` |
| `threshold.CPP::FindSubSeqTS` / `FindSubSeqTS2` | Build the information-rich walk, parent-match counts, and initial ascent/descent bounds | One equality-profile pass reused from MaxChi/CHIMAERA |
| `threshold.CPP::Seq3PVals` / `Get3SeqPvalC` | Sum the finite hypergeometric maximum-excursion distribution | Bounded on-demand exact dynamic program |
| `Module3.bas::GetTSPVal` | Read the exact table, use `SiegmundDiscrete` outside it, and use a scaled table fallback when the approximation leaves `(0,1)` | Exact cache plus supplied approximation and bounded scaled-exact fallback |
| `threshold.CPP::CheckwrapC` | Extend the excursion through the origin, move the beginning to the next informative site, and convert an origin-wrapping interval for linear input | Source-shaped circular/linear coordinate conversion |
| `Module5.bas::SiegmundDiscrete`, `ApproxNu`, `ApproxNormCDF`, `ApproxNormPDF` | Large-profile discrete-boundary approximation | Literal source-shaped numerical functions |
| `Module3.bas::SwapRound` / `TSXOver` | Choose the lower-P orientation, swap parent roles/counts/bounds, apply low-information exits, Dunn–Šidák correction, and emit a call | Method-labelled candidate entering shared cyclic reconciliation |
| `Module3.bas::CheckSplit3Seq` / `SubPVal` | After earlier events, trim a tract at missing/erased runs and recompute a sub-tract excursion probability | Active later-round split/re-probability path |
| `Module2.bas` late dispatch / `Module3.bas::TSXOver(1)` | Rerun three target rotations, evaluate both Findall orientations, and duplicate each accepted record with inverted parents/interval | Representative and finalized-list non-coordinate-changing evidence |

## Walk and role mapping

For a call `TSXOver` receives `Seq1` and `Seq2` as candidate parents and `Seq3` as the target.
`FindSubSeqTS` adds one when the target matches `Seq1` and subtracts one when it matches `Seq2`.
All-different, monomorphic, gap/missing, and parent-equal sites do not enter the walk. The DLL
returns `Y - 1`, and `TSXOver` exits when that last index is below three, so the automated route
requires at least four information-rich sites.

The maximum descent identifies a tract matching parent two: parent one is provisionally major and
parent two provisionally minor. If the reverse walk has a strictly smaller P-value, `SwapRound`
exchanges both parents, counts, excursions, probabilities, and coordinates. An exact tie retains
the original descent. These roles remain provisional and enter the same detectable, distance,
phylogenetic, and BURT review workflow as other methods.

`TSXOver` calls `GetTSPVal` before `CheckwrapC`. Consequently a circular prefix extension can
increase the boundary excursion `nK` without recalculating its probability. The port preserves
both values: `probabilityExcursion` is the pre-wrap magnitude passed to `Seq3PVals`/`GetTSPVal`,
while `maximumExcursion` is the post-`CheckwrapC` magnitude used by the later source rejection gates
and tract evidence. Collapsing these two values would make some origin-spanning calls too
significant.

The shared MaxChi equality slots are pair 1–2, pair 1–3, and pair 2–3. Consequently the target-to-
parent-one slots are `{0,2,1}` and target-to-parent-two slots are `{1,0,2}` for targets `{0,1,2}`.
The static source contract pins these mappings because confusing member indices with pair-profile
indices silently produces the wrong walks.

## Exact probability and bounded fallback

`Seq3PVals(m,n,k)` is the probability that a uniformly ordered finite walk with `m` plus steps and
`n` minus steps reaches a maximum drawdown of at least `k`. The supplied desktop implementation
stores a four-dimensional recursive table in `Single` values and then materializes a three-
dimensional lookup table. The browser computes the same finite distribution directly: the compact
state is the number of plus steps used and the current drop below the running maximum. States that
reach `k` leave the survival set. State probabilities and each accumulation remain `float`, matching
the supplied `Single` table while reducing WASM memory and extended-precision work. The final
`1 - survival` result is cached by `(m,n,k)` and widened only for the shared result contract.

One evaluation is limited to eight million estimated transitions and dimensions of at most 4096.
Beyond that bound the supplied `SiegmundDiscrete`/`ApproxNu` route is used. If that route is not in
the open probability interval, the source's scaled-table idea is retained with a bounded scaled
exact evaluation, its post-truncation `onM/nM` (or `onN/nN`) power adjustment, and the supplied
`1e-300` floor when a positive scaled value underflows during exponentiation. Every signal records `exactProbability` and
`siegmundFallback`; the UI and checkpoint never present an approximation as exact.

## Correction and rejection gates

After the lower-P orientation wins, the active `TSXOver` path exits when `nN > 0` and `nK = 1`, or
when `nN - nM = nK`. A literal zero probability is also not emitted because the later source gate
requires the unbounded product `xpvalue > 0`. When project correction is enabled, 3SEQ uses the
supplied Dunn–Šidák value `1 - (1 - p)^M` above `p = 1e-15`, not the Bonferroni multiplication
used by the other currently active discovery families. At or below that threshold the source uses
the separate product `p × M`; the port retains that literal branch and also uses the product if the
stable exponential form returns zero. Query/reference mode supplies its group-pair × query
opportunity count as `M`.

## Performance and browser workflow

RDP categories, pair similarities, MaxChi equality tracks, CHIMAERA inputs, GENECONV categories,
3SEQ walks, and the triplet missing/erasure map are prepared in one alignment-byte pass. Each 3SEQ
target then filters the retained non-monomorphic positions without reading the alignment again.
Walk construction and breakpoint discovery are linear in retained sites. Exact probability work is
bounded and cached across triplets with identical `(m,n,k)` states; the cache is capped at 8192
entries. Plot reconstruction is on demand and downsampled after all three target walks are built.

Ordinary and later cyclic passes can emit 3SEQ candidates into the shared strongest-first event
scheduler. Signals retain target, direction, counts, excursion, raw/corrected probability, exact or
fallback route, coordinates, and provisional roles through JSON, restore, CSV, and review.

After at least one event is present, `FindSubSeqTS2` writes `XPosDiff[x]` as the number of retained
walk sites at or before alignment coordinate `x`; the compact port obtains the same mapping with
`upper_bound`. Only a candidate that first clears the normal corrected threshold enters
`CheckSplit3Seq`. The first missing/erased state encountered on the beginning and ending sides
defines two possible trimmed sub-tracts. `SubPVal` finds each sub-tract's maximum-minus-minimum walk
range and passes that excursion through `GetTSPVal`. The selected orientation is split first; if it
becomes worse than the unsplit opposite orientation, that orientation is split too and replaces it
only on a strictly lower probability. Correction and the source call gate are then evaluated again.
The exported `missingDataSplitApplied` records that this later-round branch encountered a missing
state, even if neither candidate boundary ultimately moved.

## Late Findall corroboration

The supplied late dispatcher invokes `TSXOver(1)` in all three target rotations for the event
representatives and finalized list candidates. The port's `threeseq_recheck_prepared` performs the
same target-order scans over the already prepared equality and missing/erasure profile.

`TSXOver` still requires the initially lower-P orientation to pass its low-information and first
corrected threshold before entering the Findall loop. With an earlier event present,
`FindallFlag = 1` forces `CheckSplit3Seq` for both orientations even when the first split remains
better. The first orientation retains the `p > 10^-15` Dunn–Šidák / smaller-tail product branch;
after `SwapRound`, the source uses its open-interval Dunn–Šidák expression for the second
orientation. Each orientation is tested
independently against the source gate.

For every accepted orientation, the active `Module3.bas` path calls `UpdateXOList3` once for the
ordinary record and once more for a copy whose parents and beginning/ending are exchanged (with
single-boundary warning direction swapped). The compact evidence therefore records two source list
entries per qualifying orientation, while retaining only the statistically best orientation as a
review summary. It does not materialize or merge those records into the browser's authoritative
event catalogue, and it does not move reconciled coordinates.

## Fidelity boundaries

- The active `FindSubSeqTS2`/`CheckSplit3Seq` inclusive position map, endpoint inclusion after a
  wrap, alignment coordinate zero/one aliases, and the distinction between encountering a missing
  state and actually moving a boundary still require authorized saved-output fixtures.
- The exact dynamic program is a mathematical replacement for the supplied four-dimensional
  lookup construction. Native saved-output fixtures must confirm float rounding, table-limit
  transitions, the scaled fallback, equal-orientation ties, and zero underflow.
- Linear first/last coordinates, origin wrapping, all-different sites, and the source's beginning-
  increment rule require authorized native fixtures.
- The manual's 100-permutation graphical envelopes are a separate manual display mode and are not
  silently treated as automated discovery probabilities.
- The representative/finalized-list Findall evidence is active, but native list-cap replacement,
  duplicate/event-catalogue interaction, warning-copy behavior, and full late method-stack event
  reconstruction remain explicit golden/full-parity boundaries.

The kernel, cyclic scheduler, evidence, restore, CSV, and review plot are source-shaped and active,
but remain unvalidated until authorized native saved-output fixtures cover these boundaries.
