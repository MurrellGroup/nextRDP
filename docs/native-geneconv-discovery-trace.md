# Supplied-source GENECONV discovery trace

## Scope and provenance

This trace covers the ordinary automated triplet GENECONV path that should become the next primary
method family in the browser port. It was derived only from the supplied RDP5 manual, VB source,
DNA5 DLL source, and companion DNA DLL source. No alternate RDP or GENECONV implementation was
consulted.

The layout-preserving manual page is PDF page 23, printed section 8.2. Section 3.3 on printed pages
4-5 supplies the automated settings and defaults. The manual says that RDP5 automated GENECONV is
triplet-only: every sequence triplet is treated as an independent alignment, monomorphic sites are
removed, all sequence pairs are screened for unusually long identical or unusually discordant
runs, and high-scoring fragments receive permutation and/or BLAST-like Karlin-Altschul
probabilities. The project-wide correction multiplies pairwise KA probabilities by the number of
pairwise comparisons in the manual's presentation. The active ordinary triplet implementation is
more specific: `MakeMCCorrection` calculates both `GCMCCorrection = choose(NS, 2)` and the shared
triplet factor `MCCorrection = choose(NS, 3)` (or the query-vs-reference analysis-list count), while
`GCXoverD(0)` thresholds and stores its automated calls with `MCCorrection`. The browser therefore
uses the active triplet/analysis-plan correction for discovery and exposes the uncorrected KA value
separately; it does not silently substitute the pair-count display factor.

## Active dispatch map

| Supplied source | Active responsibility | Port target |
| --- | --- | --- |
| `MainForm22.frm`, around 22880-23030 | Ordinary automated scan progress and `AlistGC2` shortlist dispatch | Combined bounded triplet worker scheduler |
| `Module2.bas`, around 19768-19871 | Later analysis-list GENECONV pass and `GCXoverD(0)` redo | Primary and downstream re-identification scheduling |
| `MathFuncsDll.cpp::AlistGC2`, around 22858-22940 | Batch triplets through `GCXoverDP2`, storing best probabilities and redo flags | Fused GENECONV pre-screen |
| `MathFuncsDll.cpp::GCXoverDP2`, around 21452-21605 | Build polymorphic strings, score six fragment tracks, calculate KA parameters/probabilities, and return the threshold result | Capacity-retained GENECONV kernel |
| `MathFuncsDll.cpp::FindSubSeqGCAP6`, around 7904-7960 | Expand compressed alignment bytes and label pair-12, pair-13, pair-23, or all-different polymorphic positions | Direct alignment-byte category preparation |
| `MathFuncsDll.cpp::GetFragsP`, around 18288-18545 | Build signed inner and outer concordant/discordant runs, including circular runs | Six compact run lists |
| `MathFuncsDll.cpp::GetMaxFragScoreP`, around 9791-9950 | Extend each positive fragment across penalized mismatch runs and retain its best ending | Source-order fragment maxima, optimized without quadratic rescans |
| `MathFuncsDll.cpp::CalcKMaxP`, around 9685-9790 | Solve source lambda/K, calculate the corrected critical score, and reject weak tracks | Safe bounded Newton solver plus source thresholds |
| `MathFuncsDll.cpp::GCCalcPValP2`, around 21751-21835 | Convert qualifying maximum scores to KA tail probabilities | Raw pairwise probability evidence |
| `MathFuncsDll.cpp::GCGetHiPValP`, `DelPValsP`, `MakeDeleteArrayP`, around 7261-7360 | Select the next lowest P-value and enforce the configured overlap count | Stable P/track/fragment ordering and compact coverage counter |
| `Module2.bas::MakeMCCorrection` and `Module31.bas::GCXoverD`, around 11885-11935 and 30181-30835 | Calculate pair and triplet/analysis-list factors; map fragment bounds, assign six-track roles, emit events, and run later endpoint checks | Active scan-plan correction plus method-labelled signal candidates entering shared reconciliation |

## Six fragment tracks and provisional roles

`FindSubSeqGCAP6` retains a polymorphic position only when all three sequence states are present and
not all equal. Its category is pair 1-2 equal, pair 1-3 equal, pair 2-3 equal, or all different.
`GetFragsP` converts those categories into three inner and three outer signed fragment tracks.
The active `GCXoverD` role mapping is literal:

| Track | Fragment class | Provisional recombinant | Minor parent | Major parent |
| ---: | --- | ---: | ---: | ---: |
| 0 | inner pair 1-2 | 1 | 2 | 3 |
| 1 | inner pair 1-3 | 1 | 3 | 2 |
| 2 | inner pair 2-3 | 2 | 3 | 1 |
| 3 | outer sequence 1 | 1 | 2 | 3 |
| 4 | outer sequence 2 | 2 | 1 | 3 |
| 5 | outer sequence 3 | 3 | 2 | 1 |

These are discovery roles only. As with RDP, MaxChi, and CHIMAERA, emitted GENECONV signals must
enter the shared detectable/distance/phylogenetic reconciliation before a final role call.

## Defaults and probability scopes

The supplied defaults initialize automated triplet scanning, ignored indels, G-scale 1, minimum
aligned fragment length 1, minimum polymorphisms 2, minimum pair score 2, and one overlapping
fragment. In the active automated `GCXoverD` body the three minimum-fragment predicates are
commented out, while the overlap limiter remains active. The port must expose that distinction
rather than claiming that inactive desktop controls affect ordinary triplet calls.

For each track, the source derives mismatch penalty `floor(V * G / D) + 1`, solves lambda and K,
and evaluates `1 - exp(-exp(-(lambda * score - log(K * V))))` with a high-score underflow branch.
The raw KA probability and active triplet/analysis-plan-corrected probability must remain separate.
An underflowed exact zero is retained by `GCGetHiPValP` and the normal `GCXoverD` emission branch;
only plotting applies a positive display floor.
In RDP5's stored `XOverList`, `GCXoverD(0)` writes `PVals * MCCorrection`; that corresponds to the
WebRDP value labelled **Project corrected**, not the separately displayed **Raw KA P**. A raw value
near `1e-18` and an RDP5/XOverList value near `1e-14` are therefore expected when the initial scan
plan contains roughly `10^4` correction opportunities.

GENECONV discovery maps its fragment endpoints and calculates both probabilities before the event
enters shared role reconciliation. `PolishBP`/BURT runs only after an anchor event has been chosen.
It can change the event coordinates shown later, but does not recalculate the stored discovery
fragment score, raw KA probability, or project-corrected probability.
The distinct `GCMCCorrection` pair-count value used by parts of the desktop presentation is a
documented native display boundary, not the ordinary `GCXoverD(0)` discovery threshold. Permutation
probabilities, pairwise manual scans, and alternative indel-block handling are separate workflows
and must not be silently folded into the ordinary automated kernel.

## Performance plan

The direct source extends every positive run across later signed runs, which can be quadratic in the
number of polymorphic runs. The browser port can preserve the exact stop rule and latest-ending tie
rule in `O(R log R)` per track: use signed prefix scores, next-strictly-smaller prefix indices, and a
rightmost-maximum range query. One category-run pass builds all six tracks. For circular input it
keeps the initial run and extends the matching terminal category run through that initial run,
exactly matching `GetFragsP`; it does not repeat the full signed-run list. Category preparation, run
construction, probability evaluation, overlap filtering, and coordinate mapping remain linear.
Buffers should retain capacity across triplets, and no complete alignment or `O(T)` triplet list
should be materialized.

## Ordinary-kernel late corroboration

The supplied later analysis-list pass again shortlists triplets through `AlistGC2` and dispatches
`GCXoverD(0)`. The browser therefore reuses the same ordinary six-track kernel for the event
representatives and each finalized nonrepresentative distance-list row. MaxChi prepares the shared
equality/missing-data workspace once; GENECONV adds no alignment-byte pass and retains only the
stable best admitted fragment plus compact workload and numerical evidence. This is
non-coordinate-changing corroboration. It is not a substitute for the separate permutation,
manual-pair, alternative-indel, or full method-stack event-reconstruction workflows.

## Fidelity boundaries for implementation

- Validate compressed-table expansion against direct alignment categories, especially missing
  states, ambiguous symbols, and the ignored-indel default.
- Preserve VB/DLL single-precision narrowing in mismatch extension and high-score probability
  branches where it changes thresholds.
- The supplied `GCXoverDP2` contains `NDiff[3] == 1` (comparison) where tracks 4 and 5 use an
  assignment when an outer-track discordance count is zero. WebRDP deliberately applies the
  intended zero guard to all three outer tracks rather than allowing a divide-by-zero/undefined
  conversion on track 3. This can differ only on that degenerate track-3 profile and is recorded as
  an explicit native defect fix, not attributed to BURT.
- Compare linear and circular origin runs, first/last-run duplication, equal P-value ordering,
  one-overlap deletion, and coordinate aliases around the first/last polymorphic site.
- The traced `Module31.bas::GCXoverD` consumes the stable lowest-P list, while supplied
  `Module30.bas`/`Module5.bas` sibling paths add a `LastMPV`/`ZZZXZ` exit after more than ten
  repeated selections. The browser currently follows the traced Module31 route; authorized
  release fixtures must establish which sibling governs each desktop scan mode before that bailout
  can be claimed or safely generalized.
- `GetFragsP` writes a zero-score terminal sentinel that can win `GetMaxFragScoreP`'s latest-equal
  endpoint update. Linear and circular native fixtures must establish its coordinate effect rather
  than treating it as an ordinary biological run.
- Keep permutation, manual pair, alternative indel, minimum-fragment filtering, and full late
  method-stack modes separately labelled until each supplied path is traced and ported.
- The kernel, cross-layer restore contract, and review plot are now source-shaped and active, but
  remain unvalidated until authorized native fixtures cover every boundary above. The browser
  overlap tree measures intended polymorphic-position coverage rather than reproducing the
  supplied helper's fragment-count-bound indexing quirk. A bounded bisection fallback on the
  Newton path prevents an unresponsive worker; `geneconvNumericalFallbackTracks` shows whether a
  dataset actually used that fallback.
