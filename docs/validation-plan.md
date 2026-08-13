# Validation plan

No compilation, bundling, type-checking, preview, or project runtime execution was performed in
sessions 1–13. Session 7 resolved lockfile metadata without installing dependencies or running
package scripts; sessions 8–13 used only source/PDF inspection and static text-contract checks. The
first runtime checkpoint must use a disposable, resource-limited environment and tiny inputs.

## Static/build gate for a later authorized session

1. Configure TypeScript with strict checking and resolve every diagnostic.
2. Compile the C++ target with Emscripten warnings enabled.
3. Confirm every generated export matches `rdp_api.h` and the worker module interface, especially
   the `v1alpha15` analysis-scheme/reference-group/method-aware
   signal/MaxChi/CHIMAERA/GENECONV/3SEQ-discovery/event restore argument order,
   co-recombinant group pointers, and bounded signal-plot/event-alignment/event-tree JSON entry
   points.
4. Check the production bundle contains relative URLs and both `.mjs` and `.wasm` assets.
5. Load the smallest three-sequence fixture before any multi-event or production input.

The GitHub Pages workflow encodes steps 1–4 as remote gates: `npm ci`, a header/CMake/worker
ABI-version-schema check, strict TypeScript checking, the single-worker Emscripten/Vite build, and
`scripts/verify-pages-output.mjs`. The artifact verifier requires
top-level `index.html`, `.nojekyll`, `wasm/rdp-core.mjs`, and a nonempty binary with the WebAssembly
magic header; it rejects root-relative HTML assets and symbolic/hard links. This workflow has not
yet been executed, so its first run remains part of the authorized runtime checkpoint.

## Parser corpus

Use equivalent three-to-ten-sequence alignments in FASTA, GDE, CLUSTAL, PHYLIP sequential,
PHYLIP interleaved, NEXUS, and MEGA. Compare names, normalized sequences, missing-site counts,
pair identities, auto-mask set, and triplet count. Include malformed lengths, quoted NEXUS names,
ambiguity codes, gaps, duplicate names, and concatenation markers.

## Automated query-vs-reference parity corpus

Capture the native `ReferenceList`, `RefNum`, `QuerySeqs`, `RefSeqs`, `AnalysisList`, `TripListLen`,
and `MCCorrection` from the supplied `MakeAnalysisListQvR` path, then compare the browser plan,
cursor order, per-round correction, emitted signal roles, and saved assignments.

- One query with two singleton reference groups; assert the sole triplet and reference/reference/
  query source order.
- Multiple queries and multiple records in each of three groups; compare the reference-pair outer
  loop, query inner loop, exact cross-group record-triplet total, and the distinct
  `choose(RefNum, 2) × QNum` correction.
- Same-group reference pairs, query/query pairs, three-reference triples, and triples containing no
  query; assert none enters RDP, MaxChi, CHIMAERA, GENECONV, or 3SEQ discovery.
- Documented `REF-A<name>`/`REF-B<name>` FASTA prefixes, mixed case and separators, ungrouped `REF`
  records, explicit numeric edits, and all-queries reset. Compare the editable browser mapping with
  native prompts; do not treat heuristic auto-group differences as numerical-kernel mismatches.
  Assert that browser detection preserves visible curation, while an explicit “Enable all” reproduces
  the native auto-assignment routine's mask-clearing side effect when desired.
- Masked, disabled, ambiguity-heavy, and short query/reference records around the
  `max(5, RDP window)` eligibility boundary; compare active groups, workload, correction, and the
  dedicated no-eligible-constrained-triplets error/termination.
- A signal whose query is recombinant, a signal whose first reference is recombinant, and one whose
  second reference is recombinant; verify the scheduling constraint never becomes a role constraint.
- An erased query tract and an erased reference tract that each create a retained fragment; verify
  inherited group identity, distinct-origin exclusion, next-round cursor order, exact working
  workload, and the port's unique-query-origin correction adaptation against native current-record
  behavior.
- Correction/event-zero repair, rejection, pending-prefix save, and completed save/reload fixtures;
  assert every rebuild refreshes roles/cursors/totals and `v1alpha15` restores the same scheme and
  full reference-group vector.
- `v1alpha1`–`v1alpha10` imports must remain fully exploratory regardless of absent or stray fields;
  v1alpha11–v1alpha15 malformed group arrays must normalize safely or fail at the core eligibility gate.
- Counts around the native cap and synthetic overflow-scale arithmetic; verify saturation occurs
  before multiplication and the cap is applied without changing the lazy scheduled-work count.

## Primary RDP parity corpus

For each case, record the supplied desktop application’s active mask, information-rich categories,
rolling counts, candidate start/end, local p-value, round correction factor, corrected p-value, and
provisional roles, then compare with WASM at full numeric precision.

- One clear linear recombinant tract.
- One tract crossing a circular origin.
- Events beginning/ending on the first or last information-rich site.
- Endpoint fixtures that verify the port's inclusive `ModSeqNumY` linear and wrapped erasure,
  fragment re-entry, CheckEnds attribution, and final FASTA output against the active source path.
- Even and odd window settings.
- Exactly 29, 30, and 31 information-rich sites around the default window boundary.
- Earlier erased tracts at exactly 29, 30, and 31 RDP information-rich triplet positions from a
  later breakpoint, including long invariant raw-coordinate gaps. Compare the distinct supplied
  beginning and ending `CheckEnds` ranges, strict `>0`/`<LSS` gates, circular sentinel omissions,
  linear edge warnings, immediate contact, and the requirement that erased-event attribution affect
  a current representative.
- Input missing runs of 9, 10, and 11 `SeqNum < 50` characters, including the native `Y-10`
  look-back, an isolated ambiguity that should not set `MissingData`, and a missing run overlapping a
  prior erased tract; compare reason decomposition without changing the aggregate `SBPFlag` result,
  then verify BURT sees the accumulated triplet mask before the current event is erased.
- Confirm `CheckEnds` is not applied to event zero and begins only after the supplied
  `SEventNumber > 0` guard, including a first event near a linear edge with input missing data.
- Equal category counts with unequal raw pair identities.
- Highest average support immediately below and at 0.7.
- Region lengths 168, 169, and 170.
- Gaps/ambiguity inside and outside a candidate tract.
- Mask sets leaving 3, 4, and many active sequences.
- P-values immediately around corrected `0.05`.
- Active triplet counts around the supplied correction cap.

## Cyclic/fragment parity corpus

- Two non-overlapping events whose first-pass probability order is known; confirm only the strongest
  is reconciled before the second full pass.
- A weaker first-pass signal that disappears after the strongest tract is erased.
- A later event detectable only through a re-entered fragment.
- Direct and fragment-assisted copies of the same original triplet emitted in different working
  orders; confirm canonical deduplication retains the strongest call and provenance.
- Multiple fragments from one origin; confirm no same-origin pair enters a numerical triplet and the
  Bonferroni opportunity count excludes those combinations.
- Fragments immediately below/at `max(5, window, ceil(length/100))` usable states.
- Duplicate same-origin fragments, the 255th/256th/257th retained copy, and visible cap status.
- Alignments of length 99,999 and 100,000 to verify the supplied fragment-re-entry cutoff.
- Termination by no significant signal, no newly erased sites, and fewer than three active origins.
- Project save/reload where event 2’s anchor uses a fragment from event 1; compare replayed evidence
  and each signal’s correction factor.

## Secondary-phase parity corpus

For each native example, capture strongest-signal order, event support, all three role hypotheses,
distance-correlation inputs/results, tree membership, role contributions, and state before/after an
accepted correction or rejection.

- Two signals sharing exactly two original triplet members at overlap just below, at, and above 0.3.
- Linear and origin-spanning tract pairs with identical symmetric overlap.
- A transitive `FindSets` chain that exercises each two-roles-imply-the-third closure.
- Competing anchors where the stronger event claims a support signal first.
- A masked sequence with corrected-significant evidence and one with trace-only evidence.
- A role swap followed by later-event re-identification.
- A breakpoint correction that changes which later signals survive.
- An automatic group with one false inclusion and one missed descendant; edit both directions,
  accept, rebuild, and verify that erasure/project replay/final FASTA use the manual group while the
  automatic baseline remains unchanged.
- A rejected early event: verify its tract is restored, its record remains rejected/fixed, and all
  later events are rediscovered without re-emitting that same fixed call.
- Four flanking boundaries that hit the 60th information-rich site and four that stop at the other
  event edge; compare every coordinate with `MakeBPosLR`/`MakeSDMP2`.
- Direct correlations immediately around `P = 0.05`, zero-variance vectors, and a candidate failing
  the strict >10-site `MakeGoodC` gate.
- Each of the three swaps and both cyclic relabellings; compare selected inversion class.
- Every `RCorrWarn` dominant/triangle branch, including both-breakpoints-warning XOR behavior.
- All six `MakeINList` outside/inside closest-pair mappings plus the unchanged-pair no-map case;
  compare every active `MakeACOR` inequality.
- First-two-correlation values around `0.95`/`0.98` that distinguish native `corc == 2` from `>= 2`.
- Positive-only, inverse-only, and mixed candidates around `r = 0.83`; verify `StripDupInv` removes
  only inverse-only membership while retaining diagnostics.
- One candidate/pair with direct `r` just below, at, and above `0.83` in one versus multiple role
  lists; compare the opening `FinalTrim` duplicate counts and verify a duplicated pair is cleared
  for every occurrence of that sequence in every role list, even where that occurrence was already
  warning/inversion-cleared or did not itself exceed `0.83`.
- `FinalTrim` matrix-score fixtures for each strict-less, tied-less-or-equal, greater-than-one, and
  greater-than-both branch in `OKSeq` 7 and 8; include the raw-tree `+1` versus collapsed-tree `+2`
  tied rewards and the repeated-closest-pair `0.5` role modifier.
- Whole-tract and breakpoint JC fixtures producing every `OKSeq` 9/12/13 positive value and
  `-0.25`/`-0.5`/`-1` penalty; cover `INList(0)` weaker-positive suppression, `INList(0/2)` positive
  halving, same-pair scaling, `RCorrWarn`, and the strict `<3` breakpoint-distance gate.
- `FindActualEvents` fixtures with overlap immediately below/at/above one third, multiple qualifying
  direct events, inverse-only parental entries, linear/wrapping/equal-endpoint tracts, and candidate
  tract intersection. Compare chosen signal/tract, all seven `MakeMatchMatX2P` cells, the 30-site and
  `0.75` saturation gates, every `OKSeq 14` reward/penalty, and the bare-`CompMat` index quirk.
- `CheckPatternX` fixtures for each `1` and `0.5` informative-state branch, missing states, tied
  representatives, all three circular spans, and equal endpoints; compare per-region pattern counts
  and the three `OKSeq 3` shares.
- `FinalTrim` `OKSeq 6` fixtures with initial one-entry lists, multiple inverse-only swap removals,
  nearest raw/collapsed distances at the source scale cutoff, correlations at `0.83`/`0.95`/`0.99`,
  found/unfound direct events, one-side versus two-side rejection, every four-region breakpoint-veto
  contribution, and one-removal-per-iteration convergence.
- `ConsensusOK` score fixtures for every `OKSeq` 0–6 multiplier combination, current versus other
  representative cells, `RCorrX`, both `OKSeq 15` states, and positive/negative/fractional 7–14
  subtotals. Compare the exact `NS As Long` half-to-even result after each assignment and the
  zero-valued sub-one branch separately from the raw matrix subtotal.
- `FinalTrim` final-list fixtures covering no nonrecombinant candidates, a nonempty candidate list
  with no source-scale neighbour, both ascending expansion passes, every INList selected-role gate,
  swap-last deletion, and the inherited INList(2) loop index after INList(1) shrinks.
- `ConsensusOK` rebuild fixtures for both primary thresholds, class 1/2 dominance, collapsed-tree
  admission and identical-parent veto, the literal role-zero direct-matrix fallback, six-distance
  equivalence, low-score direct widening, straggler collection, and one-empty-role restoration.
- Shared selected-tree cleanup fixtures for each anchor-outlier side, both parent orderings, raw
  equality mismatch, sixfold direct cutoffs, cross-region direction checks, representative ranks at
  0.75/0.95 and a 0.5 rank delta, swap-last deletion order, mixed `FAMatSmall <= HDF` admission,
  and the always-true strict four-matrix fallback.
- Primary-RDP post-group fixtures with the representative skip, an unavailable profile, no emitted
  signal, a noncandidate-recombinant emission, candidate traces below/above 30% event overlap, and
  corrected P values around the ordinary cutoff. Compare the widened native `LowestProb` boundary
  with `LowP * 100000`, emitted record counts, best tract selection, and circular wrapping.
- `CalcMatchY` fixtures with tract-variable counts 3/5/29/30 to cover VB half-to-even rounding and
  the 15-site cap; compare all four 40-variable-site walks, pseudo breakpoint offsets, the 160-site
  map bound, every signed `MakeVarMap2` state, circular rolling sums, and `OKSeq 17` products.
- Six-checkpoint fixtures immediately below/at/above 0.6/0.65/0.7/0.75/0.8/0.9; compare raw
  `OKSeq 18`, both `BreaksExist` flags, zero-sentinel sampling, and each `ConsensusOK` topology-order
  rejection with and without a candidate outside the bounded tree panel.
- Six-region JC matrices with gaps, saturation, exactly 9/10 comparable positions, and both whole
  tract partitions.
- Four-to-ten-taxon NJ ties, negative limbs, and known splits; compare native patristic matrices.
- Bootstrap branches at 4/10, 5/10, and 6/10 support.
- A candidate satisfying each paired-tree affinity check, plus one that passes raw trees but fails
  after collapse.
- Active tree panels of 99, 100, and 101 records to verify full-panel/fallback labeling.
- For each usable tree, verify the on-demand edge list has `nodeCount - 1` unique undirected edges,
  leaf nodes map to the saved working-panel order, internal support/collapse flags match the numeric
  summary, and fragment leaves retain their original sequence and source-event IDs after reload.
- Every detectable/distance/phylogenetic truth-table combination for the two-of-three group.
- Role fixtures where each displayed metric wins a different role; compare scores, full/half
  contributions, the collapsed-tree special condition, mapped weights, and the subset winner
  separately from the full desktop winner.

The primary RDP milestone is complete only when every mismatch is fixed or documented as an
approved representation/performance change.

## MaxChi discovery and confirmation parity corpus

Capture `LenXoverSeq`, `HWindowWidth`, `CriticalDiff`, all three pair profiles, banned positions,
`MaxX`/`MaxY`, initial and grown χ², `WinWin`, raw tail, within-triplet `xMPV`, and final corrected
probability from the supplied desktop/DLL path. MaxChi discovery is active in the source checkpoint;
compare these values before calling that path parity validated or using it for production decisions.

- Fixed 70-site windows with 6, 7, 69, 70, 71, 92, and 94 usable variable sites, including both
  `MakeWindowSizeP` fallback branches and `HWindowWidth <= CriticalDiff` rejection.
- All three match-pair tracks as the unique strongest peak, plus exact χ² ties to verify native
  boundary-first/pair-first ordering.
- Match differences immediately below, at, and above `CriticalDiff`; verify the DLL's strict
  comparison and the `LowestProb / 6` screening threshold.
- Peak positions `MaxX = 0`, `1`, `HWindowWidth`, `LenXoverSeq - 1`, and positions wrapping the
  origin; retain the source `MaxX = 0 → 1` growth quirk separately from the reported peak.
- `MakeBanWinP` fixtures with missing data before the first variable site, on the first/last
  variable site, between adjacent variable sites, and across the origin. Compare both `BanWin[0]`
  and `BanWin[LenXoverSeq]`, plus the rule allowing a half-window to end on but not traverse a
  missing-data boundary.
- Repeat the same boundary cases using only a prior erased tract, only native ten-character
  `MissingData`, and their union; verify unrelated prior events do not enter the triplet mask.
- Linear sequences with otherwise strongest peaks crossing each end, beside circular controls.
- Growth fixtures where the best window shrinks below the initial half-window, stays equal, grows,
  encounters `MDMap`, and exhausts `MaxFailCount`. Compare denominator selection in `xMPV`.
- χ² values around the `NormalZ` 5.9999999 branch and the `ChiPVal2P` zero fallback, retaining full
  numeric precision rather than comparing formatted UI values.
- Bonferroni factors 1, 2, the active round triplet count, and the supplied cap; verify raw,
  within-triplet, and corrected probabilities remain distinct fields.
- Eleven-position smoothing fixtures that isolate the supplied `-5 ... +6` sum divided by 11,
  including the zero-padded `LenXoverSeq` cell and positions on both sides of the circular origin.
- Multiple-peak fixtures with equal raw χ² across boundaries/pairs; compare strict `FindMChiP`
  ordering, the chosen grow/side/optimization result, and the next surviving peak after each
  completed `DestroyPeaks` basin.
- Rejected grown peaks that exercise `DestroyPeakP`, exactly two versus three consecutive wasted
  attempts, an accepted hit that resets `WasteOfTime`, and 100 versus 101 `Redox` attempts.
- Left- and right-tract fixtures for `FindSide`, `OptLeftBPMC`, and `OptRightBPMC`, including
  one-step missing boundaries, numerical (not circular-shortest) inclusion of `pMaxX`, linear invalid
  wraps, and converted alignment coordinates.
- Combined rounds where RDP and MaxChi each supply the strongest corrected signal, an exact
  cross-method tie falls through stable signal order, both methods group one event, and one method's
  first-round signal disappears after the other's tract is erased.
- Representative and finalized-list rechecks with representative skip, unavailable profile, no
  critical peak, uncorrected trace, and corrected hit. Confirm that this separately labelled recheck
  never replaces detected or manually edited event coordinates.
- Confirm method-aware profile JSON reports pair identity for RDP and χ² for MaxChi, and that CSV
  records detection methods plus the complete MaxChi anchor attempt/side/peak/window/flank/filter
  diagnostics.

## CHIMAERA discovery parity corpus

Capture each `YP` rotation's `LXOS`, `XDP` position map, binary score string, `HWindowWidth`,
`CriticalDiff`, `BanWin`/`MDMap`, every raw χ² value, `MaxX`, initial and grown χ²/window, flank
scores, tract boundaries, candidate roles, raw tail, within-triplet probability, and final corrected
probability from the supplied `AlistChi` → `FastRecCheckChim` → `CXoverA` path. CHIMAERA discovery
is active in this source checkpoint; compare these values before calling it parity validated.

- One fixture for each target rotation where only its target/parent-one equality is present, only its
  target/parent-two equality is present, both alternate across the alignment, and neither is present.
  Confirm monomorphic and all-different sites are excluded and alignment coordinates match
  `FindSubSeqDP3/6` exactly.
- Fixed 60-site windows with 6, 7, 59, 60, 61, and fallback-sized information-rich profiles,
  including `HWindowWidth <= CriticalDiff` rejection and both `MakeWindowSizeP` fallback branches.
- Target profiles whose raw score strings are complements. χ² profiles must match while reported
  parent-one inside/outside rates and provisional parent ordering reverse appropriately.
- Exact raw χ² ties at several positions independently in each target. Confirm lower position wins,
  target rotations remain in source order, and destruction of one target never mutates another.
- Missing-data and earlier-erasure runs before/on/between/after target information-rich sites,
  including positions that are variable for MaxChi but absent from one CHIMAERA target. Compare the
  target-specific `XPosDiff` mapping, both native boundary aliases, and linear end bans.
- Circular and linear peaks around coordinate one, raw position zero, both alignment ends, and a
  wrapped tract. Compare reported peak, growth coordinate, tract endpoints, and wrap flag.
- Growth fixtures that retain, enlarge, or reject the starting window; hit a missing boundary; and
  exhaust `MaxFailCount`. Compare raw, within-target (`Len/H × 3`), and project-corrected scopes.
- Left/right `FindSide` fixtures and `CXoverA` breakpoint optimization fixtures with exact flank ties,
  one-step missing boundaries, and source peak inclusion in the completed destruction interval.
- Completed and rejected smoothed basins around the origin, including the supplied twelve-term/
  eleven-divisor behavior, exactly two/three wasted attempts, an accepted call resetting the counter,
  and 100/101 attempts independently for each target.
- Parent-order fixtures where the target matches parent one more inside, more outside, or equally.
  Compare the preliminary target/parents separately from the shared late role recommendation.
- Combined rounds where CHIMAERA supplies the strongest signal, ties RDP or MaxChi, joins their
  support event, disappears after another method's tract erasure, or becomes detectable only through
  fragment re-entry. Compare stable method order and all workload counters.
- Confirm on-demand plot JSON contains one target/parent-one trace, retains both breakpoints plus the
  chosen peak under thinning, and labels later erased/fragment profiles as reconstructions.
- Confirm CSV and project JSON retain the full CHIMAERA anchor trace; save/reload must reproduce the
  method, target, evidence, counters, event anchor/support IDs, and review state exactly.
- Representative and finalized-list `FastRecCheckChim` rechecks covering all three target winners,
  target ties, representative skip, unavailable profile, usable profile without a critical peak,
  uncorrected trace, and corrected hit. Compare target count, information-rich length, windows,
  peak coordinate, χ², and all probability scopes against the supplied secondary `AlistChi`
  calls; assert the result never changes detected or manually edited event coordinates.
- For a finalized-list row, instrument triplet preparation and assert MaxChi plus CHIMAERA consume
  one shared alignment-byte/equality/missing-data pass. Project JSON, review badges/tooltips, and CSV
  must agree on both rechecks and retain the related-method/non-independence caution.

## GENECONV discovery parity corpus

Capture `LenXoverSeq`, all six `NDiff` values, every signed `FragSt`/`FragEn`/`FragScore`
run, `MissPen`, `MaxScorePos`, `FragMaxScore`, `HiFragScore`, lambda, K, critical score,
raw KA P-value, lowest-P selection order, overlap counters, mapped breakpoints, and preliminary
roles from the supplied `FindSubSeqGCAP6` → `GetFragsP` → `GetMaxFragScoreP` →
`CalcKMaxP` → `GCCalcPValP2` → `GCXoverD` path. Ordinary ignored-indel GENECONV discovery
is active in this source checkpoint; compare these values before calling it parity validated.

- One fixture uniquely favoring each of inner tracks 0–2 and outer tracks 3–5. Include
  all-different sites and confirm the literal six-way recombinant/minor/major mapping.
- Direct compact-state categories versus supplied compressed `FSSGC` expansion for A/C/G/T,
  ambiguity codes, isolated and long missing runs, gaps, and all-equal sites. Confirm the active
  ignored-indel default and polymorphic coordinate map.
- Linear first/last positive fragments and circular first/last category runs that match, differ,
  merge as positive outer runs, or cover the origin. Confirm the initial run remains present, the
  matching terminal run includes its origin copy, the full signed list is not repeated, and compare
  fragment run order, endpoints, signed scores, and sentinels.
- Every positive fragment start followed by scores that remain positive, hit exactly zero, become
  negative by one, recover after a nonnegative step, and tie the prior maximum. Confirm the source
  stops only below zero and the latest ending wins exact score ties.
- Compare the optimized prefix/next-lower/rightmost-range implementation with a literal bounded
  quadratic reference loop over the supplied run arrays for every fixture.
- Mismatch scales 1, 2, and a high valid value with discordant counts near float/integer rounding
  boundaries. Record the DLL's single-precision `LTG`/`MPen` narrowing and every resulting
  score/critical change.
- `NDiff` values 0, 1, 2, 3, one track equal to all polymorphic sites, and cases immediately
  around the native `MinDiff < 3 && MaxDiff > MinDiff * 10` skew filter.
- Lambda/K fixtures where native Newton converges immediately, needs many iterations, oscillates,
  crosses one, or becomes non-finite. The ordinary converged route must match at full precision;
  fallback cases must be labelled and compared as the documented browser safety adaptation.
- KA scores below/equal/above zero, immediately around 32 and 700, and extreme values exercising
  the source underflow division. Compare raw values before any project correction or display floor,
  and confirm an underflowed zero is retained rather than promoted to a minimum positive P-value.
- Capture both `MCCorrection` and `GCMCCorrection` for unconstrained and query-vs-reference runs.
  Confirm ordinary `GCXoverD(0)` detection uses the triplet/analysis-list `MCCorrection`, while
  native screens that divide by the pair-count `GCMCCorrection` remain a separately reported
  presentation boundary; the browser's corrected discovery P must follow the active call.
- Critical values below/at/above four and integer fragment scores just below/equal/above the
  noninteger threshold. Confirm both strict `score > 3` and strict `score > critval` gates.
- Exact equal raw P-values within one track and across tracks. Confirm stable track 0→5, then source
  fragment order; compare lowest-P iteration with `GCGetHiPValP`.
- Nonoverlapping, touching, nested, crossing, and origin-wrapped fragments with overlap allowances
  1, 2, and 3. Compare the port's polymorphic-position coverage interpretation with
  `DelPValsP`/`MakeDeleteArrayP`, recording the supplied fragment-count-bound quirk separately.
- More than ten selection-loop visits with distinct and exactly equal raw P values. Capture the
  active `LastMPV`/`ZZZXZ` escape behavior separately from the intended overlap allowance, including
  the repeated visit to a just-accepted fragment.
- Linear and circular maxima whose latest equal score lands on the zero-score terminal sentinel;
  compare `MaxScorePos`, the `LenXoverSeq` endpoint alias, and the resulting alignment coordinate.
- Breakpoints before the first and after the last polymorphic site, adjacent polymorphic sites,
  origin wrapping, and one-site alignment gaps. Compare `XDiffPos(start - 1) + 1` and
  `XDiffPos(end + 1) - 1` aliases exactly.
- Combined cyclic rounds where GENECONV is strongest, ties each other active method, groups with
  cross-method support, disappears after tract erasure, or appears only through retained fragment
  re-entry. Confirm exact corrected/local-P ties use the supplied RDP → GENECONV → MaxChi →
  CHIMAERA → 3SEQ method-major order, then compare roles, correction factors, and all five GENECONV workload counters.
- Confirm the on-demand plot pairs tracks 0/5, 1/4, and 2/3 by colour; reports
  `-log10(raw KA P)`; preserves selected endpoints under thinning; and labels later-round
  reconstruction honestly.
- Confirm CSV/project JSON retain complete GENECONV anchor evidence and `v1alpha15` save/reload
  reproduces method code 3, trace values, counters, signal/event IDs, and review state. Pre-v13
  projects must restore GENECONV disabled even if stray new fields are present.
- For the event representatives and every finalized nonrepresentative row, compare the ordinary
  six-track recheck's profile/skew status, stable best fragment, provisional recombinant, raw and
  corrected P, overlap rejection count, and numerical fallbacks. Assert MaxChi preparation is
  shared, the recheck cannot move reconciled coordinates, and representative rows are skipped.
- Confirm the active automated minimum length/polymorphism/pair-score predicates remain reported
  inactive. Test permutation, manual-pair, alternative-indel, and full late event-reconstruction
  inputs only after those separately supplied workflows are ported; never infer their results from
  ordinary KA discovery or ordinary-kernel corroboration.

## 3SEQ discovery parity corpus

Use authorized desktop saved output from the supplied `FindSubSeqTS`/`FindSubSeqTS2` →
`Seq3PVals`/`GetTSPVal` → `CheckwrapC` → `TSXOver` route. Do not substitute a third-party 3SEQ
implementation as an oracle.

- Exercise all three target rotations with parent-one matches, parent-two matches, parent-equal,
  all-different, input missing, and earlier-erasure states at first/middle/last coordinates. Confirm
  the equality-slot maps `{0,2,1}` and `{1,0,2}`, source target order 3→1→2, and the effective
  four-information-rich-site floor.
- Build walks whose maximum descent wins, maximum ascent wins, and raw probabilities tie exactly.
  Confirm the strict lower-P swap, provisional major/minor mapping, counts, direction, excursion,
  and coordinate selection.
- Exercise `CheckwrapC` with `BE < EN`, `BE >= EN`, a new wrapped maximum, a new wrapped descent,
  beginning at the last informative site, and both circular/linear topology. Verify the source's
  post-wrap beginning advance and linear complement conversion independently. Assert that exact/
  fallback probability still uses the pre-wrap excursion while exported boundary evidence retains
  the potentially larger post-wrap value.
- Compare the compact exact hypergeometric DP with `Seq3PVals` table values for dense small
  `(m,n,k)` grids, including zero dimensions, `k=0`, `k=1`, `k=n`, unreachable excursions, exact
  float underflow, and values immediately on both sides of the project cutoff.
- Cross the eight-million-transition and 4096-dimension browser bounds. Compare
  `SiegmundDiscrete`, `ApproxNu`, normal CDF/PDF, the invalid-approximation scaled-table branch,
  power adjustment, and the source `10^-300` underflow behavior. Every exported call must label
  exact versus fallback honestly.
- Verify the post-orientation exits `nN > 0 && nK = 1` and `nN - nM = nK`, literal zero rejection,
  no-correction mode, and Dunn–Šidák correction for both exploratory triplet counts and
  query/reference group-pair × query opportunities. Test immediately above, at, and below
  `p = 10^-15`, plus subtraction underflow and `corrected == 1`/small-product calls where the
  source product comparison matters.
- In combined cyclic fixtures, make 3SEQ strongest, tie each active method, merge into an existing
  cross-method event, disappear after erasure, and appear through a retained fragment. Confirm the
  method-major tie order RDP → GENECONV → MaxChi → CHIMAERA → 3SEQ.
- Validate the active `FindSubSeqTS2` inclusive `XPosDiff`/`upper_bound` mapping around informative
  and missing coordinates. Test left/right missing splits, unchanged versus moved bounds,
  orientation re-selection, corrected-P re-gating, wrapped missing runs, and coordinate-zero/one
  aliases.
- Confirm the on-demand plot assigns target 0/1/2 to the three colour slots, retains signed extrema
  and both breakpoints after thinning, scales below and above zero, and labels a later-round
  original-alignment reconstruction.
- Confirm `v1alpha15` save/reload reproduces method code 4, all evidence, exact/fallback/split flags, four
  workload counters, signal/event IDs, and review state. Every pre-v14 project must restore 3SEQ
  disabled even if stray fields are present; unknown method labels must fail rather than become RDP.
- For event representatives and every finalized nonrepresentative row, compare `TSXOver(1)` target
  order, initial lower-P gate, both Findall orientations, forced split calls, the distinct correction
  branches, qualifying count, and the two `XOverList` records per accepted orientation. Exercise
  exact ties, one/two/no qualifying orientations, wrapped intervals, missing-state trimming, list
  capacity/replacement, and representative skipping. Assert the compact best-orientation summary
  cannot move reconciled coordinates or enter the authoritative discovery catalogue.
- Keep manual 100-permutation envelopes and full late event-catalogue reconstruction out of parity
  claims until those separately supplied paths are implemented.

## Primary BootScan discovery and cache

- Compare `SEQBOOT2` weights under seed 3 and several nondefault seeds, including replicate zero,
  both warm-up calls, and the unused final generated replicate that advances the next site's random
  stream. Compare `FastBootDist` float JC values and the saturated sentinel with gaps/ambiguities.
- Exercise strict `GetPltVal` ties; all three whole-sequence distance orderings; 40%, 50%, 70%, and
  exact-cutoff support; `OverlapNum` VB6 rounding; linear edges; circular padded windows; full-origin
  traversal; and equal beginning/ending rejection through `ScanBSPlots`, `FindBeginBS`, and
  `FindEndBS`.
- Compare `FindDaughter` provisional roles and both breakpoints before later `CentreBP` and BURT
  movement. Include a boundary on an invariant site to verify `XPosDiff`'s preceding-score behavior.
- Compare raw `MakeScoresBS`/`ProbCalc` probability for tract lengths 3, 169, 170, and longer;
  round-to-even match scaling; full-information-profile background probability; factorial versus
  log-domain accumulation; exact zero/underflow; and project correction. Confirm BURT does not
  recalculate the saved discovery p-value.
- Scan four or more sequences with repeated pairs and assert cache misses equal unique exact pairs,
  later triplets hit those profiles, peak bytes stay within 64 MiB, in-flight `shared_ptr` profiles
  survive FIFO eviction, and every erasure/re-entry round invalidates stale entries. Compare event
  results with a cache-disabled build.
- Confirm `v1alpha17` saves/restores method code 5, discovery evidence, all cache/work counters,
  event support, and review plot; pre-v17 imports must leave primary BootScan disabled. Exercise
  public-ABI and production-WASM Pages fixtures in both exploratory and query/reference schedules.
- Keep tree, similarity, alternative substitution-model, permutation, manual, and full late
  catalogue modes outside parity claims until their supplied paths are separately implemented.

## Workflow acceptance

- File data never leave the browser.
- The UI stays responsive during each scan round and cancellation completes between bounded batches.
- Review decisions are recorded in event order; later calls remain inspectable while blocked.
- Accepting an unchanged event advances directly. Correcting or rejecting one requires downstream
  reconciliation before later decisions or final alignment exports.
- A direct worker/API attempt to decide or edit a later event, request an unsolicited rebuild, or
  export any of the four event-derived final FASTA variants early is rejected even if UI controls
  are bypassed; the two curation-only FASTA partitions remain independently available.
- Accepted-event FASTA fixtures cover overlapping and wrapped tracts: remove current co-group rows,
  remove the inclusive union of event columns, mask group tracts, and split aligned fragments. Also
  verify the explicit all-rows-removed and all-columns-removed errors.
- Curation FASTA fixtures partition every original row exactly once between enabled-only and
  masked/disabled-only outputs, preserve alignment coordinates, and report an empty excluded-row
  selection without applying the accepted-event readiness gate. Exercise all three downloads before
  any scan exists, confirm full FASTA ignores curation, and repeat after project restore.
- Bulk curation fixtures reapply the supplied auto-mask, enable all, mask all, and disable all on a
  dataset larger than the 500 rendered-row cap; confirm disjoint state sets, primary triplet counts,
  blocked continuation below three enabled rows, and matching pre-scan FASTA partitions.
- A full-circle candidate run must not become a primary call when its mapped beginning and ending
  are equal (`FastRecCheckP` requires `EN != BE`); a manually entered/restored equal-endpoint
  circular event must still traverse the full alignment for erasure and accepted-event exports.
- Reloading a `v1alpha15` project reproduces analysis scheme, every reference-group assignment,
  settings, signal methods, MaxChi, CHIMAERA, GENECONV, and 3SEQ discovery traces,
  per-signal correction factors, cumulative triplets/rounds, all method workload counters, fragment
  provenance, event anchors, edits, decisions, and any pending invalidation marker. Assert that
  replaying saved events does not increase the authoritative cumulative count, and that the saved
  terminal reason plus final-round processed/total counts survive reload.
- Reloading during a pending correction restores only the valid event prefix, remaps every retained
  support/anchor signal ID, preserves manual group membership, and re-identifies—not replays—the tail.
- Loading `v1alpha1`–`v1alpha10` projects supplies exploratory query/reference defaults. Projects
  through `v1alpha9` keep MaxChi disabled, and every `v1alpha1`–`v1alpha11` project keeps CHIMAERA
  disabled; every pre-v13 project also keeps GENECONV disabled, and every pre-v14 project keeps
  3SEQ disabled, preserving saved discovery
  semantics. All older schemas deterministically rebuild the current evidence tier.
- Sequence-state fixtures distinguish enabled, masked, and disabled rows: only enabled rows affect
  primary triplet/correction counts; masked rows receive trace and late set evidence; disabled rows
  receive neither but remain present in saved tree leaves. Verify disjoint state replay and reject
  disabled manual event roles or co-recombinant membership.
- Breakpoint inspection returns recombinant/major/minor rows first, deduplicates group/evidence
  rows, never exceeds the requested row cap, and does not expose sequence sites outside the two
  requested windows.
- Linear breakpoint windows clip cleanly at sites 1 and `L`; circular windows wrap coordinates in
  order without repeating a coordinate, including alignments shorter than the requested context.
- Editing roles, breakpoints, or co-group membership refreshes an open inspector from the immutable
  original alignment rather than stale working-fragment rows.
- Profiles longer than 2,000 points retain the nearest samples for both breakpoints, the selected
  method peak, and each applicable maximum. CHIMAERA retains only its target/parent-one trace;
  GENECONV retains the three inner/outer KA fragment envelopes; 3SEQ retains all three signed
  target walks plus their extrema.
  First-round input-row signals are labelled exact; erased or
  fragment-assisted later-round signals are labelled original-alignment reconstructions.
- Event JSON, on-demand alignment JSON, CSV, and review badges agree on both uncertain boundaries,
  qualifying prior event IDs, immediate adjacency, the configured RDP window, and the nearest
  information-rich-site count; none labels the manual review interval as a statistical confidence
  interval.
- Validate active BURT confidence metadata against no/short information profiles, one and multiple
  HMM switches, inside/outside signed interval matches, linear clamps, circular wrapped 95%/99%
  ranges, input-missing repositioning, the three-usable-site revert, role-swap recomputation, and
  repeatability under the supplied seed. Preserve the DLL's `0 ... 20` inclusive training starts
  and strict `> 0.995`/`> 0.999` posterior gates.
- With breakpoint polishing disabled, verify that no event is marked attempted, primary/manual
  coordinates remain unchanged, project replay retains the setting, and evidence reports
  `unavailableReason: "disabled"`.
- Tree inspection returns all six saved regions, never transfers a distance matrix, distinguishes
  unusable-region JC fallback, and refreshes role/group labels after correction. Toggling weak-
  branch collapse changes only display length; it does not alter stored topology or evidence.
- JSON preserves every set, tree/fallback marker, role score/contribution, and complete group without
  claiming full native consensus parity.
- Tract-masked FASTA preserves sequence count/length and modifies only accepted current groups.
- Fragment FASTA preserves alignment length and reconstructs each processed sequence when its
  ordered fragments are overlaid on its remainder.
