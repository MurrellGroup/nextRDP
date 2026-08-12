# Validation plan

No compilation, bundling, type-checking, preview, or project runtime execution was performed in
sessions 1–9. Session 7 resolved lockfile metadata without installing dependencies or running
package scripts; sessions 8–9 used only source/PDF inspection and static text-contract checks. The
first runtime checkpoint must use a disposable, resource-limited environment and tiny inputs.

## Static/build gate for a later authorized session

1. Configure TypeScript with strict checking and resolve every diagnostic.
2. Compile the C++ target with Emscripten warnings enabled.
3. Confirm every generated export matches `rdp_api.h` and the worker module interface, especially
   the `v1alpha9` signal/event restore argument order, co-recombinant group pointers, and bounded
   event-alignment and event-tree JSON entry points.
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

## MaxChi confirmation parity corpus

Capture `LenXoverSeq`, `HWindowWidth`, `CriticalDiff`, all three pair profiles, banned positions,
`MaxX`/`MaxY`, initial and grown χ², `WinWin`, raw tail, within-triplet `xMPV`, and final corrected
probability from the supplied desktop/DLL path. Compare those values before allowing MaxChi evidence
to affect discovery, coordinates, group membership, or role voting.

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
- Representative and finalized-list rechecks with representative skip, unavailable profile, no
  critical peak, uncorrected trace, and corrected hit. Confirm no MaxChi value changes RDP event
  coordinates in this checkpoint.
- A future discovery fixture with multiple peaks must remain pending until smoothing,
  `DestroyPeakP`, retry order, and event construction are ported from the supplied source.

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
- Reloading a `v1alpha9` project reproduces settings, signals, per-signal correction factors,
  fragment provenance, event anchors, edits, decisions, and any pending invalidation marker.
- Reloading during a pending correction restores only the valid event prefix, remaps every retained
  support/anchor signal ID, preserves manual group membership, and re-identifies—not replays—the tail.
- Loading `v1alpha1`–`v1alpha8` projects supplies conservative defaults and deterministically rebuilds
  the current evidence tier.
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
