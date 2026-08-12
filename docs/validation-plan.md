# Validation plan

No compilation, bundling, type-checking, preview, or runtime execution was performed in sessions
1–5. The first runtime checkpoint must use a disposable, resource-limited environment and tiny inputs.

## Static/build gate for a later authorized session

1. Configure TypeScript with strict checking and resolve every diagnostic.
2. Compile the C++ target with Emscripten warnings enabled.
3. Confirm every generated export matches `rdp_api.h` and the worker module interface, especially
   the `v1alpha6` signal/event restore argument order and co-recombinant group pointers.
4. Check the production bundle contains relative URLs and both `.mjs` and `.wasm` assets.
5. Load the smallest three-sequence fixture before any multi-event or production input.

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
- Endpoint fixtures that distinguish inclusive `ModSN` erasure from the port’s current
  endpoint-exclusive tract helper.
- Even and odd window settings.
- Exactly 29, 30, and 31 information-rich sites around the default window boundary.
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
  lists; compare the opening `FinalTrim` duplicate counts and filtered pair flags.
- Six-region JC matrices with gaps, saturation, exactly 9/10 comparable positions, and both whole
  tract partitions.
- Four-to-ten-taxon NJ ties, negative limbs, and known splits; compare native patristic matrices.
- Bootstrap branches at 4/10, 5/10, and 6/10 support.
- A candidate satisfying each paired-tree affinity check, plus one that passes raw trees but fails
  after collapse.
- Active tree panels of 99, 100, and 101 records to verify full-panel/fallback labeling.
- Every detectable/distance/phylogenetic truth-table combination for the two-of-three group.
- Role fixtures where each displayed metric wins a different role; compare scores, full/half
  contributions, the collapsed-tree special condition, mapped weights, and the subset winner
  separately from the full desktop winner.

The primary RDP milestone is complete only when every mismatch is fixed or documented as an
approved representation/performance change.

## Workflow acceptance

- File data never leave the browser.
- The UI stays responsive during each scan round and cancellation completes between bounded batches.
- Review decisions are recorded in event order; later calls remain inspectable while blocked.
- Accepting an unchanged event advances directly. Correcting or rejecting one requires downstream
  reconciliation before later decisions or final alignment exports.
- A direct worker/API attempt to decide or edit a later event, request an unsolicited rebuild, or
  export either final FASTA early is rejected even if UI controls are bypassed.
- Reloading a `v1alpha6` project reproduces settings, signals, per-signal correction factors,
  fragment provenance, event anchors, edits, decisions, and any pending invalidation marker.
- Reloading during a pending correction restores only the valid event prefix, remaps every retained
  support/anchor signal ID, preserves manual group membership, and re-identifies—not replays—the tail.
- Loading `v1alpha1`–`v1alpha5` projects supplies conservative defaults and deterministically rebuilds
  the current evidence tier.
- JSON preserves every set, tree/fallback marker, role score/contribution, and complete group without
  claiming full native consensus parity.
- Tract-masked FASTA preserves sequence count/length and modifies only accepted current groups.
- Fragment FASTA preserves alignment length and reconstructs each processed sequence when its
  ordered fragments are overlaid on its remainder.
