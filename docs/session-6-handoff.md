# Session 6 handoff

This checkpoint is source-only. No compiler, Emscripten build, npm build, type-checker, preview,
server, or runtime test was invoked.

## Completed in this phase

1. Traced the active `ConsensusOK` and remaining `FinalTrim` source paths. The native final score
   consumes pattern, tree-position, breakpoint-distance, method-recheck, and `CalcMatchY` families.
   This phase now maps the base through `OKSeq 6`; final-list membership 15 is still unavailable, so
   the browser continues to prohibit the partial late score from pruning a co-recombinant group.
2. Corrected the opening `FinalTrim` duplicate-correlation implementation. The count still uses
   warning/inversion-cleared direct correlations strictly above `0.83`, but once a sequence/pair is
   duplicated its copied value is now cleared at every occurrence in all three RLists, exactly as
   the source's separate counting and clearing loops require.
3. Added `RdpScanner::event_alignment_json` and `_rdp_get_event_alignment_json`. It reads only the
   immutable original alignment and emits two bounded windows centred on the current beginning and
   ending breakpoints.
4. Implemented row priority matching the manual's “show relevant sequences” review intent:
   recombinant, major parent, minor parent, manually current group, automatic group, masked trace,
   detectable set, distance-correlation set, and phylogenetic-correlation set. Rows are unique,
   representatives cannot be displaced by the cap, and omitted-row counts remain explicit.
5. Added circular coordinate wrapping and linear endpoint clipping. Requested context is clamped to
   5–100 sites per side, short circular alignments avoid repeated coordinates, and rows are clamped
   to 3–64.
6. Added a lazy worker/client request. The browser transfers no breakpoint sequence data until the
   inspector is opened, and never receives the complete alignment through this path.
7. Added a responsive React alignment inspector with both breakpoint panels, 1-based coordinate
   rulers, selectable ±15/30/60 context, sticky sequence labels, role/group/trace badges, and
   nucleotide colours for major-parent matches, minor-parent matches, shared-parent states,
   recombinant-only states, and gaps/ambiguities.
8. An open inspector refreshes after role, breakpoint, or co-group edits. Selecting another event
   closes and clears the previous alignment slice so stale sequence context is not shown.
9. Added persistent breakpoint-erasure context. If either current boundary or its immediate
   neighbour falls inside a prior erased tract affecting one of the three representatives, the
   event JSON, CSV, review header, and inspector identify the earlier event(s) and tell the reviewer
   to treat that boundary as uncertain, matching the manual's deleted-region warning.
10. The inspector also finds the closest expected parent-informative state on each side—major to
    minor at the beginning, minor to major at the end—and displays the resulting manual review
    bracket explicitly as something distinct from a statistical confidence interval.
11. Added result-level `lateConsensus` and `breakpointInspection` metadata. It states that late
   consensus is diagnostic, lists implemented and pending stages, and records the alignment bounds.
12. Advanced the engine to `0.6.0-session-6` and project JSON to
    `org.rdp-web.project/v1alpha7`; import accepts `v1alpha1` through `v1alpha7`.
13. Retained compact topology edge lists, branch lengths, bootstrap support/collapse flags, and
    original/working/fragment provenance for the six NJ trees already built during event
    reconciliation. Ordinary results and project JSON still retain summaries only.
14. Added `RdpScanner::event_trees_json`, `_rdp_get_event_trees_json`, and a lazy worker/client path.
    The request transfers shared leaf metadata plus six `O(K)` edge lists; it never transfers the
    `O(K²)` JC or patristic matrices and never rebuilds a tree during review.
15. Added a responsive paired tree inspector for whole-event, 5′, and 3′ comparisons. It labels
    current roles/groups, masked traces, retained fragments and bootstrap percentages; weak branch
    collapse is reversible, and arbitrary rooting is explicitly a display choice.
16. Added result-level `treeInspection` capability metadata and static validation cases for
    topology/leaf/provenance correspondence.
17. Ported the exact score branches for `FinalTrim` `OKSeq` 7, 8, 9, 12, and 13 using the already
    available collapsed/raw outside/inside patristic matrices and four breakpoint JC matrices.
    This includes the source's asymmetric tied rewards, negative penalties, warning and saturation
    gates, repeated-closest-pair role modifier, and distinct-pair positive downweighting.
18. Added per-candidate matrix scores, availability/fallback flags, and the raw native matrix
    subtotal to JSON, CSV, and ordered review. It is explicitly diagnostic and cannot prune or add
    group membership while final `OKSeq 15` and `ConsensusOK` list rebuilding remain missing.
19. Ported `CalcMatchY` from the supplied `Module3.bas`, `MakeVarMap2`, `MakeCntHit2`, and
    `MakeLenFrag` paths. The implementation preserves the four 40-variable-site scans, the
    three-alignment-length and 160-variable-site bounds, VB `CLng` half-to-even rounding, signed
    map, single-precision circular smoother, regional product (`OKSeq 17`), and six-checkpoint
    standard breakpoint class (`OKSeq 18`).
20. Added the first active `ConsensusOK` pass over `OKSeq 18`: raw outside/inside patristic topology
    must keep the candidate closest to its proposed role and preserve the representative ordering.
    Raw and topology-filtered classes, breakpoint-presence flags, fallback state, and all checkpoint
    values are retained in JSON, CSV, and the review table without changing group membership.
21. Replaced cyclic parent-role arithmetic with the source `CompMat` ordering. This fixes role 1's
    parent order across correlation, matrix, topology, and displayed hypothesis fields while keeping
    intentionally symmetric vote paths unchanged.
22. Traced the active VB assignments for `OKSeq` 10 and 11. Element 10 is explicitly zeroed for all
    roles/candidates and 11 has no active assignment, so both are now recorded as source-zero slots.
23. Ported the active DLL `FindActualEvents` lookup and `MakeMatchMatX2P` interval distances used by
    `OKSeq 14`. The port retains pre-`StripDupInv` membership, inverse-parent rejection, strict
    greater-than-one-third symmetric overlap, best-event selection, circular/equal-endpoint walks,
    tract intersection, 30-site/0.75 saturation, the bare-`CompMat` sequence-index quirk, and all
    signed score branches. Per-pair valid/difference prefixes avoid nine repeated alignment scans.
24. Completed the active 7–14 matrix subtotal and exposed the seven `MakeMatchMat` inputs, selected
    event tract/overlap, saturation, raw subtotal, and source-index disclosure in JSON, CSV, and UI.
25. Ported `CheckPatternX`'s three region-by-role shares for `OKSeq 3`, plus `ConsensusOK` inputs
    0–5 and `RCorrX`. Corrected correlation P values, event overlap, set/list/duplicate membership,
    informative pattern shares, representative sentinel, and other-representative zero cells are
    individually auditable.
26. Ported the `FinalTrim` `OKSeq 6` nearest-nonrecombinant fixed-point pass. It reconstructs the
    native post-`StripDupInv` swap-last order, raw/collapsed outside/inside limits, strict correlation
    tickets, event-found gate, iterative deletion, and paired four-breakpoint JC veto.
27. Added a non-pruning CScore preview through `OKSeq` 0–6, `RCorrX`, and the complete active matrix
    family. The port preserves the source's `NS As Long` narrowing after every addition—including
    its zero-valued sub-one fallback—and labels `OKSeq 15` as the sole neutral membership input.

## Fidelity boundary

- Primary RDP cyclic discovery, bounded fragment re-entry, all three early evidence sets, manual
  group/role/breakpoint correction, ordered downstream rebuild, project replay, and accepted-event
  alignment variants remain implemented but uncompiled and unvalidated.
- The alignment inspector is a faithful browser representation of the manual's relevant-sequence
  review workflow, not a port of the desktop drawing code. Its underlying bases and coordinates are
  exact slices of the original input alignment; its colouring is a transparent UI derivation from
  the three current representative rows.
- The graphical tree inspector is likewise a browser representation rather than a port of the
  desktop drawing routines. Its topology, branch lengths, bootstrap support, weak-branch state,
  and leaf order come directly from the event trees used by the implemented reconciliation path;
  its chosen root, ordering, colours, and SVG geometry are display-only.
- The active `FinalTrim` matrix family 7–14, nearest-nonrecombinant membership 6, `CalcMatchY`
  17/18, `CheckPatternX`, and the `ConsensusOK` score path through 0–6 are mapped. Post-6 list
  expansion, selected-role special pruning, final membership 15, and the `ConsensusOK` list
  rebuild/straggler/score-and-prune stages remain pending. None of these diagnostics deletes the
  browser's current two-of-three membership yet.
- Native breakpoint uncertainty beyond the erased-tract adjacency warning remains open.

## Static inspection performed

Only source text was inspected. The ABI/header/export/worker names, project-schema acceptance list,
JSON field names, and documentation were reviewed without executing project code. See
`docs/validation-plan.md` for the deliberately deferred build and runtime gates.

## Next fidelity phase

Continue immediately after the mapped OKSeq 6 pass: port both distance-only RList expansions, the
order-sensitive selected-role special pruning, and final membership OKSeq 15. Then reproduce the
`ConsensusOK` score/`CalcMatchY` acceptance branches, equivalence and straggler expansion, and
empty-role rollback before allowing the late stack to change automatic membership. Breakpoint
uncertainty ranges remain a separate review-fidelity task.
