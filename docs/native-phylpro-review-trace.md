# Supplied PHYLPRO event-review trace

This trace covers the ordinary graphical PHYLPRO route used while reviewing a selected event. It
was derived only from the supplied RDP5 manual, `Module30.bas`, and `threshold.CPP`.
No alternate RDP implementation was consulted, and the supplied desktop/native sources were not compiled.

## Workflow boundary

The manual presents PHYLPRO as a recombinant-identification and breakpoint-review plot. It also
states that its proposed permutation test was not implemented, so PHYLPRO cannot provide the fast
significance calculation needed for exploratory automated discovery. The browser therefore creates
this profile only when a reviewer opens it for a reconciled event. It does not emit a signal, alter
breakpoints, affect event order, or invent a p-value.

The default review request uses the supplied ordinary settings: total window 60, pairwise missing
observations ignored, and the zero-distance self observation excluded. The reviewer may change
those three plot inputs without rerunning discovery. The project does not serialize the generated
points; they are deterministic and are rebuilt from the saved alignment, current event roles and
request options.

## Supplied call path

| Supplied routine | Browser mapping | Retained behavior |
| --- | --- | --- |
| `FindSubSeqPP` | `eligible_columns` | Build an alignment-coordinate map after the selected gap policy, retaining polymorphic columns |
| `PXoverD` | `phylpro_profile` | Form equal left/right half-windows around each mapped partition and evaluate the three displayed role sequences |
| `MakePDstMat` | initial `TargetWindows` population | Count pairwise Hamming differences, ignoring a comparison whenever either state is missing |
| `UpdatePDstMat` | four rolling `add_position` updates | Remove the outgoing left site, move the central site from right to left, and add the incoming right site |
| `PPRegression` | `source_pearson` | Correlate the two distance vectors with the supplied double-precision sums and single-precision output; return 1 for zero variance |
| `PXoverD(CFlag = 0)` | `event_phylpro_json` / `EventPhylproInspector` | Return and display recombinant, major-parent and minor-parent curves only |

Coordinates in the API are one-based alignment positions. Circular scans emit one value per
eligible mapped column and reproduce the supplied wraparound half-windows. For a linear analysis,
the browser deliberately emits only partitions with two complete half-windows instead of importing
the opposite alignment edge. An oversized request is capped to half of the eligible profile.

## Bounded target-row optimization

The supplied DLL allocates and rolls both halves of an all-pairs `N × N` distance matrix even when
the ordinary event plot reads only three rows. A Pearson coefficient for a selected row depends on
that target-to-context row and not on distances between two other context sequences. The browser
therefore stores only the recombinant, major-parent, and minor-parent rows: a three-target
optimization.

For `L` eligible columns and `N` context sequences, this changes the profile from `O(LN²)` rolling
matrix work and `O(N²)` matrix memory to `O(LN)` work and `O(N)` working memory. The initial and
rolling counts, observation order, self exclusion, sums, variance guards, and final float narrowing
are otherwise unchanged; the final `AA / BB` value is not clamped. `verify-phylpro-core.cpp`
compares every optimized point against an
independent full-window recomputation across circular/linear topology, both gap policies, self
inclusion, disabled context, and window capping.

The JSON response is downsampled only after the complete profile has been evaluated. It always
retains both ends, all three global minima, and samples nearest both event breakpoints; calculation
minima and telemetry continue to describe the full profile.

## Explicit source repair and adaptations

The commented per-column `ReDim Match(255)` is not treated as a defect. RDP's `SeqNum` encoding
stores A/C/G/T as 66/68/72/85 (the character code plus one) and missing data as 46; `FindSubSeqPP`
explicitly clears exactly those five counters on each column. The browser's fresh four-state set is
the compact normalized equivalent of that reset.

One supplied `FindSubSeqPP` indexing detail is treated as a defect rather than copied silently:

- after compacting eligible sequences, the supplied active path fills both `SubMaskSeq(compact)`
  and `RevSeq(original)` but indexes `RevSeq(compact)` in its scan/DLL call. The browser uses an
  explicit compact-to-original context vector, so a disabled row in the middle cannot substitute
  a different sequence or duplicate row zero.

Masked original sequences remain PHYLPRO context, while disabled originals are excluded. Synthetic
cyclic fragments are not inserted into this review plot: roles and context are reconstructed from
the immutable original alignment. The manual's conceptual prose describes all alignment columns,
whereas the supplied active `FindSubSeqPP` calculation constructs a polymorphic-column map. The
browser follows that active intended calculation and labels the difference in the UI/API.

As throughout this port, the compact browser alignment normalizes IUPAC ambiguity and gap-like
symbols to missing state zero and U to T. PHYLPRO then ignores a pair containing state zero, or
removes its whole column under the strip policy. The literal DLL snippets test encoded gap 46 only
and can compare other raw character codes as categorical states. Ambiguity-heavy PHYLPRO output is
therefore an explicit native-golden boundary rather than a parity claim.

This compact-index repair can change a profile relative to execution of the defective statement in
this particular source snapshot. It is not a statistical improvement and does not create
significance; it makes the documented distance-profile context deterministic and auditable. Authorized
saved desktop curves are still required before calling the route native-parity validated.

## Deliberate non-claims

- No PHYLPRO permutation p-value is implemented, because the supplied ordinary RDP5 workflow has
  no active significance test.
- PHYLPRO is not an eighth exploratory discovery method and does not enter the cyclic shortlist.
- The browser does not claim native validation for equal-value rendering, obscure ambiguity codes,
  the compact-index defect above, or manual/printing-only plot variants.
- No alternate reference implementation was used.
