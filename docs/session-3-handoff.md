# Session 3 handoff

This is a source-only checkpoint. Nothing in the supplied desktop source or browser port was
compiled, bundled, served, previewed, or executed.

## Delivered in this phase

1. Every event now carries three role hypotheses. Each of the anchor triplet's sequences is treated
   in turn as the presumed recombinant and the other two are used as the hypothesis parents.
2. `MakeBPosLR`'s four circular boundary walks are represented with the supplied `VSN = 60` value.
   The region builder then follows the exact start/end arrangement passed to `MakeSDMP2`, including
   its asymmetric boundary inclusion.
3. Each candidate is reduced over five regions to the three `MakeSDMP2` pattern fractions: matches
   parent one, matches parent two, or differs where the two parents agree.
4. `FillRmat`'s three matrix pairs are reproduced without allocating general matrices: first
   breakpoint, second breakpoint, and mean outside profile versus the bounded tract.
5. Each pair is a fixed six-value Pearson correlation. With six observations, the two-sided t-test
   probability has four degrees of freedom and is evaluated directly as
   `1 - 1.5 × |r| + 0.5 × |r|³`.
6. Candidate eligibility records non-gap coverage across the two source boundary spans and applies
   the supplied `MakeGoodC` strict `> 10` rule. A positive direct-polarity correlation with
   uncorrected `P < 0.05` enters the distance-correlation set.
7. The currently mapped detectable-signal set and distance-correlation set are intersected. This is
   exposed as a `twoSetConsensus` subset, never as a complete co-recombinant set.
8. The review UI shows all role hypotheses and detailed evidence for the current one. Project JSON
   preserves coefficients, p-values, site counts, overlap eligibility, strongest matrix pair, and
   detectable support; CSV carries the principal set memberships.
9. Project schema `v1alpha3` is emitted. Imports from `v1alpha1`, `v1alpha2`, and `v1alpha3` restore
   the saved alignment/signals and deterministically recompute the current evidence.

## Fidelity boundary

The manual-level direct distance-correlation rule is present, but the entire native secondary
stack is not. In particular:

- `CalCR`'s inverse category permutations, `RCorrWarn`, `AcceptableCoR`, and `MakeRList` aggregate
  filters remain explicit deltas.
- Detectable-set collection is anchored to the two hypothesis parents and does not yet reproduce
  every later native expansion/rescoring pass.
- Six Jukes–Cantor distance matrices, bootstrapped neighbor-joining trees, support collapse, and the
  phylogenetic-correlation set are absent.
- Because the third evidence set is absent, the general “present in any two of three” grouping is
  not complete. The current intersection is conservative but may omit sequences that would be
  supported by detectable + phylogenetic or distance + phylogenetic evidence.
- The full recombinant/parent identification score battery is absent, so roles remain provisional.

## Next phase

Port the supplied phylogenetic-correlation path in layers: source-faithful Jukes–Cantor matrices for
the six existing coordinate regions; deterministic bootstrap sampling; neighbor joining and
sub-50% support collapse; tree-side grouping; then the complete two-of-three set combination. Keep
the matrix and tree work in the analysis worker and reuse compact triangular distances so the UI
thread and memory profile remain suitable for static browser hosting.

Before claiming parity, add golden fixtures that compare boundary coordinates, the fifteen pattern
fractions per candidate, all nine role/pair correlations, native inversion/warning decisions, and
the three evidence sets against the supplied application.
