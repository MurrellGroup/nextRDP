# Supplied-source SISCAN discovery and confirmation trace

This note records the source path used for the browser port. It was traced only from the supplied
RDP5 Visual Basic project, `dnaDLLSource`, `dna5DLLSource`, and manual. No alternate RDP implementation was consulted.

## Native entry points and defaults

The ordinary scanner calls `Module3b.bas::SSXoverC`. The equivalent longer implementation remains
in `Module3.bas`, while `Module30.bas` contains the manual/plotting path. Default setup in
`Module5.bas` is:

- `SSFastFlag = 1`;
- `SSGapFlag = 0` (strip gaps);
- `SSVarPFlag = 2` (one-, two-, and three-variable categories);
- `SSOutlyerFlag = 1` (nearest fourth sequence);
- `SSWinLen = 200`, `SSStep = 20`;
- `SSNumPerms = 1000` for the final region and `SSNumPerms2 = 100` for scanning;
- random seed 3 in the standard configuration.

The desktop program leaves primary SISCAN discovery off in the normal automated selection but
uses it as an ordinary later confirmation method. The browser therefore defaults primary SISCAN
off and fixed-region confirmation on; both are explicit settings.

## Call mapping

| Supplied source | Browser implementation | Retained behavior |
| --- | --- | --- |
| `Module3b.bas::SSXoverC` | `siscan_discover` | Window screen, pair-switch region construction, boundary shrink, final region scoring, role geometry, and threshold gate |
| `threshold.CPP::MakeDistanceBakB` plus RDP tree preparation and `threshold.CPP::GetSSOL` | `build_source_wpgma_context` and `nearest_source_outlier` | Direct similarity matrix, source-shaped WPGMA/cophenetic context, and nearest eligible fourth sequence |
| `threshold.CPP::Get3Score` | `build_triplet_categories` | Gap-stripped triplet pattern categories |
| `threshold.CPP::GetPScores2` | `pattern_category` and `count_patterns` | Fourth-sequence conversion into the exact 16 partition categories |
| `Module2.bas::SetUpSiScan` | compile-time category tables in `siscan.cpp` | Exact `Seq34Conv`, partition groups, and summed groups used by default variable-pattern mode |
| `threshold.CPP::MakeVRand` | `ensure_vertical_random_prefix` | Microsoft CRT `rand()` stream and flat byte-prefix template |
| `threshold.CPP::DoPerms3` / `dna5 MathFuncsDll.cpp::DoPerms3P` | `permute_patterns` | Category-constrained vertical randomization with the scan/final permutation counts |
| `threshold.CPP::MakeZValue2` | `source_z_score` | Population variance and zero-variance behavior for partition and summed scores |
| `threshold.CPP::DoSums` | summed-score construction in `permute_patterns` | Exact 15 summed categories from the 16 partition counts |
| `threshold.CPP::QuickCheckB` | `quick_check_window` | Fast-screen control flow, including the supplied missing-braces quirk |
| `threshold.CPP::FindMaxZ` | `strongest_window_pair` / `strongest_region_score` | Highest-pair comparison and eligible alternative partition/summed score |
| `threshold.CPP::ShrinkRegionC` | `shrink_region` | Conversion from consecutive switched windows to nucleotide bounds |
| `Module30.bas::NormalZ` | `source_normal_z` | Supplied polynomial normal-tail approximation |
| late `SSXoverC` calls from `FinalTrim` | `siscan_recheck` | Nearest-outlier full fixed-region score without changing reconciled coordinates |

## Probability stages

The browser preserves the desktop stages as separate evidence rather than collapsing them:

1. `normalTailPValue = NormalZ(maxZ)`;
2. `regionLengthAdjustedPValue = normalTailPValue × alignmentLength / regionLength`;
3. `windowAdjustedPValue = regionLengthAdjustedPValue × alignmentLength / SSWinLen`;
4. when project correction is active, multiply by the frozen initial scan-plan opportunity count;
5. clamp the reported probability to `[10^-300, 1]` so a finite nonzero audit value survives
   ordinary floating-point underflow.

This separation is important when comparing very small values. The port evaluates the polynomial
and products in `double` and retains a positive subnormal where possible. A desktop build may round
or underflow at a different stage; such a numerical difference is not treated as evidence that the
breakpoints changed.

## Caches and cyclic invalidation

The desktop supplies the distance/tree context and vertical-randomization template to repeated
`SSXoverC` calls. The browser keeps the same lifetimes explicitly:

- the direct and WPGMA/cophenetic matrices are built once per cyclic round and reused by every
  SISCAN triplet;
- the seeded MakeVRand byte stream is an extensible flat prefix, so repeated windows and later
  rounds reuse identical bytes rather than regenerating them;
- accepting an event invalidates only the state-dependent distance/tree context;
- starting another analysis in the same browser context resets cache counters while retaining the
  deterministic prefix, so progress never includes telemetry from the preceding dataset;
- `v1alpha19` project restoration round-trips the WPGMA build/comparison/merge and generated-random
  counts as audit data; cached matrices themselves are deliberately reconstructed when work resumes;
- the existing XOverList/BestXOList-style triplet shortlist skips SISCAN entirely for an unchanged
  triplet whose method summary can be replayed.

This keeps the expensive context at `O(N²L) + O(N² log N)` per affected round rather than per
triplet. Per-window work remains proportional to retained category sites times the requested
permutations.

## Browser safety normalization

RDP5 can re-enter fragments created by previous erasure. The browser excludes any proposed fourth
sequence whose original-sequence identity matches one of the three triplet members. This is the
`TraceSub` intent made explicit and prevents a fragment from serving as its own outlier. Disabled
original sequences are also ineligible. The source scans linear SISCAN windows even for circular
datasets; circularity is retained for event geometry and wrapped fixed-region confirmation, not by
inventing an extra origin-crossing scan window.

The desktop manual plot exposes all 15 partition categories and 9 category sums. The ordinary
browser review endpoint is intentionally smaller: for each of the three sister pairs and each
window it selects the eligible partition/summed value with greatest absolute Z and preserves that
value's sign. This provides a signed diagnostic without claiming full manual-plot parity.

## Validation boundary

`scripts/verify-siscan-core.cpp` deterministically plants a sister-pair switch and checks outlier
selection, roles, region coverage, probability stages, cached WPGMA reuse, bit-stable random-prefix
replay (including the first 20 Microsoft-CRT seed-3 bytes), fixed-region confirmation, review
plotting, round invalidation, and same-context analysis restart. It is a port regression,
not native parity evidence. Authorized native saved-output fixtures are still required for exact
window maps, category Z scores, final bounds, and every probability stage.
