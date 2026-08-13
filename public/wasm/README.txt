The generated Emscripten module is written here by scripts/build-wasm.sh.

Session 20 adds the supplied event-tree calculation path and keeps the
supplied native projects read-only.
The source exports lazy breakpoint-alignment and six-region event-tree review endpoints.
Six-region event trees now use source-shaped single-precision Clearcut NJ,
Microsoft-CRT SEQBOOT2 resampling, retained-base-tree support pseudocounts,
50% node collapse, and the midpoint-rooted ultrametric rank matrices consumed by
TreePhPr, TreeSubPhPr, TreeSubDist, TrpScore, and reconciliation. Graphs retain
the source-serialized five-decimal branch lengths while analytical comparisons
use those ranked matrices.
The primary-RDP late-consensus source now actively maps OKSeq 0-18, FinalTrim,
ConsensusOK, shared selected-tree cleanup, and finalized-list RDP rechecks.
Source-shaped ordinary-triplet MaxChi discovery is active alongside RDP, including
the grow/side/optimized-tract and smoothed multi-peak destruction/retry lifecycle.
Source-shaped ordinary-triplet CHIMAERA discovery is also active: all three target
rotations build their information-rich binary profiles, scan raw chi-square peaks,
grow and optimize tracts, and enter the same strongest-first cyclic scheduler. Its
profile preparation reuses the cached MaxChi variable-site alignment-byte pass.
Source-shaped ordinary-triplet GENECONV discovery is active with six signed
inner/outer fragment tracks, ignored indels, mismatch penalties, lambda/K and
Karlin-Altschul probabilities, stable lowest-P overlap selection, restore, review,
plot, and export coverage. It reuses that same cached equality pass; prefix/range
queries replace quadratic fragment extension and a bounded numerical fallback
prevents the supplied Newton start from hanging a worker.
Source-shaped primary BootScan discovery is active with seeded SEQBOOT2 weights,
strict closest-pair votes, MakeScoresBS binomial probabilities, and a bounded
64 MiB shared pair/window/bootstrap profile cache. The cache avoids recalculating
pair work across triplets in one round and is invalidated after cyclic erasure.
Fully exploratory and automated one-query/two-cross-group-reference schedules are
active; the latter uses a lazy catalogue and preserves the supplied group-pair x
query correction separately from exact record-triplet progress.
The separate BootScan, MaxChi, three-target CHIMAERA, ordinary six-track GENECONV,
and two-orientation 3SEQ event/list corroboration paths remain active and share one
prepared triplet pass without changing reconciled coordinates.
CheckEnds uncertainty and BURT/BenHMM statistical confidence are active in source.
The linked host regression covers primary BootScan discovery, cache reuse,
cyclic reconciliation, and plot output. Authorized native golden validation remains pending.
