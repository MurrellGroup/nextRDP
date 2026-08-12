# Port status — session 5

This status is deliberately conservative. “Implemented” means source exists; it does not mean
compiled or native-parity validated in this checkpoint.

| Area | State | Notes |
| --- | --- | --- |
| Manual workflow mapping | Implemented | Dataset, settings, cyclic scan, ordered review/repair, and accepted-event exports are represented |
| Local alignment loading | Implemented | FASTA, GDE, CLUSTAL/MUSCLE, sequential/interleaved PHYLIP, NEXUS, and MEGA |
| Dataset diagnostics and auto-mask | Implemented | Diversity, missing data, variable/informative sites, warnings, and the supplied closest-pair masking workflow |
| Primary RDP triplet screen | Implemented, unvalidated | C++20/WASM source; information-rich windows, tract boundaries, probabilities, and correction cap |
| Non-blocking scan control | Implemented | Dedicated worker, 512-combination batches, round-aware progress, and cancellation |
| Cyclic event discovery | Implemented, unvalidated | Strongest unexplained event is reconciled and erased before a fresh complete eligible-triplet pass |
| Source-shaped fragment re-entry | Implemented, bounded | Gap-padded working fragments retain original/event provenance; same-origin triplets are forbidden; 100,000-site source cutoff and 256-fragment browser cap are explicit |
| Detectable-signal reconciliation | Implemented, unvalidated | Two shared original sequences, strict >30% symmetric tract overlap, and iterative three-role closure |
| Masked-sequence trace checks | Implemented, unvalidated | Relaxed follow-up RDP profiles; corrected significance is retained separately |
| Distance-pattern correlation set | Implemented, unvalidated | Five regions, six-value Pearson, all five relabellings, `RCorrWarn`, `MakeGoodC`, active `MakeACOR`, `MakeRList` dual-r override/aggregate, `StripDupInv`, and first `FinalTrim` duplicate-correlation cleanup |
| Event phylogenetic matrices | Implemented, unvalidated | Six Jukes–Cantor matrices with compact direct/reference caching |
| Bootstrap neighbour joining | Implemented, unvalidated | Deterministic ten-replicate event path; internal branches below 50% support collapse to zero length |
| Phylogenetic-correlation set | Implemented, unvalidated | Paired collapsed/raw tree affinity; labelled distance fallback outside the bounded tree panel |
| Complete two-of-three group | Implemented, unvalidated | Union of sequences present in any two of detectable, distance, and phylogenetic sets |
| Recombinant/parent recommendation | Implemented, partial | Nine auditable metrics; eight use mapped native full/half weights, while unported native method families keep `nativeWeightParity` false |
| Event review and correction | Implemented | Ordered accept/reject, role recommendation, manual roles/breakpoints, manual complete-group membership, and explicit downstream invalidation |
| Downstream re-identification | Implemented, unvalidated | Corrected accepted events are re-erased; rejected tracts are restored; fixed earlier calls remain auditable while all later calls are rediscovered |
| Project/CSV export | Implemented | Reloadable project schema `v1alpha6`; round correction factors, fragment provenance, all evidence sets, role metrics, traces, and decisions are retained |
| Project import/resume | Implemented, unvalidated | Schemas `v1alpha1`–`v1alpha6` restore primary state; pending repairs retain/remap only the valid signal/event prefix before downstream discovery resumes |
| Recombination-free FASTA | Implemented, unvalidated | Accepted tracts are replaced by gaps for current complete co-recombinant groups; core and UI reject incomplete-review export |
| Mosaic-fragment FASTA | Implemented, unvalidated | Event-ordered tract-masked originals plus aligned fragment-only records; core and UI reject incomplete-review export |
| Exact late native consensus | Partial | Active `MakeACOR`/`MakeRList`/`StripDupInv` and initial `FinalTrim` duplicate cleanup are mapped; complete `ConsensusOK` scoring and later tree/distance pruning remain |
| GENECONV / MAXCHI and other methods | Not yet ported | Queued after primary RDP native parity |
| Graphical tree/alignment review | Not yet ported | Current UI exposes numeric region, bootstrap, membership, and profile diagnostics |
| Compilation and runtime testing | Not performed | Explicit project constraint for this checkpoint |
| Native-vs-WASM golden parity suite | Designed only | See `docs/validation-plan.md` |

## Recommended next checkpoint

1. Port the remaining evidence fields and threshold branches needed by `ConsensusOK`, then add them to JSON as
   separately auditable late-consensus diagnostics before allowing them to prune a set.
2. Continue the active `FinalTrim` path after duplicate positive-correlation cleanup: nearest
   non-recombinant tree boundaries and breakpoint-relative distance checks.
3. Add breakpoint uncertainty/range diagnostics for events adjacent to previously erased tracts.
4. Add the remaining recombinant-identification method families without hiding the current
   per-method scores or native full/half contributions.
5. Only after explicit permission to execute code, establish small golden primary/secondary results
   against the supplied desktop application before calling the port parity validated.
