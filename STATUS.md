# Port status — session 9

This status is deliberately conservative. “Implemented” means source exists; it does not mean
compiled or native-parity validated in this checkpoint.

| Area | State | Notes |
| --- | --- | --- |
| Manual workflow mapping | Implemented | Dataset, settings, cyclic scan, ordered review/repair, and accepted-event exports are represented |
| Local alignment loading | Implemented | FASTA, GDE, CLUSTAL/MUSCLE, sequential/interleaved PHYLIP, NEXUS, and MEGA |
| Dataset diagnostics and sequence curation | Implemented, unvalidated | Diversity, missing data, variable/informative sites, warnings, supplied closest-pair auto-mask, separate enabled/masked/disabled states, and auto-mask/enable-all/mask-all/disable-all controls that operate beyond the 500-row rendering cap; masked rows re-enter secondary evidence and all rows remain tree context |
| Primary RDP triplet screen | Implemented, unvalidated | C++20/WASM source; information-rich windows, tract boundaries, probabilities, and correction cap |
| Non-blocking scan control | Implemented | Dedicated worker, 512-combination batches, round-aware progress, and cancellation |
| Cyclic event discovery | Implemented, unvalidated | Strongest unexplained event is reconciled and erased before a fresh complete eligible-triplet pass |
| Source-shaped fragment re-entry | Implemented, bounded | Gap-padded working fragments retain original/event provenance; same-origin triplets are forbidden; 100,000-site source cutoff and 256-fragment browser cap are explicit |
| Detectable-signal reconciliation | Implemented, unvalidated | Two shared original sequences, strict >30% symmetric tract overlap, and iterative three-role closure |
| Masked-sequence trace checks | Implemented, unvalidated | Relaxed follow-up RDP profiles; corrected significance is retained separately |
| Distance-pattern correlation set | Implemented, unvalidated | Five regions, six-value Pearson, all five relabellings, `RCorrWarn`, `MakeGoodC`, `MakeACOR`, `MakeRList`, `StripDupInv`, full active RFF=0 `FinalTrim`, and `ConsensusOK` list rebuilding |
| Nearest-nonrecombinant membership | Implemented, active | `FinalTrim` `OKSeq 6` fixed-point pass preserves post-inversion list order, native thresholds, event-found gate, raw/collapsed tree bounds, and paired breakpoint-distance veto before both final expansions |
| Late matrix scoring | Implemented, active | `FinalTrim` `OKSeq` 7–14 supplies collapsed/raw tree position, whole-tract/breakpoint JC, source-zero 10/11, and detected-event distance to completed CScore and active grouping |
| Regional/breakpoint matching | Implemented, active | `CalcMatchY` `OKSeq` 17/18 preserves bounded flank reconstruction, signed mapping, smoothing, checkpoints, standard thresholds, and the `ConsensusOK` raw-tree topology filter; unavailable/fallback cases remain explicit |
| FinalTrim/ConsensusOK membership | Implemented through selected-tree cleanup, unvalidated | Both expansions, selected-role pruning, `OKSeq 15`, complete score/rebuild/fallback, shared raw/direct outlier cleanup, and strict inlier admission determine the distance set |
| Primary-RDP post-group recheck | Implemented, active, unvalidated | Finalized nonrepresentative rows are rerun against both role representatives with the native `LowP * 100000` lift; emitted/candidate/overlap/significant counts and best tract probabilities are retained |
| MaxChi strongest-peak recheck | Implemented, active, unvalidated | Supplied `FastRecCheckMC2` variable-site profiles, critical-difference gate, rolling χ² scan, `MissingData`/erasure and linear-edge bans, source `ChiPVal2P`, window growth, and three probability levels corroborate each event triplet and finalized nonrepresentative distance-list row |
| MaxChi exploratory discovery | Not yet ported | Primary RDP remains the only event-discovery/coordinate authority; smoothing, multi-peak destruction, event construction, and native retry scheduling remain explicit work |
| Event phylogenetic matrices | Implemented, unvalidated | Six Jukes–Cantor matrices with compact direct/reference caching |
| Bootstrap neighbour joining | Implemented, unvalidated | Deterministic ten-replicate event path; internal branches below 50% support collapse to zero length |
| Phylogenetic-correlation set | Implemented, unvalidated | Paired collapsed/raw tree affinity; labelled distance fallback outside the bounded tree panel |
| Complete two-of-three group | Implemented, unvalidated | Union of sequences present in any two of detectable, finalized `ConsensusOK` distance, and phylogenetic sets |
| Recombinant/parent recommendation | Implemented, partial | Nine auditable metrics; eight use mapped native full/half weights, while unported native method families keep `nativeWeightParity` false |
| Event review and correction | Implemented | Ordered accept/reject, role recommendation, manual roles/breakpoints, manual complete-group membership, and explicit downstream invalidation |
| Graphical breakpoint alignment | Implemented, unvalidated | On-demand original-alignment windows at both breakpoints; role/group/evidence row priority, circular wrapping, and parent-aware base colouring; bounded to 64 rows and 100 flanking sites per side |
| Native RDP `CheckEnds` breakpoint uncertainty | Implemented, active, unvalidated | Reconstructs the erased-triplet information map, source-shaped input `MissingData`, asymmetric beginning/ending DLL ranges, literal wrap comparisons, and linear edge gates; JSON/CSV/UI split native warning, erased-event IDs, immediate contact, input-missing, and edge reasons |
| BURT/BenHMM statistical breakpoint confidence | Implemented, active, unvalidated | Default-enabled `PolishBP(20)` setting → supplied three-state/three-symbol 21-start seeded Viterbi training → forward/reverse posterior → strict 99%/95% source ranges → signed matching, repositioning, missing-data repair, and revert guards; disabled runs preserve coordinates and results stay distinct from manual parent-state brackets |
| Graphical tree review | Implemented, unvalidated | Lazy side-by-side SVG review of all three paired regional comparisons; exact saved NJ edge lengths, bootstrap labels/collapse state, current roles/groups, and retained-fragment provenance are exposed without rebuilding trees |
| Downstream re-identification | Implemented, unvalidated | Corrected accepted events are re-erased; rejected tracts are restored; fixed earlier calls remain auditable while all later calls are rediscovered |
| Project/CSV export | Implemented | Reloadable project schema `v1alpha9`; MaxChi/RDP rechecks, final-trim/rebuild stages, round factors, fragment provenance, evidence sets, roles, traces, and decisions are retained |
| Project import/resume | Implemented, unvalidated | Schemas `v1alpha1`–`v1alpha9` restore primary state; pending repairs retain/remap only the valid signal/event prefix before downstream discovery resumes |
| Checkpoint-loss protection | Implemented, unvalidated | Completed/edited analyses visibly become dirty; tab exit and destructive dataset/settings replacement warn until project JSON is downloaded |
| Sequence-curation FASTA | Implemented, unvalidated | Full, enabled-only, and masked-or-disabled-only downloads work immediately after alignment loading, preserve complete aligned rows, use the current disjoint UI state, and report an empty complementary selection explicitly |
| Recombinant-sequence removal FASTA | Implemented, unvalidated | Omits every sequence in an accepted current co-recombinant group; core and UI reject incomplete review, stale reconciliation, and an empty result |
| Recombinant-column removal FASTA | Implemented, unvalidated | Deletes the inclusive union of every accepted event tract from all rows in linear time; core and UI reject incomplete review, stale reconciliation, and an empty result |
| Recombination-free FASTA | Implemented, unvalidated | Accepted tracts are replaced by gaps for current complete co-recombinant groups; core and UI reject incomplete-review export |
| Mosaic-fragment FASTA | Implemented, unvalidated | Event-ordered tract-masked originals plus aligned fragment-only records; core and UI reject incomplete-review export |
| GitHub Pages deployment | Configured, unvalidated | Default-branch/manual Actions workflow uses current official Node-24-generation actions, a Node 20 project build, locked npm dependencies, Emscripten 5.0.1, ABI/version/schema and TypeScript checks, single-worker WASM/Vite build, hidden-file-aware artifact verification, and official Pages deployment |
| Exact late native consensus | RDP complete plus MaxChi corroboration, unvalidated | `OKSeq` 0–18, expansions, selected-role pruning, complete CScore/rebuild/fallback, shared selected-role tree cleanup, finalized-list RDP signal rechecks, and bounded MaxChi strongest-peak rechecks are active; the full native MaxChi scheduler and other method rechecks remain |
| GENECONV and other methods | Not yet ported | Queued after primary RDP and MaxChi confirmation parity work |
| Compilation and runtime testing | Not performed locally | The workflow is configured to compile/check on GitHub only when the user runs it; this source checkpoint did not invoke it |
| Native-vs-WASM golden parity suite | Designed only | See `docs/validation-plan.md` |

## Recommended next checkpoint

1. Port MaxChi smoothing, ordered multi-peak destruction/retry, event construction, and exploratory
   scheduling while preserving RDP as an independently auditable method stream.
2. Add focused native golden fixtures for both MaxChi confirmation and active `PolishBP`/BURT,
   especially missing-data boundary windows, linear ends, `MaxX = 0`, grown-window correction,
   Microsoft-CRT random starts, signed wrapped ranges, and missing-data relocation.
3. Add the remaining recombinant-identification method families and their post-group rechecks
   without hiding the current per-method scores or native full/half contributions.
4. Only after explicit permission to execute code, establish small golden primary/secondary results
   against the supplied desktop application before calling the port parity validated.
