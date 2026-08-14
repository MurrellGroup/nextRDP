import { AlertTriangle, Check, Clock3, Cpu, LoaderCircle, Pause, Play, ShieldCheck } from "lucide-react";

import type { ScanOptions, ScanProgress } from "../lib/types";

interface ScanStepProps {
  options: ScanOptions;
  sequenceCount: number;
  tripletCount: number;
  correctionTestCount: number;
  querySequenceCount: number;
  referenceSequenceCount: number;
  referenceGroupCount: number;
  progress: ScanProgress;
  running: boolean;
  error: string;
  hasResults: boolean;
  onStart: () => void;
  onCancel: () => void;
  onBack: () => void;
  onReview: () => void;
}

const integer = new Intl.NumberFormat();

function formatDuration(milliseconds: number): string {
  const safe = Math.max(0, Number.isFinite(milliseconds) ? milliseconds : 0);
  if (safe < 1000) return `${Math.round(safe)} ms`;
  const totalSeconds = safe / 1000;
  if (totalSeconds < 60) return `${totalSeconds.toFixed(3)} s`;
  const minutes = Math.floor(totalSeconds / 60);
  const seconds = totalSeconds - minutes * 60;
  if (minutes < 60) return `${minutes}m ${seconds.toFixed(3).padStart(6, "0")}s`;
  const hours = Math.floor(minutes / 60);
  return `${hours}h ${String(minutes % 60).padStart(2, "0")}m ${seconds.toFixed(3).padStart(6, "0")}s`;
}

export function ScanStep({
  options,
  sequenceCount,
  tripletCount,
  correctionTestCount,
  querySequenceCount,
  referenceSequenceCount,
  referenceGroupCount,
  progress,
  running,
  error,
  hasResults,
  onStart,
  onCancel,
  onBack,
  onReview,
}: ScanStepProps) {
  const percentage = Math.round(progress.fraction * 1000) / 10;
  const indeterminate = running && progress.processedTriplets === 0;
  const timing = progress.timing;
  const recentRounds = timing?.rounds.slice(-3) ?? [];
  const currentCorrectionTests = progress.state === "idle"
    ? correctionTestCount
    : progress.correctionTests;
  const currentTotalTriplets = progress.state === "idle"
    ? tripletCount
    : progress.totalTriplets;
  const currentSequenceCount = progress.state === "idle"
    ? sequenceCount
    : progress.activeWorkingSequenceCount;
  const currentQueryCount = progress.state === "idle"
    ? querySequenceCount
    : progress.queryWorkingSequenceCount;
  const currentReferenceCount = progress.state === "idle"
    ? referenceSequenceCount
    : progress.referenceWorkingSequenceCount;
  const currentReferenceGroupCount = progress.state === "idle"
    ? referenceGroupCount
    : progress.activeReferenceGroupCount;
  const discoveryMethods = [
    "RDP",
    ...(options.geneconvEnabled ? ["GENECONV"] : []),
    ...(options.bootscanPrimaryEnabled ? ["BootScan"] : []),
    ...(options.maxChiEnabled ? ["MaxChi"] : []),
    ...(options.chimaeraEnabled ? ["CHIMAERA"] : []),
    ...(options.siscanPrimaryEnabled ? ["SISCAN"] : []),
    ...(options.threeSeqEnabled ? ["3SEQ"] : []),
  ];
  return (
    <section className="step-page scan-page" aria-labelledby="scan-title">
      <header className="page-heading">
        <div>
          <span className="eyebrow">03 · Scan</span>
          <h1 id="scan-title">Detect recombination signals</h1>
          <p>
            The engine runs away from the interface thread. You can stop between bounded triplet
            batches without freezing the dataset view.
          </p>
        </div>
        <div className="privacy-note">
          <ShieldCheck size={18} />
          <span>
            <strong>Worker-isolated WASM</strong>
            No upload or server round-trip.
          </span>
        </div>
      </header>

      {error ? (
        <div className="notice notice-red" role="alert">
          <AlertTriangle size={19} />
          <div>
            <strong>Scan stopped</strong>
            <p>{error}</p>
          </div>
        </div>
      ) : null}

      {hasResults && progress.cycleTermination === "user-stopped" ? (
        <div className="notice notice-blue" role="status">
          <Check size={19} />
          <div>
            <strong>Scan stopped cleanly</strong>
            <p>
              The unfinished cyclic round was discarded. Events completed in earlier rounds were
              retained and finalized for review.
            </p>
          </div>
        </div>
      ) : null}

      <div className="scan-console content-card">
        <div className="scan-summary">
          <div>
            <span className="eyebrow">
              {options.analysisMode === "query-reference" ? "Query vs reference plan" : "Exploratory plan"}
            </span>
            <h2>{progress.scanRound > 1 ? `Cyclic discovery pass ${progress.scanRound}` : "Primary triplet screen"}</h2>
          </div>
          <span className={`run-state run-${progress.state}`}>
            {running ? <LoaderCircle className="spin" size={16} /> : hasResults ? <Check size={16} /> : <Pause size={16} />}
            {running ? "Running" : hasResults ? "Complete" : "Ready"}
          </span>
        </div>

        <div className="scan-facts">
          <div>
            <span>Active working rows</span>
            <strong>{integer.format(currentSequenceCount)}</strong>
          </div>
          <div>
            <span>{options.analysisMode === "query-reference" ? "Constrained triplets" : "Unique triplets"}</span>
            <strong>{integer.format(tripletCount)}</strong>
          </div>
          <div>
            <span>Analysis scheme</span>
            <strong>{options.analysisMode === "query-reference" ? "Query vs reference" : "Fully exploratory"}</strong>
          </div>
          {options.analysisMode === "query-reference" ? (
            <div>
              <span>Input roles</span>
              <strong>
                {integer.format(currentQueryCount)} Q · {integer.format(currentReferenceCount)} R · {integer.format(currentReferenceGroupCount)} groups
              </strong>
            </div>
          ) : null}
          <div>
            <span>Topology</span>
            <strong>{options.circular ? "Circular" : "Linear"}</strong>
          </div>
          <div>
            <span>Discovery methods</span>
            <strong>{discoveryMethods.join(" + ")}</strong>
          </div>
          <div>
            <span>Correction</span>
            <strong>
              {options.correction === "bonferroni"
                ? options.threeSeqEnabled
                  ? `${integer.format(currentCorrectionTests)} opportunities · × for RDP family · Dunn–Šidák for 3SEQ`
                  : `Bonferroni × ${integer.format(currentCorrectionTests)}`
                : "None"}
            </strong>
          </div>
          <div>
            <span>Breakpoint confidence</span>
            <strong>{options.polishBreakpoints ? "BURT enabled" : "Preserve detected calls"}</strong>
          </div>
          <div>
            <span>Compute</span>
            <strong>
              {progress.state === "idle"
                ? `${options.cpuThreads} CPU${options.cpuThreads === 1 ? "" : "s"} requested`
                : `${progress.execution.activeThreads} CPU${progress.execution.activeThreads === 1 ? "" : "s"} · ${progress.execution.mode === "wasm-pthreads" ? "WASM threads" : "single worker"}`}
            </strong>
          </div>
        </div>

        <div className="progress-block" aria-live="polite">
          <div className="progress-copy">
            <span>
              Round {progress.scanRound} · {integer.format(progress.processedTriplets)} / {integer.format(currentTotalTriplets)} triplets
            </span>
            <strong>{indeterminate ? "Working…" : `${percentage}%`}</strong>
          </div>
          <div
            className={`progress-track${running ? " is-active" : ""}${indeterminate ? " is-indeterminate" : ""}`}
            role="progressbar"
            aria-valuemin={0}
            aria-valuemax={100}
            aria-valuenow={indeterminate ? undefined : percentage}
            aria-valuetext={indeterminate ? `Round ${progress.scanRound} is starting` : `${percentage}%`}
          >
            <span style={indeterminate ? undefined : { width: `${Math.min(100, Math.max(0, percentage))}%` }} />
          </div>
          {timing ? (
            <div className="timing-panel" aria-label="Analysis timing">
              <div className="timing-total">
                <Clock3 size={15} />
                <span>{running ? "Elapsed" : "Total run time"}</span>
                <strong>{formatDuration(timing.totalMs)}</strong>
              </div>
              <dl className="timing-phases">
                <div><dt>Setup</dt><dd>{formatDuration(timing.setupMs)}</dd></div>
                <div><dt>Primary scan</dt><dd>{formatDuration(timing.primaryMs)}</dd></div>
                <div><dt>Cyclic rescans</dt><dd>{formatDuration(timing.cyclicRescanMs)}</dd></div>
                <div><dt>Final evidence</dt><dd>{formatDuration(timing.reconciliationMs)}</dd></div>
              </dl>
              {running && timing.currentRoundMs > 0 ? (
                <div className="timing-rounds">
                  <span>Current round {progress.scanRound}: {formatDuration(timing.currentRoundMs)}</span>
                  {recentRounds.map((round) => (
                    <span key={round.round}>R{round.round}: {formatDuration(round.elapsedMs)}</span>
                  ))}
                </div>
              ) : timing.rounds.length > 0 ? (
                <div className="timing-rounds">
                  <span>{timing.rounds.length} timed round{timing.rounds.length === 1 ? "" : "s"}</span>
                  {recentRounds.map((round) => (
                    <span key={round.round}>R{round.round}: {formatDuration(round.elapsedMs)}</span>
                  ))}
                </div>
              ) : null}
            </div>
          ) : null}
          <div className="progress-meta">
            <span>
              {integer.format(progress.signalCount)} signals · {integer.format(progress.eventCount)} event candidates
            </span>
            <span>
              {integer.format(progress.cumulativeTriplets)} cumulative triplets · RDP {options.windowSites}
              {options.maxChiEnabled ? ` · MaxChi ${options.maxChiWindowSites}` : ""}
              {options.chimaeraEnabled ? ` · CHIMAERA ${options.chimaeraWindowSites}` : ""}
              {options.geneconvEnabled
                ? ` · GENECONV G${options.geneconvMismatchScale}/overlap ${options.geneconvMaxOverlaps}`
                : ""}
              {options.bootscanPrimaryEnabled
                ? ` · BootScan ${options.bootscanWindowSites}/${options.bootscanStepSites} × ${options.bootscanBootstrapReplicates}`
                : ""}
              {options.siscanPrimaryEnabled
                ? ` · SISCAN ${options.siscanWindowSites}/${options.siscanStepSites} × ${options.siscanScanPermutations}/${options.siscanPValuePermutations}`
                : ""}
              {options.threeSeqEnabled ? " · 3SEQ exact/Siegmund" : ""}
            </span>
          </div>
          {progress.tripletSummariesReused > 0 || progress.fragmentSequencesPruned > 0 ? (
            <div className="progress-meta">
              <span>
                {integer.format(progress.cleanTripletsPruned)} clean triplets permanently pruned · {integer.format(progress.cachedSignalsReused)} cached signals replayed
              </span>
              <span>
                {integer.format(progress.methodScansSkipped)} method scans skipped · {integer.format(progress.fragmentSequencesPruned)} event-free fragments removed
              </span>
            </div>
          ) : null}
          {progress.invalidScheduleTripletsSkipped > 0 ? (
            <div className="progress-meta">
              <span>
                {integer.format(progress.tripletKernelEvaluations)} triplet kernel evaluations
              </span>
              <span>
                {integer.format(progress.invalidScheduleTripletsSkipped)} same-origin combinations skipped inside WASM batches
              </span>
            </div>
          ) : null}
          {options.maxChiEnabled ? (
            <div className="progress-meta">
              <span>
                {integer.format(progress.maxChiProfilesScanned)} MaxChi profiles · {integer.format(progress.maxChiPeakAttempts)} raw peaks examined
              </span>
              <span>{integer.format(progress.maxChiCandidatesFound)} corrected MaxChi candidates</span>
            </div>
          ) : null}
          {options.chimaeraEnabled ? (
            <div className="progress-meta">
              <span>
                {integer.format(progress.chimaeraProfilesScanned)} CHIMAERA target profiles · {integer.format(progress.chimaeraPeakAttempts)} raw peaks examined
              </span>
              <span>{integer.format(progress.chimaeraCandidatesFound)} corrected CHIMAERA candidates</span>
            </div>
          ) : null}
          {options.geneconvEnabled ? (
            <div className="progress-meta">
              <span>
                {integer.format(progress.geneconvFragmentsScored)} GENECONV positive starts · {integer.format(progress.geneconvQualifiedFragments)} above KA critical score
              </span>
              <span>
                {integer.format(progress.geneconvCandidatesFound)} corrected candidates · {integer.format(progress.geneconvOverlapRejections)} overlap rejections
              </span>
            </div>
          ) : null}
          {options.threeSeqEnabled ? (
            <div className="progress-meta">
              <span>
                {integer.format(progress.threeSeqProfilesScanned)} 3SEQ target walks · {integer.format(progress.threeSeqExactEvaluations)} exact tails
              </span>
              <span>
                {integer.format(progress.threeSeqApproximateEvaluations)} Siegmund fallbacks · {integer.format(progress.threeSeqCandidatesFound)} threshold-passing candidates
              </span>
            </div>
          ) : null}
          {options.bootscanPrimaryEnabled ? (
            <div className="progress-meta">
              <span>
                {integer.format(progress.bootscanProfilesScanned)} BootScan triplet profiles · {integer.format(progress.bootscanCandidateRegionsScored)} support regions scored
              </span>
              <span>
                {integer.format(progress.bootscanCandidatesFound)} corrected candidates · {integer.format(progress.bootscanPairProfileCacheHits)} / {integer.format(progress.bootscanPairProfilesRequested)} pair-profile cache hits
              </span>
            </div>
          ) : null}
          {options.siscanPrimaryEnabled ? (
            <>
              <div className="progress-meta">
                <span>
                  {integer.format(progress.siscanProfilesScanned)} SISCAN triplet profiles · {integer.format(progress.siscanWindowsScored)} windows · {integer.format(progress.siscanCandidateRegionsScored)} pair-switch regions
                </span>
                <span>
                  {integer.format(progress.siscanCandidatesFound)} corrected candidates · {integer.format(progress.siscanPermutationDraws)} vertical-permutation draws
                </span>
              </div>
              <div className="progress-meta">
                <span>
                  {integer.format(progress.siscanContextBuilds)} WPGMA context build{progress.siscanContextBuilds === 1 ? "" : "s"} · {integer.format(progress.siscanContextPairComparisons)} direct pair comparisons · {integer.format(progress.siscanContextTreeMerges)} tree merges
                </span>
                <span>{integer.format(progress.siscanRandomValuesGenerated)} cached MakeVRand values</span>
              </div>
            </>
          ) : null}
        </div>

        {progress.maxChiPeakLimitTriplets > 0 ? (
          <div className="notice notice-amber" role="status">
            <AlertTriangle size={18} />
            <p>
              {integer.format(progress.maxChiPeakLimitTriplets)} triplet
              {progress.maxChiPeakLimitTriplets === 1 ? " reached" : "s reached"} the supplied
              100-peak retry bound while positive raw peaks remained. Any additional peaks for those
              triplets were deliberately not explored.
            </p>
          </div>
        ) : null}

        {progress.chimaeraPeakLimitTargets > 0 ? (
          <div className="notice notice-amber" role="status">
            <AlertTriangle size={18} />
            <p>
              {integer.format(progress.chimaeraPeakLimitTargets)} CHIMAERA target profile
              {progress.chimaeraPeakLimitTargets === 1 ? " reached" : "s reached"} the supplied
              100-peak retry bound while positive raw peaks remained. Later peaks in those target
              profiles were deliberately not explored.
            </p>
          </div>
        ) : null}

        {progress.geneconvNumericalFallbackTracks > 0 ? (
          <div className="notice notice-amber" role="status">
            <AlertTriangle size={18} />
            <p>
              {integer.format(progress.geneconvNumericalFallbackTracks)} GENECONV track
              {progress.geneconvNumericalFallbackTracks === 1 ? " used" : "s used"} the bounded
              root-bracketing fallback after the supplied Newton start became unstable. The
              fallback preserves a responsive worker and is retained in exported diagnostics.
            </p>
          </div>
        ) : null}

        <ol className="phase-list">
          <li className={hasResults || progress.phase !== "primary" ? "is-done" : running ? "is-active" : ""}>
            <span>{hasResults || progress.phase !== "primary" ? <Check size={16} /> : running ? <LoaderCircle className="spin" size={16} /> : "1"}</span>
            <div>
              <strong>Detect {discoveryMethods.join(", ")} signals</strong>
              <small>
                RDP information-rich tracts
                {options.maxChiEnabled
                  ? ", MaxChi pairwise raw χ² peaks with source-shaped tract inference"
                  : ""}
                {options.chimaeraEnabled
                  ? ", and CHIMAERA target-rotated parent-match profiles"
                  : ""}
                {options.geneconvEnabled
                  ? ", and GENECONV six-track KA fragment scores"
                  : ""}
                {options.siscanPrimaryEnabled
                  ? ", and SISCAN nearest-outlier vertical-permutation Z scores"
                  : ""}
                {options.threeSeqEnabled
                  ? ", and 3SEQ target-rotated exact hypergeometric random walks"
                  : ""}
                {options.correction === "bonferroni"
                  ? " with method-specific project correction and bounded peak retries."
                  : " with uncorrected significance and bounded peak retries."}
              </small>
            </div>
          </li>
          <li className={hasResults || progress.phase === "reconciliation" ? "is-done" : progress.phase === "cyclic-rescan" ? "is-active" : ""}>
            <span>{hasResults || progress.phase === "reconciliation" ? <Check size={16} /> : progress.phase === "cyclic-rescan" ? <LoaderCircle className="spin" size={16} /> : "2"}</span>
            <div>
              <strong>Erase and re-screen cyclically</strong>
              <small>
                Strongest event, three-set co-group, tract erasure, fragment re-entry, then XOverList-style summary reuse for unchanged triplets and fresh kernels for affected rows.
              </small>
            </div>
          </li>
          <li className={hasResults ? "is-done" : progress.phase === "reconciliation" ? "is-active" : ""}>
            <span>{hasResults ? <Check size={16} /> : progress.phase === "reconciliation" ? <LoaderCircle className="spin" size={16} /> : "3"}</span>
            <div>
              <strong>Finalize event evidence</strong>
              <small>
                {options.polishBreakpoints
                  ? "Polish detected breakpoints with BURT, preserve event order, and prepare the review hypothesis."
                  : "Preserve detected breakpoints, retain event order, and prepare the review hypothesis."}
              </small>
            </div>
          </li>
        </ol>

        <div className="scan-controls">
          <span className="compute-status">
            <Cpu size={15} /> {options.cpuThreads} CPU{options.cpuThreads === 1 ? "" : "s"} selected
          </span>
          {!running && !hasResults ? (
            <button
              className="button button-primary button-large"
              type="button"
              onClick={onStart}
              disabled={sequenceCount < 3 || tripletCount === 0}
            >
              <Play size={18} fill="currentColor" /> Start scan
            </button>
          ) : null}
          {running ? (
            <button
              className={progress.phase === "primary" || progress.phase === "cyclic-rescan" ? "button button-danger button-large" : "button button-disabled button-large"}
              type="button"
              onClick={onCancel}
              disabled={progress.phase === "reconciliation"}
            >
              {progress.phase === "reconciliation" ? <LoaderCircle className="spin" size={18} /> : <Pause size={18} />}
              {progress.phase === "reconciliation" ? "Finalizing evidence sets" : "Stop and review completed events"}
            </button>
          ) : null}
          {hasResults ? (
            <button className="button button-primary button-large" type="button" onClick={onReview}>
              Review events
            </button>
          ) : null}
        </div>
      </div>

      <footer className="step-actions">
        <button className="button button-quiet" type="button" onClick={onBack} disabled={running}>
          Back to settings
        </button>
        <span>Closing this tab stops the in-memory analysis.</span>
      </footer>
    </section>
  );
}
