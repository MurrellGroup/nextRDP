import { AlertTriangle, Check, LoaderCircle, Pause, Play, ShieldCheck } from "lucide-react";

import type { ScanOptions, ScanProgress } from "../lib/types";

interface ScanStepProps {
  options: ScanOptions;
  sequenceCount: number;
  tripletCount: number;
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

export function ScanStep({
  options,
  sequenceCount,
  tripletCount,
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

      <div className="scan-console content-card">
        <div className="scan-summary">
          <div>
            <span className="eyebrow">Exploratory plan</span>
            <h2>{progress.scanRound > 1 ? `Cyclic RDP pass ${progress.scanRound}` : "Primary RDP screen"}</h2>
          </div>
          <span className={`run-state run-${progress.state}`}>
            {running ? <LoaderCircle className="spin" size={16} /> : hasResults ? <Check size={16} /> : <Pause size={16} />}
            {running ? "Running" : hasResults ? "Complete" : "Ready"}
          </span>
        </div>

        <div className="scan-facts">
          <div>
            <span>Active sequences</span>
            <strong>{integer.format(sequenceCount)}</strong>
          </div>
          <div>
            <span>Unique triplets</span>
            <strong>{integer.format(tripletCount)}</strong>
          </div>
          <div>
            <span>Topology</span>
            <strong>{options.circular ? "Circular" : "Linear"}</strong>
          </div>
          <div>
            <span>Correction</span>
            <strong>{options.correction === "bonferroni" ? "Bonferroni" : "None"}</strong>
          </div>
        </div>

        <div className="progress-block" aria-live="polite">
          <div className="progress-copy">
            <span>
              Round {progress.scanRound} · {integer.format(progress.processedTriplets)} / {integer.format(progress.totalTriplets || tripletCount)} triplets
            </span>
            <strong>{percentage}%</strong>
          </div>
          <div
            className="progress-track"
            role="progressbar"
            aria-valuemin={0}
            aria-valuemax={100}
            aria-valuenow={percentage}
          >
            <span style={{ width: `${Math.min(100, Math.max(0, percentage))}%` }} />
          </div>
          <div className="progress-meta">
            <span>
              {integer.format(progress.signalCount)} signals · {integer.format(progress.eventCount)} event candidates
            </span>
            <span>{integer.format(progress.cumulativeTriplets)} cumulative comparisons · window {options.windowSites}</span>
          </div>
        </div>

        <ol className="phase-list">
          <li className={hasResults || progress.phase !== "primary" ? "is-done" : running ? "is-active" : ""}>
            <span>{hasResults || progress.phase !== "primary" ? <Check size={16} /> : running ? <LoaderCircle className="spin" size={16} /> : "1"}</span>
            <div>
              <strong>Detect primary RDP signals</strong>
              <small>Information-rich triplets, rolling pair identities, binomial tail tests.</small>
            </div>
          </li>
          <li className={hasResults || progress.phase === "reconciliation" ? "is-done" : progress.phase === "cyclic-rescan" ? "is-active" : ""}>
            <span>{hasResults || progress.phase === "reconciliation" ? <Check size={16} /> : progress.phase === "cyclic-rescan" ? <LoaderCircle className="spin" size={16} /> : "2"}</span>
            <div>
              <strong>Erase and re-screen cyclically</strong>
              <small>Strongest event, three-set co-group, tract erasure, fragment re-entry, then a complete fresh triplet pass.</small>
            </div>
          </li>
          <li className={hasResults ? "is-done" : progress.phase === "reconciliation" ? "is-active" : ""}>
            <span>{hasResults ? <Check size={16} /> : progress.phase === "reconciliation" ? <LoaderCircle className="spin" size={16} /> : "3"}</span>
            <div>
              <strong>Finalize event evidence</strong>
              <small>Stop when no signal remains, preserve event order, and prepare the review hypothesis.</small>
            </div>
          </li>
        </ol>

        <div className="scan-controls">
          {!running && !hasResults ? (
            <button className="button button-primary button-large" type="button" onClick={onStart}>
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
              {progress.phase === "reconciliation" ? "Finalizing evidence sets" : "Stop after current batch"}
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
