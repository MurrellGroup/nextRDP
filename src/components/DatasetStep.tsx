import { AlertTriangle, FileText, Search, ShieldCheck, UploadCloud } from "lucide-react";
import { useMemo, useRef, useState } from "react";

import type { DatasetSummary, SequenceAnalysisState } from "../lib/types";
import { Metric } from "./Metric";

interface DatasetStepProps {
  engineReady: boolean;
  engineMessage?: string;
  dataset: DatasetSummary | null;
  filename: string;
  fileSize: number;
  masked: Set<number>;
  disabled: Set<number>;
  busy: boolean;
  onLoad: (file: File) => void;
  onSequenceStateChange: (index: number, state: SequenceAnalysisState) => void;
  onAllSequenceStatesChange: (
    action: "auto-mask" | "enable-all" | "mask-all" | "disable-all",
  ) => void;
  onExportFullAlignment: () => void;
  onExportEnabledSequences: () => void;
  onExportMaskedOrDisabledSequences: () => void;
  onContinue: () => void;
}

const integer = new Intl.NumberFormat();
const percent = new Intl.NumberFormat(undefined, { style: "percent", maximumFractionDigits: 1 });

function chooseThree(count: number): number {
  return count < 3 ? 0 : (count * (count - 1) * (count - 2)) / 6;
}

function bytes(value: number): string {
  if (value < 1024) return `${value} B`;
  if (value < 1024 * 1024) return `${(value / 1024).toFixed(1)} KB`;
  return `${(value / (1024 * 1024)).toFixed(1)} MB`;
}

export function DatasetStep({
  engineReady,
  engineMessage,
  dataset,
  filename,
  fileSize,
  masked,
  disabled,
  busy,
  onLoad,
  onSequenceStateChange,
  onAllSequenceStatesChange,
  onExportFullAlignment,
  onExportEnabledSequences,
  onExportMaskedOrDisabledSequences,
  onContinue,
}: DatasetStepProps) {
  const input = useRef<HTMLInputElement>(null);
  const [dragging, setDragging] = useState(false);
  const [filter, setFilter] = useState("");
  const activeSequenceCount = Math.max(
    0,
    (dataset?.sequenceCount ?? 0) - masked.size - disabled.size,
  );
  const activeTripletCount = chooseThree(activeSequenceCount);

  const visibleSequences = useMemo(() => {
    if (!dataset) return [];
    const needle = filter.trim().toLowerCase();
    const matches = needle
      ? dataset.sequences.filter((sequence) => sequence.name.toLowerCase().includes(needle))
      : dataset.sequences;
    return matches.slice(0, 500);
  }, [dataset, filter]);

  const chooseFile = (files: FileList | null) => {
    const file = files?.item(0);
    if (file) onLoad(file);
  };

  return (
    <section className="step-page" aria-labelledby="dataset-title">
      <header className="page-heading">
        <div>
          <span className="eyebrow">01 · Dataset</span>
          <h1 id="dataset-title">Start with the alignment</h1>
          <p>
            RDP analyses aligned nucleotides. Inspect diversity and missing data before choosing
            a scan—not after a surprising result appears.
          </p>
        </div>
        <div className="privacy-note">
          <ShieldCheck size={18} />
          <span>
            <strong>Local by design</strong>
            Your alignment stays in this browser.
          </span>
        </div>
      </header>

      {!engineReady && engineMessage ? (
        <div className="notice notice-amber" role="status">
          <AlertTriangle size={19} />
          <div>
            <strong>WASM engine not loaded</strong>
            <p>{engineMessage}</p>
          </div>
        </div>
      ) : null}

      <input
        ref={input}
        type="file"
        accept=".fas,.fasta,.fa,.aln,.phy,.phylip,.nex,.nexus,.meg,.mega,.gde,.txt,.json,application/json"
        onChange={(event) => {
          chooseFile(event.target.files);
          event.currentTarget.value = "";
        }}
        hidden
      />
      <button
        type="button"
        className={`drop-zone${dragging ? " is-dragging" : ""}${dataset ? " is-compact" : ""}`}
        onClick={() => input.current?.click()}
        onDragEnter={(event) => {
          event.preventDefault();
          setDragging(true);
        }}
        onDragOver={(event) => event.preventDefault()}
        onDragLeave={() => setDragging(false)}
        onDrop={(event) => {
          event.preventDefault();
          setDragging(false);
          chooseFile(event.dataTransfer.files);
        }}
        disabled={!engineReady || busy}
      >
        <span className="drop-icon">
          {dataset ? <FileText size={25} /> : <UploadCloud size={27} />}
        </span>
        <span>
          <strong>{busy ? "Reading alignment…" : dataset ? "Replace alignment" : "Drop an alignment here"}</strong>
          <small>
            {dataset
              ? `${filename} · ${bytes(fileSize)}`
              : "or resume .rdpweb.json · FASTA, GDE, CLUSTAL, PHYLIP, NEXUS, MEGA"}
          </small>
        </span>
      </button>

      {dataset ? (
        <>
          <div className="metrics-grid metrics-five">
            <Metric label="Sequences" value={integer.format(dataset.sequenceCount)} />
            <Metric label="Alignment" value={integer.format(dataset.alignmentLength)} detail="nucleotide columns" />
            <Metric label="Variable sites" value={integer.format(dataset.variableSiteCount)} />
            <Metric label="Mean identity" value={dataset.meanPairIdentity == null ? "—" : percent.format(dataset.meanPairIdentity)} />
            <Metric
              label="Triplets"
              value={integer.format(activeTripletCount)}
              detail={`${masked.size} masked · ${disabled.size} disabled`}
            />
          </div>

          {dataset.warnings.map((warning) => (
            <div className="notice notice-amber" key={warning}>
              <AlertTriangle size={18} />
              <p>{warning}</p>
            </div>
          ))}

          <div className="content-card sequence-card">
            <div className="card-heading sequence-heading">
              <div>
                <span className="eyebrow">Sequence curation</span>
                <h2>Choose the exploratory set</h2>
                <p>
                  Enabled rows enter every screen. Masked rows skip the primary triplet catalogue
                  but remain in secondary checks and trees; disabled rows remain only as tree context.
                </p>
              </div>
              <div className="sequence-heading-controls">
                <label className="search-box">
                  <Search size={16} />
                  <span className="sr-only">Filter sequence names</span>
                  <input
                    value={filter}
                    onChange={(event) => setFilter(event.target.value)}
                    placeholder="Filter sequences"
                  />
                </label>
                <div className="sequence-bulk-actions" aria-label="Bulk sequence curation">
                  <button
                    className="button button-secondary button-compact"
                    type="button"
                    onClick={() => onAllSequenceStatesChange("auto-mask")}
                    title="Restore the supplied RDP closest-pair auto-mask recommendation"
                    disabled={busy}
                  >
                    Auto-mask
                  </button>
                  <button
                    className="button button-quiet button-compact"
                    type="button"
                    onClick={() => onAllSequenceStatesChange("enable-all")}
                    disabled={busy}
                  >
                    Enable all
                  </button>
                  <button
                    className="button button-quiet button-compact"
                    type="button"
                    onClick={() => onAllSequenceStatesChange("mask-all")}
                    disabled={busy}
                  >
                    Mask all
                  </button>
                  <button
                    className="button button-quiet button-compact"
                    type="button"
                    onClick={() => onAllSequenceStatesChange("disable-all")}
                    disabled={busy}
                  >
                    Disable all
                  </button>
                </div>
              </div>
            </div>
            <div className="sequence-table" role="table" aria-label="Alignment sequences">
              <div className="sequence-row sequence-header" role="row">
                <span role="columnheader">Analysis</span>
                <span role="columnheader">Sequence</span>
                <span role="columnheader">Valid sites</span>
                <span role="columnheader">Missing</span>
              </div>
              {visibleSequences.map((sequence) => {
                const isMasked = masked.has(sequence.index);
                const isDisabled = disabled.has(sequence.index);
                const analysisState: SequenceAnalysisState = isDisabled
                  ? "disabled"
                  : isMasked
                    ? "masked"
                    : "enabled";
                return (
                  <div
                    className={`sequence-row${isMasked ? " is-masked" : ""}${isDisabled ? " is-disabled" : ""}`}
                    key={sequence.index}
                    role="row"
                  >
                    <span>
                      <select
                        className={`sequence-state sequence-state-${analysisState}`}
                        aria-label={`Analysis state for ${sequence.name}`}
                        value={analysisState}
                        disabled={busy}
                        onChange={(event) => onSequenceStateChange(
                          sequence.index,
                          event.target.value as SequenceAnalysisState,
                        )}
                      >
                        <option value="enabled">Enabled</option>
                        <option value="masked">Masked</option>
                        <option value="disabled">Disabled</option>
                      </select>
                    </span>
                    <strong title={sequence.name}>{sequence.name}</strong>
                    <span>{integer.format(sequence.validSites)}</span>
                    <span>{percent.format(sequence.missingFraction)}</span>
                  </div>
                );
              })}
            </div>
            {dataset.sequences.length > 500 ? (
              <p className="table-footnote">Showing the first 500 matching sequences for browser rendering performance.</p>
            ) : null}
            <div className="sequence-export-actions">
              <span>
                Save the full alignment or either curated partition now; no scan is required.
              </span>
              <div>
                <button
                  className="button button-quiet button-compact"
                  type="button"
                  onClick={onExportFullAlignment}
                  disabled={busy}
                >
                  Save full FASTA
                </button>
                <button
                  className="button button-quiet button-compact"
                  type="button"
                  onClick={onExportEnabledSequences}
                  disabled={busy || activeSequenceCount === 0}
                >
                  Save enabled FASTA
                </button>
                <button
                  className="button button-quiet button-compact"
                  type="button"
                  onClick={onExportMaskedOrDisabledSequences}
                  disabled={busy || masked.size + disabled.size === 0}
                >
                  Save masked / disabled FASTA
                </button>
              </div>
            </div>
            <div className="manual-note">
              <strong>RDP5 dataset check</strong>
              <span>
                The manual’s approximate undetectable-distance threshold is currently{" "}
                {dataset.recommendedMinimumDistance.toFixed(4)} for this dataset.
              </span>
            </div>
          </div>

          <footer className="step-actions">
            <span>
              {activeSequenceCount < 3
                ? "Enable at least three unmasked sequences to continue."
                : `${activeSequenceCount} sequences will enter the primary scan.`}
            </span>
            <button
              className="button button-primary"
              type="button"
              onClick={onContinue}
              disabled={busy || activeSequenceCount < 3}
            >
              Set analysis options
            </button>
          </footer>
        </>
      ) : null}
    </section>
  );
}
