import { AlertTriangle, FileText, Search, ShieldCheck, UploadCloud } from "lucide-react";
import { useMemo, useRef, useState } from "react";

import type { DatasetSummary } from "../lib/types";
import { Metric } from "./Metric";

interface DatasetStepProps {
  engineReady: boolean;
  engineMessage?: string;
  dataset: DatasetSummary | null;
  filename: string;
  fileSize: number;
  masked: Set<number>;
  busy: boolean;
  onLoad: (file: File) => void;
  onMaskChange: (index: number, masked: boolean) => void;
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
  busy,
  onLoad,
  onMaskChange,
  onContinue,
}: DatasetStepProps) {
  const input = useRef<HTMLInputElement>(null);
  const [dragging, setDragging] = useState(false);
  const [filter, setFilter] = useState("");
  const activeSequenceCount = Math.max(0, (dataset?.sequenceCount ?? 0) - masked.size);
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
            <Metric label="Triplets" value={integer.format(activeTripletCount)} detail={`${masked.size} masked`} />
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
                  Masked sequences are excluded from the initial triplet search, but remain in the
                  project for later trace-signal checks.
                </p>
              </div>
              <label className="search-box">
                <Search size={16} />
                <span className="sr-only">Filter sequence names</span>
                <input
                  value={filter}
                  onChange={(event) => setFilter(event.target.value)}
                  placeholder="Filter sequences"
                />
              </label>
            </div>
            <div className="sequence-table" role="table" aria-label="Alignment sequences">
              <div className="sequence-row sequence-header" role="row">
                <span role="columnheader">Use</span>
                <span role="columnheader">Sequence</span>
                <span role="columnheader">Valid sites</span>
                <span role="columnheader">Missing</span>
              </div>
              {visibleSequences.map((sequence) => {
                const isMasked = masked.has(sequence.index);
                return (
                  <label className={`sequence-row${isMasked ? " is-masked" : ""}`} key={sequence.index}>
                    <span>
                      <input
                        type="checkbox"
                        checked={!isMasked}
                        onChange={(event) => onMaskChange(sequence.index, !event.target.checked)}
                      />
                    </span>
                    <strong title={sequence.name}>{sequence.name}</strong>
                    <span>{integer.format(sequence.validSites)}</span>
                    <span>{percent.format(sequence.missingFraction)}</span>
                  </label>
                );
              })}
            </div>
            {dataset.sequences.length > 500 ? (
              <p className="table-footnote">Showing the first 500 matching sequences for browser rendering performance.</p>
            ) : null}
            <div className="manual-note">
              <strong>RDP5 dataset check</strong>
              <span>
                The manual’s approximate undetectable-distance threshold is currently{" "}
                {dataset.recommendedMinimumDistance.toFixed(4)} for this dataset.
              </span>
            </div>
          </div>

          <footer className="step-actions">
            <span>{activeSequenceCount} sequences will enter the primary scan.</span>
            <button className="button button-primary" type="button" onClick={onContinue}>
              Set analysis options
            </button>
          </footer>
        </>
      ) : null}
    </section>
  );
}
