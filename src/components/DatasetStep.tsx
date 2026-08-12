import { AlertTriangle, FileText, Info, Search, ShieldCheck, UploadCloud } from "lucide-react";
import { useEffect, useMemo, useRef, useState } from "react";

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
  referenceGroups: Map<number, number>;
  eligibleSequenceCount: number;
  exploratoryTripletCount: number;
  busy: boolean;
  onLoad: (file: File) => void;
  onSequenceStateChange: (index: number, state: SequenceAnalysisState) => void;
  onAllSequenceStatesChange: (
    action: "auto-mask" | "enable-all" | "mask-all" | "disable-all",
  ) => void;
  onReferenceGroupChange: (index: number, group: number) => void;
  onReferenceGroupsChange: (indices: number[], group: number) => void;
  onAllReferenceGroupsChange: (action: "detect" | "compact" | "clear") => void;
  onExportFullAlignment: () => void;
  onExportEnabledSequences: () => void;
  onExportMaskedOrDisabledSequences: () => void;
  onContinue: () => void;
}

const integer = new Intl.NumberFormat();
const percent = new Intl.NumberFormat(undefined, { style: "percent", maximumFractionDigits: 1 });

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
  referenceGroups,
  eligibleSequenceCount,
  exploratoryTripletCount,
  busy,
  onLoad,
  onSequenceStateChange,
  onAllSequenceStatesChange,
  onReferenceGroupChange,
  onReferenceGroupsChange,
  onAllReferenceGroupsChange,
  onExportFullAlignment,
  onExportEnabledSequences,
  onExportMaskedOrDisabledSequences,
  onContinue,
}: DatasetStepProps) {
  const input = useRef<HTMLInputElement>(null);
  const [dragging, setDragging] = useState(false);
  const [filter, setFilter] = useState("");
  const [selectedSequenceIndices, setSelectedSequenceIndices] = useState<Set<number>>(new Set());
  const [bulkReferenceGroup, setBulkReferenceGroup] = useState("");
  const curatedEnabledCount = Math.max(
    0,
    (dataset?.sequenceCount ?? 0) - masked.size - disabled.size,
  );
  const referenceGroupCount = new Set(referenceGroups.values()).size;

  const matchingSequences = useMemo(() => {
    if (!dataset) return [];
    const needle = filter.trim().toLowerCase();
    return needle
      ? dataset.sequences.filter((sequence) => sequence.name.toLowerCase().includes(needle))
      : dataset.sequences;
  }, [dataset, filter]);
  const visibleSequences = matchingSequences.slice(0, 500);
  const allMatchingSelected = matchingSequences.length > 0 && matchingSequences.every(
    (sequence) => selectedSequenceIndices.has(sequence.index),
  );
  const normalizedBulkGroup = Number(bulkReferenceGroup);
  const bulkGroupIsValid = Number.isInteger(normalizedBulkGroup) &&
    normalizedBulkGroup > 0 && normalizedBulkGroup <= 0xffff_ffff;

  useEffect(() => {
    setFilter("");
    setSelectedSequenceIndices(new Set());
    setBulkReferenceGroup("");
  }, [dataset]);

  const toggleSequenceSelection = (index: number) => {
    setSelectedSequenceIndices((current) => {
      const next = new Set(current);
      if (next.has(index)) next.delete(index);
      else next.add(index);
      return next;
    });
  };

  const toggleMatchingSelection = () => {
    setSelectedSequenceIndices((current) => {
      const next = new Set(current);
      matchingSequences.forEach((sequence) => {
        if (allMatchingSelected) next.delete(sequence.index);
        else next.add(sequence.index);
      });
      return next;
    });
  };

  const applySelectedReferenceGroup = (group: number) => {
    if (selectedSequenceIndices.size === 0) return;
    onReferenceGroupsChange([...selectedSequenceIndices], group);
  };

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
              label="Exploratory triplets"
              value={integer.format(exploratoryTripletCount)}
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
                <h2>Choose the primary analysis set</h2>
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
            <div className="reference-toolbar">
              <div>
                <strong>Query vs reference roles</strong>
                <span>
                  Leave a group blank for a query. References with the same positive group number
                  are never paired together; the scan scheme is selected on the next page. Detecting
                  names does not silently change enabled, masked, or disabled states.
                </span>
              </div>
              <div>
                <button
                  className="button button-quiet button-compact"
                  type="button"
                  onClick={() => onAllReferenceGroupsChange("detect")}
                  disabled={busy}
                  title="Infer documented REF-A&lt;name&gt; style reference groups"
                >
                  Detect REF names
                </button>
                <button
                  className="button button-quiet button-compact"
                  type="button"
                  onClick={() => onAllReferenceGroupsChange("compact")}
                  disabled={busy || referenceGroups.size === 0}
                  title="Renumber groups by first appearance without changing membership"
                >
                  Compact groups
                </button>
                <button
                  className="button button-quiet button-compact"
                  type="button"
                  onClick={() => onAllReferenceGroupsChange("clear")}
                  disabled={busy || referenceGroups.size === 0}
                >
                  All queries
                </button>
              </div>
            </div>
            {referenceGroups.size > 0 ? (
              <div className="notice notice-blue reference-detection-note">
                <Info size={17} />
                <p>
                  {referenceGroups.size.toLocaleString()} reference-labelled sequence{referenceGroups.size === 1 ? "" : "s"} across {referenceGroupCount.toLocaleString()} group{referenceGroupCount === 1 ? "" : "s"}. These assignments affect primary triplets only if you choose query vs reference on the settings page.
                </p>
              </div>
            ) : null}
            <div className="reference-selection-toolbar" aria-label="Bulk query and reference assignment">
              <span>
                {selectedSequenceIndices.size > 0
                  ? `${selectedSequenceIndices.size.toLocaleString()} selected`
                  : "Select rows, or select every row matching the current name filter."}
              </span>
              <div>
                <label className="bulk-reference-group-field">
                  <span>Reference group</span>
                  <input
                    type="number"
                    min="1"
                    max="4294967295"
                    step="1"
                    inputMode="numeric"
                    value={bulkReferenceGroup}
                    placeholder="e.g. 2"
                    disabled={busy}
                    onChange={(event) => setBulkReferenceGroup(event.target.value)}
                  />
                </label>
                <button
                  className="button button-secondary button-compact"
                  type="button"
                  onClick={() => applySelectedReferenceGroup(normalizedBulkGroup)}
                  disabled={busy || selectedSequenceIndices.size === 0 || !bulkGroupIsValid}
                >
                  Assign group
                </button>
                <button
                  className="button button-quiet button-compact"
                  type="button"
                  onClick={() => applySelectedReferenceGroup(0)}
                  disabled={busy || selectedSequenceIndices.size === 0}
                >
                  Make queries
                </button>
                <button
                  className="button button-quiet button-compact"
                  type="button"
                  onClick={() => setSelectedSequenceIndices(new Set())}
                  disabled={busy || selectedSequenceIndices.size === 0}
                >
                  Clear selection
                </button>
              </div>
            </div>
            <div className="sequence-table" role="table" aria-label="Alignment sequences">
              <div className="sequence-row sequence-header" role="row">
                <span className="sequence-analysis-cell" role="columnheader">
                  <input
                    type="checkbox"
                    checked={allMatchingSelected}
                    disabled={busy || matchingSequences.length === 0}
                    aria-label={allMatchingSelected
                      ? "Unselect every sequence matching the current filter"
                      : "Select every sequence matching the current filter"}
                    title={`${allMatchingSelected ? "Unselect" : "Select"} all ${matchingSequences.length.toLocaleString()} matching rows`}
                    onChange={toggleMatchingSelection}
                  />
                  <span>Analysis</span>
                </span>
                <span role="columnheader">Query / reference</span>
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
                    <span className="sequence-analysis-cell">
                      <input
                        type="checkbox"
                        checked={selectedSequenceIndices.has(sequence.index)}
                        disabled={busy}
                        aria-label={`Select ${sequence.name} for bulk query/reference assignment`}
                        onChange={() => toggleSequenceSelection(sequence.index)}
                      />
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
                    <label className="reference-group-field">
                      <span className="sr-only">Reference group for {sequence.name}</span>
                      <input
                        type="number"
                        min="0"
                        max="4294967295"
                        step="1"
                        inputMode="numeric"
                        value={referenceGroups.get(sequence.index) ?? ""}
                        placeholder="Query"
                        disabled={busy}
                        aria-label={`Reference group for ${sequence.name}; blank means query`}
                        onChange={(event) => onReferenceGroupChange(
                          sequence.index,
                          event.target.value === "" ? 0 : Number(event.target.value),
                        )}
                      />
                      <small>
                        {referenceGroups.has(sequence.index)
                          ? `Ref ${referenceGroups.get(sequence.index)}`
                          : "Query"}
                      </small>
                    </label>
                    <strong title={sequence.name}>{sequence.name}</strong>
                    <span>{integer.format(sequence.validSites)}</span>
                    <span>{percent.format(sequence.missingFraction)}</span>
                  </div>
                );
              })}
            </div>
            {matchingSequences.length > 500 ? (
              <p className="table-footnote">
                Showing the first 500 matching sequences for browser rendering performance. The header checkbox selects all {matchingSequences.length.toLocaleString()} matches, including rows not rendered.
              </p>
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
                  disabled={busy || curatedEnabledCount === 0}
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
              {eligibleSequenceCount < 3
                ? "Enable at least three unmasked sequences with enough usable sites to continue."
                : `${eligibleSequenceCount} sequences currently meet the primary-scan threshold.`}
            </span>
            <button
              className="button button-primary"
              type="button"
              onClick={onContinue}
              disabled={busy || eligibleSequenceCount < 3}
            >
              Set analysis options
            </button>
          </footer>
        </>
      ) : null}
    </section>
  );
}
