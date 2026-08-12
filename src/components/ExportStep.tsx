import {
  ArrowLeft,
  CheckCircle2,
  Dna,
  Download,
  FileArchive,
  FileJson,
  FileSpreadsheet,
  Scissors,
} from "lucide-react";

import type { ScanResults } from "../lib/types";

interface ExportStepProps {
  results: ScanResults;
  filename: string;
  onCsv: () => void;
  onProject: () => void;
  checkpointDirty: boolean;
  checkpointSaving: boolean;
  onFullAlignment: () => void;
  onEnabledSequences: () => void;
  onMaskedOrDisabledSequences: () => void;
  onRecombinantSequencesRemoved: () => void;
  onRecombinantColumnsRemoved: () => void;
  onRecombinationFree: () => void;
  onFragmented: () => void;
  onBack: () => void;
}

export function ExportStep({
  results,
  filename,
  onCsv,
  onProject,
  checkpointDirty,
  checkpointSaving,
  onFullAlignment,
  onEnabledSequences,
  onMaskedOrDisabledSequences,
  onRecombinantSequencesRemoved,
  onRecombinantColumnsRemoved,
  onRecombinationFree,
  onFragmented,
  onBack,
}: ExportStepProps) {
  const accepted = results.events.filter((event) => event.reviewState === "accepted").length;
  const reviewed = results.events.filter((event) => event.reviewState !== "unreviewed").length;
  const alignmentsReady = results.finalAlignmentReady;
  const excludedSequenceCount =
    results.maskedSequenceIndices.length + results.disabledSequenceIndices.length;
  return (
    <section className="step-page export-page" aria-labelledby="export-title">
      <header className="page-heading">
        <div>
          <span className="eyebrow">05 · Export</span>
          <h1 id="export-title">Save a reproducible snapshot</h1>
          <p>
            Keep the alignment, scan settings, event sets, trace evidence, and review decisions together;
            use CSV only as the human-readable summary.
          </p>
        </div>
        <div className="dataset-chip">
          <FileArchive size={18} />
          <span>
            <strong>{filename}</strong>
            {reviewed} of {results.events.length} events reviewed
          </span>
        </div>
      </header>

      <div className="export-summary content-card">
        <div>
          <CheckCircle2 size={22} />
          <span>
            <strong>Cyclic three-set analysis captured</strong>
            {results.events.length} events across {results.scanRounds} full passes · {results.workingFragmentSequenceCount} working fragments · {results.signals.length} retained signals · {accepted} accepted
          </span>
        </div>
        <span className="fidelity-badge">Session 9 snapshot</span>
      </div>

      <div className="export-grid">
        <article className="export-card">
          <span className="export-icon export-json"><FileJson size={24} /></span>
          <div>
            <span className="eyebrow">Recommended</span>
            <h2>RDP Web project</h2>
            <p>
              Reloadable JSON with the alignment, signals, three role hypotheses, correlation
              evidence, automatic/current groups, edits, and review state.
            </p>
          </div>
          <button
            className="button button-primary"
            type="button"
            disabled={checkpointSaving}
            onClick={onProject}
          >
            <Download size={17} />
            {checkpointSaving
              ? "Saving checkpoint…"
              : checkpointDirty
                ? "Download latest .rdpweb.json"
                : "Download .rdpweb.json"}
          </button>
        </article>

        <article className="export-card">
          <span className="export-icon export-csv"><FileSpreadsheet size={24} /></span>
          <div>
            <span className="eyebrow">Tabular summary</span>
            <h2>Event table</h2>
            <p>
              One row per event with all three evidence sets, automatic/current groups, early
              FinalTrim diagnostics, role votes, traces, and review decision.
            </p>
          </div>
          <button className="button button-secondary" type="button" onClick={onCsv}>
            <Download size={17} /> Download .csv
          </button>
        </article>

        <article className="export-card">
          <span className="export-icon export-fasta"><Dna size={23} /></span>
          <div>
            <span className="eyebrow">Loaded alignment</span>
            <h2>Entire alignment</h2>
            <p>
              Save every normalized aligned row as FASTA without applying sequence curation or event decisions.
            </p>
          </div>
          <button className="button button-secondary" type="button" onClick={onFullAlignment}>
            <Download size={17} /> Download full .fasta
          </button>
        </article>

        <article className="export-card">
          <span className="export-icon export-sequences"><Dna size={23} /></span>
          <div>
            <span className="eyebrow">Sequence curation</span>
            <h2>Enabled sequences only</h2>
            <p>
              Save only rows that entered the primary exploratory screen; masked and disabled rows are omitted unchanged.
            </p>
          </div>
          <button className="button button-secondary" type="button" onClick={onEnabledSequences}>
            <Download size={17} /> Download enabled-only .fasta
          </button>
        </article>

        <article className={excludedSequenceCount > 0 ? "export-card" : "export-card is-locked"}>
          <span className="export-icon export-columns"><Dna size={23} /></span>
          <div>
            <span className="eyebrow">Sequence curation</span>
            <h2>Masked or disabled only</h2>
            <p>
              Save the complementary curated rows: masked secondary-analysis sequences plus disabled tree-context sequences.
            </p>
          </div>
          <button
            className={excludedSequenceCount > 0 ? "button button-secondary" : "button button-disabled"}
            type="button"
            disabled={excludedSequenceCount === 0}
            onClick={onMaskedOrDisabledSequences}
          >
            <Download size={17} /> Download excluded rows .fasta
          </button>
        </article>

        <article className={alignmentsReady ? "export-card" : "export-card is-locked"}>
          <span className="export-icon export-sequences"><Dna size={23} /></span>
          <div>
            <span className="eyebrow">Accepted events</span>
            <h2>Remove recombinant sequences</h2>
            <p>
              Omit every sequence in an accepted current co-recombinant group and retain all other aligned records unchanged.
            </p>
          </div>
          <button
            className={alignmentsReady ? "button button-secondary" : "button button-disabled"}
            type="button"
            disabled={!alignmentsReady}
            onClick={onRecombinantSequencesRemoved}
          >
            <Download size={17} /> Download sequence-filtered .fasta
          </button>
        </article>

        <article className={alignmentsReady ? "export-card" : "export-card is-locked"}>
          <span className="export-icon export-columns"><Scissors size={23} /></span>
          <div>
            <span className="eyebrow">Accepted events</span>
            <h2>Remove recombinant columns</h2>
            <p>
              Delete the union of all columns covered by accepted event tracts from every sequence, producing a shorter alignment.
            </p>
          </div>
          <button
            className={alignmentsReady ? "button button-secondary" : "button button-disabled"}
            type="button"
            disabled={!alignmentsReady}
            onClick={onRecombinantColumnsRemoved}
          >
            <Download size={17} /> Download column-filtered .fasta
          </button>
        </article>

        <article className={alignmentsReady ? "export-card" : "export-card is-locked"}>
          <span className="export-icon export-fasta"><Dna size={23} /></span>
          <div>
            <span className="eyebrow">Accepted events</span>
            <h2>Tract-masked alignment</h2>
            <p>
              Replace accepted recombinant tracts with gaps in each complete co-recombinant group while preserving coordinates.
            </p>
          </div>
          <button
            className={alignmentsReady ? "button button-secondary" : "button button-disabled"}
            type="button"
            disabled={!alignmentsReady}
            onClick={onRecombinationFree}
          >
            <Download size={17} /> Download masked .fasta
          </button>
        </article>

        <article className={alignmentsReady ? "export-card" : "export-card is-locked"}>
          <span className="export-icon export-fragments"><Scissors size={23} /></span>
          <div>
            <span className="eyebrow">Manual workflow</span>
            <h2>Mosaic fragments</h2>
            <p>
              Keep the tract-masked originals and append aligned fragment-only copies in accepted event order.
            </p>
          </div>
          <button
            className={alignmentsReady ? "button button-secondary" : "button button-disabled"}
            type="button"
            disabled={!alignmentsReady}
            onClick={onFragmented}
          >
            <Download size={17} /> Download fragments .fasta
          </button>
        </article>
      </div>

      {!alignmentsReady ? (
        <div className="notice notice-amber">
          <Scissors size={18} />
          <p>
            Finish every event decision and reconcile any corrected or rejected event before producing final alignment variants.
          </p>
        </div>
      ) : null}

      <div className="export-caveat">
        <strong>Interpretation boundary</strong>
        <p>
          This snapshot completes the manual’s detectable, distance-correlation, and bootstrap-tree
          evidence sets, exposes original-alignment breakpoint windows, and produces the
          sequence-curation subsets plus accepted-event alignment variants. A source-shaped MaxChi
          strongest-peak recheck now corroborates representative triplets and finalized distance
          lists, but it does not discover or move events. MaxChi exploratory discovery, the remaining
          method families, and native golden parity validation remain explicit fidelity work;
          BURT/BenHMM confidence is active but still carries that validation boundary.
        </p>
      </div>

      <footer className="step-actions">
        <button className="button button-quiet" type="button" onClick={onBack}>
          <ArrowLeft size={16} /> Continue review
        </button>
      </footer>
    </section>
  );
}
