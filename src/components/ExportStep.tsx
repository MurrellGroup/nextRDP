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
  onRecombinationFree: () => void;
  onFragmented: () => void;
  onBack: () => void;
}

export function ExportStep({
  results,
  filename,
  onCsv,
  onProject,
  onRecombinationFree,
  onFragmented,
  onBack,
}: ExportStepProps) {
  const accepted = results.events.filter((event) => event.reviewState === "accepted").length;
  const reviewed = results.events.filter((event) => event.reviewState !== "unreviewed").length;
  const alignmentsReady = results.finalAlignmentReady;
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
        <span className="fidelity-badge">Session 5 snapshot</span>
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
          <button className="button button-primary" type="button" onClick={onProject}>
            <Download size={17} /> Download .rdpweb.json
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
          evidence sets and produces the accepted-event alignment variants. Unported native role
          methods, the remaining late correlation/tree filter stack, and parity validation remain
          explicit fidelity work.
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
