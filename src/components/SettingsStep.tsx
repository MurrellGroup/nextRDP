import { Check, CircleDot, Info, LockKeyhole, ScanSearch } from "lucide-react";

import type { ScanOptions } from "../lib/types";

interface SettingsStepProps {
  options: ScanOptions;
  sequenceCount: number;
  tripletCount: number;
  onChange: (options: ScanOptions) => void;
  onBack: () => void;
  onContinue: () => void;
}

const methods = [
  {
    name: "RDP",
    description: "Original informative-site triplet scan with binomial significance.",
    state: "ready",
  },
  {
    name: "GENECONV",
    description: "Fragment scoring and Karlin–Altschul / permutation significance.",
    state: "porting",
  },
  {
    name: "MAXCHI",
    description: "Variable-site maximum χ² breakpoint scan.",
    state: "porting",
  },
  {
    name: "BOOTSCAN",
    description: "Sliding-window phylogenetic support; secondary by default.",
    state: "queued",
  },
  {
    name: "CHIMAERA",
    description: "Two-variable-site χ² scan derived from MAXCHI.",
    state: "queued",
  },
  {
    name: "SISCAN",
    description: "Sister-scanning permutation test; secondary by default.",
    state: "queued",
  },
  {
    name: "3SEQ",
    description: "Exact triplet test and recombinant tract inference.",
    state: "queued",
  },
] as const;

export function SettingsStep({
  options,
  sequenceCount,
  tripletCount,
  onChange,
  onBack,
  onContinue,
}: SettingsStepProps) {
  const set = <Key extends keyof ScanOptions>(key: Key, value: ScanOptions[Key]) => {
    onChange({ ...options, [key]: value });
  };

  return (
    <section className="step-page" aria-labelledby="settings-title">
      <header className="page-heading">
        <div>
          <span className="eyebrow">02 · Settings</span>
          <h1 id="settings-title">Design the preliminary scan</h1>
          <p>
            Defaults follow the RDP5 manual: an exploratory, corrected screen first; detailed
            checking and hypothesis refinement afterwards.
          </p>
        </div>
        <div className="dataset-chip">
          <ScanSearch size={18} />
          <span>
            <strong>{sequenceCount} active sequences</strong>
            {tripletCount.toLocaleString()} unique triplets
          </span>
        </div>
      </header>

      <div className="settings-layout">
        <div className="settings-main">
          <div className="content-card">
            <div className="card-heading">
              <span className="eyebrow">Analysis scheme</span>
              <h2>Automated exploratory analysis</h2>
              <p>
                Every eligible triplet is screened without assuming a pre-defined non-recombinant
                reference set.
              </p>
            </div>
            <div className="segmented" role="radiogroup" aria-label="Sequence topology">
              <button
                type="button"
                role="radio"
                aria-checked={!options.circular}
                className={!options.circular ? "is-selected" : ""}
                onClick={() => set("circular", false)}
              >
                <CircleDot size={16} /> Linear sequences
              </button>
              <button
                type="button"
                role="radio"
                aria-checked={options.circular}
                className={options.circular ? "is-selected" : ""}
                onClick={() => set("circular", true)}
              >
                <CircleDot size={16} /> Circular sequences
              </button>
            </div>
          </div>

          <div className="content-card">
            <div className="card-heading split-heading">
              <div>
                <span className="eyebrow">Primary methods</span>
                <h2>Signal detection panel</h2>
              </div>
              <span className="fidelity-badge">1 of 7 ported</span>
            </div>
            <div className="method-grid">
              {methods.map((method) => (
                <article className={`method-card method-${method.state}`} key={method.name}>
                  <div className="method-state">
                    {method.state === "ready" ? <Check size={16} /> : <LockKeyhole size={15} />}
                  </div>
                  <div>
                    <h3>{method.name}</h3>
                    <p>{method.description}</p>
                    <span>{method.state === "ready" ? "Included in this scan" : method.state}</span>
                  </div>
                </article>
              ))}
            </div>
          </div>
        </div>

        <aside className="settings-side">
          <div className="content-card sticky-card">
            <div className="card-heading">
              <span className="eyebrow">Statistical controls</span>
              <h2>RDP method</h2>
            </div>
            <label className="field">
              <span>Highest acceptable p-value</span>
              <input
                type="number"
                min="0.0000000001"
                max="1"
                step="0.001"
                value={options.pValueCutoff}
                onChange={(event) => set("pValueCutoff", Number(event.target.value))}
              />
              <small>Applied after the selected multiple-comparison correction.</small>
            </label>
            <label className="field">
              <span>Multiple comparisons</span>
              <select
                value={options.correction}
                onChange={(event) =>
                  set("correction", event.target.value as ScanOptions["correction"])
                }
              >
                <option value="bonferroni">Bonferroni correction</option>
                <option value="none">No correction</option>
              </select>
              <small>Bonferroni is the RDP5 default for exploratory scans.</small>
            </label>
            <label className="field">
              <span>RDP window</span>
              <div className="input-suffix">
                <input
                  type="number"
                  min="5"
                  max="1001"
                  step="1"
                  value={options.windowSites}
                  onChange={(event) => set("windowSites", Number(event.target.value))}
                />
                <span>variable sites</span>
              </div>
              <small>The supplied RDP5 default is 30.</small>
            </label>
            <div className="inline-note">
              <Info size={17} />
              <p>
                Smaller windows increase sensitivity to short tracts and noise; larger windows
                trade sensitivity for stability.
              </p>
            </div>
          </div>
        </aside>
      </div>

      <footer className="step-actions">
        <button className="button button-quiet" type="button" onClick={onBack}>
          Back to dataset
        </button>
        <button className="button button-primary" type="button" onClick={onContinue}>
          Review scan plan
        </button>
      </footer>
    </section>
  );
}
