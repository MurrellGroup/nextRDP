import { Check, CircleDot, Info, LockKeyhole, ScanSearch } from "lucide-react";

import type { ScanOptions } from "../lib/types";

interface SettingsStepProps {
  options: ScanOptions;
  sequenceCount: number;
  tripletCount: number;
  exploratoryTripletCount: number;
  queryReferenceTripletCount: number;
  querySequenceCount: number;
  referenceSequenceCount: number;
  referenceGroupCount: number;
  queryReferenceCorrectionTestCount: number;
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
    description: "Six signed fragment tracks with source mismatch scoring, Karlin–Altschul tails, and overlap filtering.",
    state: "ready",
  },
  {
    name: "BOOTSCAN",
    description: "Seeded sliding-window bootstrap distances, strict closest-pair support, and binomial tract significance.",
    state: "ready",
  },
  {
    name: "MAXCHI",
    description: "Raw χ² peak ordering, tract-side matching, peak destruction/retry, and cyclic discovery.",
    state: "ready",
  },
  {
    name: "CHIMAERA",
    description: "Target-rotated information-rich binary χ² profiles with source-shaped tract inference.",
    state: "ready",
  },
  {
    name: "SISCAN",
    description: "Sister-scanning permutation test; secondary by default.",
    state: "queued",
  },
  {
    name: "3SEQ",
    description: "Target-rotated hypergeometric random walks with exact tails and source circular tract inference.",
    state: "ready",
  },
] as const;

export function SettingsStep({
  options,
  sequenceCount,
  tripletCount,
  exploratoryTripletCount,
  queryReferenceTripletCount,
  querySequenceCount,
  referenceSequenceCount,
  referenceGroupCount,
  queryReferenceCorrectionTestCount,
  onChange,
  onBack,
  onContinue,
}: SettingsStepProps) {
  const schemeValid = options.analysisMode === "exploratory" || (
    querySequenceCount >= 1 && referenceSequenceCount >= 2 &&
    referenceGroupCount >= 2 && queryReferenceTripletCount > 0
  );
  const settingsValid = sequenceCount >= 3 && schemeValid &&
    Number.isFinite(options.pValueCutoff) &&
    options.pValueCutoff > 0 &&
    options.pValueCutoff <= 1 &&
    Number.isInteger(options.windowSites) &&
    options.windowSites >= 5 &&
    options.windowSites <= 1001 &&
    (!options.maxChiEnabled || (
      Number.isInteger(options.maxChiWindowSites) &&
      options.maxChiWindowSites >= 12 &&
      options.maxChiWindowSites <= 2000
    )) &&
    (!options.chimaeraEnabled || (
      Number.isInteger(options.chimaeraWindowSites) &&
      options.chimaeraWindowSites >= 12 &&
      options.chimaeraWindowSites <= 2000
    )) &&
    (!options.geneconvEnabled || (
      Number.isInteger(options.geneconvMismatchScale) &&
      options.geneconvMismatchScale >= 1 &&
      options.geneconvMismatchScale <= 1000 &&
      Number.isInteger(options.geneconvMaxOverlaps) &&
      options.geneconvMaxOverlaps >= 1 &&
      options.geneconvMaxOverlaps <= 100
    )) &&
    (!(options.bootscanPrimaryEnabled || options.bootscanSecondaryEnabled) || (
      Number.isInteger(options.bootscanWindowSites) &&
      options.bootscanWindowSites >= 5 &&
      options.bootscanWindowSites <= 5000 &&
      Number.isInteger(options.bootscanStepSites) &&
      options.bootscanStepSites >= 1 &&
      options.bootscanStepSites <= Math.floor(options.bootscanWindowSites / 2) &&
      Number.isInteger(options.bootscanBootstrapReplicates) &&
      options.bootscanBootstrapReplicates >= 10 &&
      options.bootscanBootstrapReplicates <= 1000 &&
      Number.isFinite(options.bootscanSupportCutoff) &&
      options.bootscanSupportCutoff >= 0.5 &&
      options.bootscanSupportCutoff <= 1 &&
      Number.isInteger(options.bootscanRandomSeed) &&
      options.bootscanRandomSeed > 0
    ));
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
            {tripletCount.toLocaleString()} {options.analysisMode === "query-reference"
              ? "constrained triplets"
              : "unique triplets"}
          </span>
        </div>
      </header>

      <div className="settings-layout">
        <div className="settings-main">
          <div className="content-card">
            <div className="card-heading">
              <span className="eyebrow">Analysis scheme</span>
              <h2>
                {options.analysisMode === "query-reference"
                  ? "Automated query vs reference"
                  : "Fully exploratory analysis"}
              </h2>
              <p>
                {options.analysisMode === "query-reference"
                  ? "Every triplet contains one query and two references drawn from different groups, as specified in the RDP5 manual."
                  : "Every eligible triplet is screened without assuming a pre-defined non-recombinant reference set."}
              </p>
            </div>
            <div className="segmented analysis-scheme" role="radiogroup" aria-label="Analysis scheme">
              <button
                type="button"
                role="radio"
                aria-checked={options.analysisMode === "exploratory"}
                className={options.analysisMode === "exploratory" ? "is-selected" : ""}
                onClick={() => set("analysisMode", "exploratory")}
              >
                <CircleDot size={16} />
                <span>
                  Fully exploratory
                  <small>{exploratoryTripletCount.toLocaleString()} triplets</small>
                </span>
              </button>
              <button
                type="button"
                role="radio"
                aria-checked={options.analysisMode === "query-reference"}
                className={options.analysisMode === "query-reference" ? "is-selected" : ""}
                onClick={() => set("analysisMode", "query-reference")}
              >
                <CircleDot size={16} />
                <span>
                  Query vs reference
                  <small>{queryReferenceTripletCount.toLocaleString()} constrained triplets</small>
                </span>
              </button>
            </div>
            {options.analysisMode === "query-reference" ? (
              <div className={`inline-note${schemeValid ? "" : " notice-amber"}`}>
                <Info size={17} />
                <p>
                  {querySequenceCount.toLocaleString()} queries · {referenceSequenceCount.toLocaleString()} references · {referenceGroupCount.toLocaleString()} groups. References in the same group are not paired. Role inference remains free to identify a reference as recombinant.
                  {schemeValid && options.correction === "bonferroni"
                    ? ` The supplied MakeAnalysisListQvR correction uses ${queryReferenceCorrectionTestCount.toLocaleString()} group-pair × query opportunities after the native cap; progress still reports every scheduled reference-record triplet.`
                    : ""}
                  {!schemeValid ? " Return to the dataset and define at least one query plus two enabled reference groups." : ""}
                </p>
              </div>
            ) : null}
            <div className="card-heading topology-heading">
              <span className="eyebrow">Sequence topology</span>
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
              <span className="fidelity-badge">RDP + GENECONV + BootScan + MaxChi + CHIMAERA + 3SEQ discovery</span>
            </div>
            <div className="method-grid">
              {methods.map((method) => (
                <article className={`method-card method-${method.state}`} key={method.name}>
                  <div className="method-state">
                    {method.state === "ready"
                      ? <Check size={16} />
                      : <LockKeyhole size={15} />}
                  </div>
                  <div>
                    <h3>{method.name}</h3>
                    <p>{method.description}</p>
                    <span>
                      {method.state === "ready"
                        ? (method.name === "MAXCHI" && !options.maxChiEnabled) ||
                          (method.name === "CHIMAERA" && !options.chimaeraEnabled) ||
                          (method.name === "GENECONV" && !options.geneconvEnabled) ||
                          (method.name === "BOOTSCAN" && !options.bootscanPrimaryEnabled) ||
                          (method.name === "3SEQ" && !options.threeSeqEnabled)
                          ? "Available · disabled for this scan"
                          : "Included in event discovery"
                        : method.state}
                    </span>
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
              <h2>Primary methods</h2>
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
                <option value="bonferroni">Project correction (method-specific)</option>
                <option value="none">No correction</option>
              </select>
              <small>
                {options.analysisMode === "query-reference"
                  ? options.correction === "bonferroni"
                    ? options.threeSeqEnabled
                      ? "The opportunity count is reference-group pairs × queries: RDP-family methods multiply by it, while 3SEQ applies its supplied Dunn–Šidák form."
                      : "Query-vs-reference mode applies the supplied reference-group-pair × query correction rule."
                    : "No correction will be applied; the source opportunity count remains in progress and result metadata."
                  : options.threeSeqEnabled
                    ? "The RDP-family methods use the project Bonferroni factor; 3SEQ uses the supplied Dunn–Šidák form with the same opportunity count."
                    : "Bonferroni is the RDP5 default for exploratory scans."}
              </small>
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
            <label className="settings-toggle">
              <input
                type="checkbox"
                checked={options.geneconvEnabled}
                onChange={(event) => set("geneconvEnabled", event.target.checked)}
              />
              <span>
                <strong>Discover events with GENECONV</strong>
                <small>
                  Runs the supplied automated triplet KA workflow with ignored indels and stable
                  lowest-P overlap selection in every cyclic pass.
                </small>
              </span>
            </label>
            <label className="field">
              <span>GENECONV mismatch scale</span>
              <div className="input-suffix">
                <input
                  type="number"
                  min="1"
                  max="1000"
                  step="1"
                  value={options.geneconvMismatchScale}
                  disabled={!options.geneconvEnabled}
                  onChange={(event) => set("geneconvMismatchScale", Number(event.target.value))}
                />
                <span>G scale</span>
              </div>
              <small>The supplied automated default is 1.</small>
            </label>
            <label className="field">
              <span>GENECONV overlapping fragments</span>
              <div className="input-suffix">
                <input
                  type="number"
                  min="1"
                  max="100"
                  step="1"
                  value={options.geneconvMaxOverlaps}
                  disabled={!options.geneconvEnabled}
                  onChange={(event) => set("geneconvMaxOverlaps", Number(event.target.value))}
                />
                <span>per site</span>
              </div>
              <small>The supplied automated default accepts one covering fragment per polymorphic site.</small>
            </label>
            <label className="settings-toggle">
              <input
                type="checkbox"
                checked={options.bootscanPrimaryEnabled}
                onChange={(event) => set("bootscanPrimaryEnabled", event.target.checked)}
              />
              <span>
                <strong>Discover events with BootScan</strong>
                <small>
                  Runs the supplied automated distance-mode BSXoverR path in every cyclic pass.
                  Pair/window bootstrap distances are reused across triplets in a bounded cache.
                </small>
              </span>
            </label>
            <label className="settings-toggle">
              <input
                type="checkbox"
                checked={options.bootscanSecondaryEnabled}
                onChange={(event) => set("bootscanSecondaryEnabled", event.target.checked)}
              />
              <span>
                <strong>Use BootScan for late corroboration</strong>
                <small>
                  Rechecks representative and final-list triplets without moving reconciled event
                  boundaries. This remains independent of primary BootScan discovery.
                </small>
              </span>
            </label>
            <label className="field">
              <span>BootScan window</span>
              <div className="input-suffix">
                <input
                  type="number"
                  min="5"
                  max="5000"
                  step="1"
                  value={options.bootscanWindowSites}
                  disabled={!options.bootscanPrimaryEnabled && !options.bootscanSecondaryEnabled}
                  onChange={(event) => set("bootscanWindowSites", Number(event.target.value))}
                />
                <span>nucleotide sites</span>
              </div>
              <small>The supplied distance-mode default is 200.</small>
            </label>
            <label className="field">
              <span>BootScan step</span>
              <div className="input-suffix">
                <input
                  type="number"
                  min="1"
                  max={Math.max(1, Math.floor(options.bootscanWindowSites / 2))}
                  step="1"
                  value={options.bootscanStepSites}
                  disabled={!options.bootscanPrimaryEnabled && !options.bootscanSecondaryEnabled}
                  onChange={(event) => set("bootscanStepSites", Number(event.target.value))}
                />
                <span>sites</span>
              </div>
              <small>The supplied default is 20.</small>
            </label>
            <label className="field">
              <span>Bootstrap replicates</span>
              <input
                type="number"
                min="10"
                max="1000"
                step="10"
                value={options.bootscanBootstrapReplicates}
                disabled={!options.bootscanPrimaryEnabled && !options.bootscanSecondaryEnabled}
                onChange={(event) => set("bootscanBootstrapReplicates", Number(event.target.value))}
              />
              <small>Replicate zero is the unresampled window, matching SEQBOOT2.</small>
            </label>
            <label className="field">
              <span>Bootstrap support cutoff</span>
              <input
                type="number"
                min="0.5"
                max="1"
                step="0.01"
                value={options.bootscanSupportCutoff}
                disabled={!options.bootscanPrimaryEnabled && !options.bootscanSecondaryEnabled}
                onChange={(event) => set("bootscanSupportCutoff", Number(event.target.value))}
              />
              <small>The manual describes approximately 70% as the default detection threshold.</small>
            </label>
            <label className="field">
              <span>BootScan random seed</span>
              <input
                type="number"
                min="1"
                max="4294967295"
                step="1"
                value={options.bootscanRandomSeed}
                disabled={!options.bootscanPrimaryEnabled && !options.bootscanSecondaryEnabled}
                onChange={(event) => set("bootscanRandomSeed", Number(event.target.value))}
              />
              <small>The supplied default is 3, using the Microsoft 15-bit rand sequence.</small>
            </label>
            <label className="settings-toggle">
              <input
                type="checkbox"
                checked={options.threeSeqEnabled}
                onChange={(event) => set("threeSeqEnabled", event.target.checked)}
              />
              <span>
                <strong>Discover events with 3SEQ</strong>
                <small>
                  Rotates all three candidate recombinants, computes exact bounded
                  hypergeometric-walk tails, and uses the supplied Siegmund fallback only for
                  larger profiles.
                </small>
              </span>
            </label>
            <label className="settings-toggle">
              <input
                type="checkbox"
                checked={options.maxChiEnabled}
                onChange={(event) => set("maxChiEnabled", event.target.checked)}
              />
              <span>
                <strong>Discover events with MaxChi</strong>
                <small>
                  Runs the supplied triplet workflow in the same strongest-first cyclic passes as RDP.
                </small>
              </span>
            </label>
            <label className="field">
              <span>MaxChi window</span>
              <div className="input-suffix">
                <input
                  type="number"
                  min="12"
                  max="2000"
                  step="2"
                  value={options.maxChiWindowSites}
                  disabled={!options.maxChiEnabled}
                  onChange={(event) => set("maxChiWindowSites", Number(event.target.value))}
                />
                <span>variable sites</span>
              </div>
              <small>The supplied default is 70, split into two half-windows.</small>
            </label>
            <label className="settings-toggle">
              <input
                type="checkbox"
                checked={options.chimaeraEnabled}
                onChange={(event) => set("chimaeraEnabled", event.target.checked)}
              />
              <span>
                <strong>Discover events with CHIMAERA</strong>
                <small>
                  Rotates each triplet member through the candidate-recombinant role and scans its
                  target-specific parent-match string in the shared cyclic scheduler.
                </small>
              </span>
            </label>
            <label className="field">
              <span>CHIMAERA window</span>
              <div className="input-suffix">
                <input
                  type="number"
                  min="12"
                  max="2000"
                  step="2"
                  value={options.chimaeraWindowSites}
                  disabled={!options.chimaeraEnabled}
                  onChange={(event) => set("chimaeraWindowSites", Number(event.target.value))}
                />
                <span>information-rich sites</span>
              </div>
              <small>The supplied default is 60, split into two half-windows.</small>
            </label>
            <label className="settings-toggle">
              <input
                type="checkbox"
                checked={options.polishBreakpoints}
                onChange={(event) => set("polishBreakpoints", event.target.checked)}
              />
              <span>
                <strong>Polish breakpoints with BURT</strong>
                <small>
                  Enabled by default in the supplied RDP5 workflow. Runs the seeded BenHMM
                  confidence pass and may reposition automatically detected breakpoints.
                </small>
              </span>
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
        <button
          className="button button-primary"
          type="button"
          onClick={onContinue}
          disabled={!settingsValid}
        >
          Review scan plan
        </button>
      </footer>
    </section>
  );
}
