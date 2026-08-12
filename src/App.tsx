import { Atom, Cpu, Menu, PanelLeftClose, ShieldCheck } from "lucide-react";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";

import { DatasetStep } from "./components/DatasetStep";
import { ExportStep } from "./components/ExportStep";
import { ReviewStep } from "./components/ReviewStep";
import { ScanStep } from "./components/ScanStep";
import { SettingsStep } from "./components/SettingsStep";
import { WorkflowNav } from "./components/WorkflowNav";
import { downloadBlob, safeStem } from "./lib/download";
import type {
  DatasetSummary,
  EventAlignmentView,
  EventTreeView,
  EventEdit,
  ReviewState,
  ScanOptions,
  ScanProgress,
  ScanResults,
  SequenceAnalysisState,
  SignalPlot,
  WorkflowStep,
} from "./lib/types";
import { RdpWorkerClient } from "./lib/wasmClient";

const initialOptions: ScanOptions = {
  circular: true,
  pValueCutoff: 0.05,
  correction: "bonferroni",
  windowSites: 30,
  polishBreakpoints: true,
  maskedSequenceIndices: [],
  disabledSequenceIndices: [],
};

const initialProgress: ScanProgress = {
  state: "idle",
  phase: "primary",
  processedTriplets: 0,
  totalTriplets: 0,
  cumulativeTriplets: 0,
  scanRound: 1,
  fixedEventCount: 0,
  signalCount: 0,
  eventCount: 0,
  cycleTermination: "not-started",
  fraction: 0,
};

function chooseThree(count: number): number {
  return count < 3 ? 0 : (count * (count - 1) * (count - 2)) / 6;
}

type EngineState =
  | { status: "loading"; message: string; threaded: false }
  | { status: "ready"; message: string; threaded: boolean }
  | { status: "error"; message: string; threaded: false };

export function App() {
  const client = useRef<RdpWorkerClient | null>(null);
  const [step, setStep] = useState<WorkflowStep>("dataset");
  const [navOpen, setNavOpen] = useState(false);
  const [engine, setEngine] = useState<EngineState>({
    status: "loading",
    message: "Initialising the analysis worker…",
    threaded: false,
  });
  const [dataset, setDataset] = useState<DatasetSummary | null>(null);
  const [filename, setFilename] = useState("");
  const [fileSize, setFileSize] = useState(0);
  const [masked, setMasked] = useState<Set<number>>(new Set());
  const [disabled, setDisabled] = useState<Set<number>>(new Set());
  const [options, setOptions] = useState<ScanOptions>(initialOptions);
  const [progress, setProgress] = useState<ScanProgress>(initialProgress);
  const [results, setResults] = useState<ScanResults | null>(null);
  const [loading, setLoading] = useState(false);
  const [running, setRunning] = useState(false);
  const [reconciling, setReconciling] = useState(false);
  const [checkpointDirty, setCheckpointDirty] = useState(false);
  const [checkpointSaving, setCheckpointSaving] = useState(false);
  const [error, setError] = useState("");
  const hasUnsavedCheckpoint = checkpointDirty && results !== null;

  useEffect(() => {
    const worker = new RdpWorkerClient();
    client.current = worker;
    const removeProgress = worker.onProgress(setProgress);
    let live = true;
    worker
      .init()
      .then(({ threaded, version }) => {
        if (!live) return;
        setEngine({
          status: "ready",
          threaded,
          message: `${threaded ? "Thread-capable" : "Single-worker"} WASM · engine ${version}`,
        });
      })
      .catch((caught: unknown) => {
        if (!live) return;
        setEngine({
          status: "error",
          threaded: false,
          message: caught instanceof Error ? caught.message : String(caught),
        });
      });
    return () => {
      live = false;
      removeProgress();
      worker.dispose();
      client.current = null;
    };
  }, []);

  useEffect(() => {
    if (!hasUnsavedCheckpoint && !running && !reconciling) return;
    const warnBeforeUnload = (event: BeforeUnloadEvent) => {
      event.preventDefault();
      event.returnValue = "";
    };
    window.addEventListener("beforeunload", warnBeforeUnload);
    return () => window.removeEventListener("beforeunload", warnBeforeUnload);
  }, [hasUnsavedCheckpoint, reconciling, running]);

  const activeSequenceCount = Math.max(
    0,
    (dataset?.sequenceCount ?? 0) - masked.size - disabled.size,
  );
  const activeTripletCount = chooseThree(activeSequenceCount);

  const confirmDiscardUnsavedAnalysis = () =>
    !hasUnsavedCheckpoint ||
    window.confirm(
      "This analysis has changes that are not in a downloaded project checkpoint. Continue and discard them?",
    );

  const enabledSteps = useMemo(() => {
    if (running) return new Set<WorkflowStep>(["scan"]);
    if (reconciling) return new Set<WorkflowStep>(["review"]);
    const enabled = new Set<WorkflowStep>(["dataset"]);
    if (dataset) {
      enabled.add("settings");
      enabled.add("scan");
    }
    if (results) {
      enabled.add("review");
      enabled.add("export");
    }
    return enabled;
  }, [dataset, reconciling, results, running]);

  const completedSteps = useMemo(() => {
    const completed = new Set<WorkflowStep>();
    if (dataset) completed.add("dataset");
    if (progress.state === "running" || results) completed.add("settings");
    if (results) completed.add("scan");
    if (results?.events.some((event) => event.reviewState !== "unreviewed")) {
      completed.add("review");
    }
    return completed;
  }, [dataset, progress.state, results]);

  const go = (next: WorkflowStep) => {
    if (!enabledSteps.has(next)) return;
    setStep(next);
    setNavOpen(false);
    window.scrollTo({ top: 0, behavior: "smooth" });
  };

  const loadAlignment = async (file: File) => {
    if (!client.current || engine.status !== "ready") return;
    if (!confirmDiscardUnsavedAnalysis()) return;
    setLoading(true);
    setError("");
    try {
      if (file.name.toLowerCase().endsWith(".json")) {
        const restored = await client.current.importProject(file);
        const restoredMask = new Set(
          restored.results?.maskedSequenceIndices ??
            restored.dataset.sequences.filter((sequence) => sequence.masked).map((sequence) => sequence.index),
        );
        const restoredDisabled = new Set(restored.results?.disabledSequenceIndices ?? []);
        restoredDisabled.forEach((index) => restoredMask.delete(index));
        setDataset(restored.dataset);
        setFilename(restored.sourceFilename);
        setFileSize(file.size);
        setMasked(restoredMask);
        setDisabled(restoredDisabled);
        setResults(restored.results);
        setCheckpointDirty(false);
        setOptions(
          restored.results
            ? {
                circular: restored.results.circular,
                pValueCutoff: restored.results.pValueCutoff,
                correction: restored.results.correction,
                windowSites: restored.results.windowSites,
                polishBreakpoints: restored.results.polishBreakpoints ?? true,
                maskedSequenceIndices: [...restoredMask],
                disabledSequenceIndices: [...restoredDisabled],
              }
            : {
                ...initialOptions,
                maskedSequenceIndices: [...restoredMask],
                disabledSequenceIndices: [...restoredDisabled],
              },
        );
        setProgress(
          restored.results
            ? {
                state: "done",
                phase: "complete",
                processedTriplets: chooseThree(
                  restored.dataset.sequenceCount - restoredMask.size - restoredDisabled.size,
                ),
                totalTriplets: chooseThree(
                  restored.dataset.sequenceCount - restoredMask.size - restoredDisabled.size,
                ),
                cumulativeTriplets: restored.results.cumulativeTriplets,
                scanRound: restored.results.scanRounds,
                fixedEventCount: 0,
                signalCount: restored.results.signals.length,
                eventCount: restored.results.events.length,
                cycleTermination: restored.results.cycleTermination,
                fraction: 1,
              }
            : {
                ...initialProgress,
                totalTriplets: chooseThree(
                  restored.dataset.sequenceCount - restoredMask.size - restoredDisabled.size,
                ),
              },
        );
        setStep(restored.results ? "review" : "dataset");
        return;
      }
      const summary = await client.current.load(file);
      setDataset(summary);
      setFilename(file.name);
      setFileSize(file.size);
      const initialMask = new Set(summary.sequences.filter((sequence) => sequence.masked).map((sequence) => sequence.index));
      setMasked(initialMask);
      setDisabled(new Set());
      setOptions({
        ...initialOptions,
        maskedSequenceIndices: [...initialMask],
        disabledSequenceIndices: [],
      });
      setProgress({ ...initialProgress, totalTriplets: summary.tripletCount });
      setResults(null);
      setCheckpointDirty(false);
      setStep("dataset");
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    } finally {
      setLoading(false);
    }
  };

  const changeSequenceState = (index: number, state: SequenceAnalysisState) => {
    if (!confirmDiscardUnsavedAnalysis()) return;
    const nextMasked = new Set(masked);
    const nextDisabled = new Set(disabled);
    nextMasked.delete(index);
    nextDisabled.delete(index);
    if (state === "masked") nextMasked.add(index);
    if (state === "disabled") nextDisabled.add(index);
    const nextActiveCount = Math.max(
      0,
      (dataset?.sequenceCount ?? 0) - nextMasked.size - nextDisabled.size,
    );
    setMasked(nextMasked);
    setDisabled(nextDisabled);
    setOptions((current) => ({
      ...current,
      maskedSequenceIndices: [...nextMasked].sort((left, right) => left - right),
      disabledSequenceIndices: [...nextDisabled].sort((left, right) => left - right),
    }));
    setProgress({ ...initialProgress, totalTriplets: chooseThree(nextActiveCount) });
    setResults(null);
    setCheckpointDirty(false);
  };

  const changeAllSequenceStates = (
    action: "auto-mask" | "enable-all" | "mask-all" | "disable-all",
  ) => {
    if (!dataset || !confirmDiscardUnsavedAnalysis()) return;
    const nextMasked = new Set<number>();
    const nextDisabled = new Set<number>();
    if (action === "auto-mask") {
      dataset.sequences.forEach((sequence) => {
        if (sequence.masked) nextMasked.add(sequence.index);
      });
    } else if (action === "mask-all") {
      dataset.sequences.forEach((sequence) => nextMasked.add(sequence.index));
    } else if (action === "disable-all") {
      dataset.sequences.forEach((sequence) => nextDisabled.add(sequence.index));
    }
    const nextActiveCount = Math.max(
      0,
      dataset.sequenceCount - nextMasked.size - nextDisabled.size,
    );
    setMasked(nextMasked);
    setDisabled(nextDisabled);
    setOptions((current) => ({
      ...current,
      maskedSequenceIndices: [...nextMasked].sort((left, right) => left - right),
      disabledSequenceIndices: [...nextDisabled].sort((left, right) => left - right),
    }));
    setProgress({ ...initialProgress, totalTriplets: chooseThree(nextActiveCount) });
    setResults(null);
    setCheckpointDirty(false);
    setError("");
  };

  const changeOptions = (next: ScanOptions) => {
    if (!confirmDiscardUnsavedAnalysis()) return;
    setOptions(next);
    setProgress({ ...initialProgress, totalTriplets: activeTripletCount });
    setResults(null);
    setCheckpointDirty(false);
    setError("");
  };

  const startScan = async () => {
    if (!client.current || !dataset) return;
    if (!confirmDiscardUnsavedAnalysis()) return;
    setError("");
    setRunning(true);
    setResults(null);
    setCheckpointDirty(false);
    const scanOptions = {
      ...options,
      maskedSequenceIndices: [...masked].sort((a, b) => a - b),
      disabledSequenceIndices: [...disabled].sort((a, b) => a - b),
    };
    setProgress({
      ...initialProgress,
      state: "running",
      totalTriplets: activeTripletCount,
    });
    try {
      const value = await client.current.scan(scanOptions);
      setResults(value);
      setCheckpointDirty(true);
      setProgress((current) => ({ ...current, state: "done", phase: "complete", fraction: 1 }));
    } catch (caught) {
      const message = caught instanceof Error ? caught.message : String(caught);
      if (!message.toLowerCase().includes("cancel")) setError(message);
      setProgress((current) => ({
        ...current,
        state: message.toLowerCase().includes("cancel") ? "cancelled" : "error",
      }));
    } finally {
      setRunning(false);
    }
  };

  const cancelScan = async () => {
    try {
      await client.current?.cancel();
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  const getPlot = useCallback(async (signalId: number): Promise<SignalPlot> => {
    if (!client.current) throw new Error("The analysis worker is not available.");
    return client.current.plot(signalId);
  }, []);

  const getEventAlignment = useCallback(
    async (eventId: number, flankSites: number, rowLimit: number): Promise<EventAlignmentView> => {
      if (!client.current) throw new Error("The analysis worker is not available.");
      return client.current.eventAlignment(eventId, flankSites, rowLimit);
    },
    [],
  );

  const getEventTrees = useCallback(async (eventId: number): Promise<EventTreeView> => {
    if (!client.current) throw new Error("The analysis worker is not available.");
    return client.current.eventTrees(eventId);
  }, []);

  const setEventReviewState = async (eventId: number, state: ReviewState) => {
    const previous = results;
    if (!previous) return;
    const changedEvent = previous.events.find((event) => event.id === eventId);
    const rejectedMarker =
      state === "rejected" && changedEvent?.tractErasedForDetection
        ? Math.min(previous.downstreamReconciliationRequiredAfter ?? eventId, eventId)
        : previous.downstreamReconciliationRequiredAfter;
    setResults({
      ...previous,
      downstreamReconciliationRequiredAfter: rejectedMarker,
      events: previous.events.map((event) =>
        event.id === eventId ? { ...event, reviewState: state } : event,
      ),
    });
    try {
      const updated = await client.current?.setEventReviewState(eventId, state);
      if (updated) {
        setResults(updated);
        setCheckpointDirty(true);
      }
    } catch (caught) {
      setResults(previous);
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  const updateEvent = async (eventId: number, edit: EventEdit) => {
    if (!client.current) return;
    setError("");
    try {
      setResults(await client.current.updateEvent(eventId, edit));
      setCheckpointDirty(true);
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
      throw caught;
    }
  };

  const updateEventGroup = async (
    eventId: number,
    sequenceIndices: number[],
    manualOverride = true,
  ) => {
    if (!client.current) return;
    setError("");
    try {
      setResults(await client.current.updateEventGroup(eventId, sequenceIndices, manualOverride));
      setCheckpointDirty(true);
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
      throw caught;
    }
  };

  const reconcileAfter = async (eventId: number) => {
    if (!client.current) return;
    setError("");
    setReconciling(true);
    try {
      setResults(await client.current.reconcileAfter(eventId));
      setCheckpointDirty(true);
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    } finally {
      setReconciling(false);
    }
  };

  const exportCsv = async () => {
    if (!client.current) return;
    try {
      const csv = await client.current.exportCsv();
      downloadBlob(new Blob([csv], { type: "text/csv;charset=utf-8" }), `${safeStem(filename)}-rdp-events.csv`);
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  const exportProject = async () => {
    if (!client.current) return;
    setCheckpointSaving(true);
    try {
      const project = await client.current.exportProject();
      downloadBlob(
        new Blob([project], { type: "application/json;charset=utf-8" }),
        `${safeStem(filename)}.rdpweb.json`,
      );
      setCheckpointDirty(false);
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    } finally {
      setCheckpointSaving(false);
    }
  };

  const exportEnabledSequences = async () => {
    if (!client.current) return;
    try {
      const fasta = await client.current.exportEnabledSequences(
        [...masked].sort((left, right) => left - right),
        [...disabled].sort((left, right) => left - right),
      );
      downloadBlob(
        new Blob([fasta], { type: "text/plain;charset=utf-8" }),
        `${safeStem(filename)}-enabled-sequences.fasta`,
      );
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  const exportFullAlignment = async () => {
    if (!client.current) return;
    try {
      const fasta = await client.current.exportEnabledSequences([], []);
      downloadBlob(
        new Blob([fasta], { type: "text/plain;charset=utf-8" }),
        `${safeStem(filename)}-full-alignment.fasta`,
      );
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  const exportMaskedOrDisabledSequences = async () => {
    if (!client.current) return;
    try {
      const fasta = await client.current.exportMaskedOrDisabledSequences(
        [...masked].sort((left, right) => left - right),
        [...disabled].sort((left, right) => left - right),
      );
      downloadBlob(
        new Blob([fasta], { type: "text/plain;charset=utf-8" }),
        `${safeStem(filename)}-masked-or-disabled-sequences.fasta`,
      );
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  const exportRecombinationFree = async () => {
    if (!client.current) return;
    try {
      const fasta = await client.current.exportRecombinationFree();
      downloadBlob(
        new Blob([fasta], { type: "text/plain;charset=utf-8" }),
        `${safeStem(filename)}-recombination-free.fasta`,
      );
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  const exportRecombinantSequencesRemoved = async () => {
    if (!client.current) return;
    try {
      const fasta = await client.current.exportRecombinantSequencesRemoved();
      downloadBlob(
        new Blob([fasta], { type: "text/plain;charset=utf-8" }),
        `${safeStem(filename)}-recombinant-sequences-removed.fasta`,
      );
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  const exportRecombinantColumnsRemoved = async () => {
    if (!client.current) return;
    try {
      const fasta = await client.current.exportRecombinantColumnsRemoved();
      downloadBlob(
        new Blob([fasta], { type: "text/plain;charset=utf-8" }),
        `${safeStem(filename)}-recombinant-columns-removed.fasta`,
      );
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  const exportFragmented = async () => {
    if (!client.current) return;
    try {
      const fasta = await client.current.exportFragmented();
      downloadBlob(
        new Blob([fasta], { type: "text/plain;charset=utf-8" }),
        `${safeStem(filename)}-mosaic-fragments.fasta`,
      );
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    }
  };

  return (
    <div className="app-shell">
      <header className="topbar">
        <button className="mobile-menu" type="button" onClick={() => setNavOpen((open) => !open)} aria-label="Toggle workflow navigation">
          {navOpen ? <PanelLeftClose /> : <Menu />}
        </button>
        <div className="brand">
          <span className="brand-mark"><Atom size={23} /></span>
          <span>
            <strong>RDP Web</strong>
            <small>Recombination detection, browser-native</small>
          </span>
        </div>
        <div className="topbar-status">
          <span className={`engine-pill engine-${engine.status}`} title={engine.message}>
            <Cpu size={15} />
            {engine.status === "ready" ? "WASM ready" : engine.status === "loading" ? "Starting engine" : "Source snapshot"}
          </span>
          {results ? (
            <span
              className={`checkpoint-pill${hasUnsavedCheckpoint ? " is-dirty" : " is-current"}`}
              aria-live="polite"
            >
              {checkpointSaving
                ? "Saving checkpoint…"
                : hasUnsavedCheckpoint
                  ? "Checkpoint needed"
                  : "Checkpoint current"}
            </span>
          ) : null}
          <span className="session-pill">Fidelity port · session 9</span>
        </div>
      </header>

      <aside className={`sidebar${navOpen ? " is-open" : ""}`}>
        <WorkflowNav
          current={step}
          enabled={enabledSteps}
          completed={completedSteps}
          onSelect={go}
        />
        <div className="sidebar-foot">
          <ShieldCheck size={17} />
          <span>
            <strong>Private analysis</strong>
            No sequence data leave the tab.
          </span>
        </div>
      </aside>

      <main className="main-content">
        {error && step !== "scan" ? (
          <div className="global-error" role="alert">
            <strong>RDP Web could not complete that action.</strong>
            <span>{error}</span>
            <button type="button" onClick={() => setError("")}>Dismiss</button>
          </div>
        ) : null}

        {step === "dataset" ? (
          <DatasetStep
            engineReady={engine.status === "ready"}
            engineMessage={engine.status === "error" ? engine.message : undefined}
            dataset={dataset}
            filename={filename}
            fileSize={fileSize}
            masked={masked}
            disabled={disabled}
            busy={loading}
            onLoad={loadAlignment}
            onSequenceStateChange={changeSequenceState}
            onAllSequenceStatesChange={changeAllSequenceStates}
            onExportFullAlignment={exportFullAlignment}
            onExportEnabledSequences={exportEnabledSequences}
            onExportMaskedOrDisabledSequences={exportMaskedOrDisabledSequences}
            onContinue={() => go("settings")}
          />
        ) : null}

        {step === "settings" && dataset ? (
          <SettingsStep
            options={options}
            sequenceCount={activeSequenceCount}
            tripletCount={activeTripletCount}
            onChange={changeOptions}
            onBack={() => go("dataset")}
            onContinue={() => go("scan")}
          />
        ) : null}

        {step === "scan" && dataset ? (
          <ScanStep
            options={options}
            sequenceCount={activeSequenceCount}
            tripletCount={activeTripletCount}
            progress={progress}
            running={running}
            error={error}
            hasResults={results !== null}
            onStart={startScan}
            onCancel={cancelScan}
            onBack={() => go("settings")}
            onReview={() => go("review")}
          />
        ) : null}

        {step === "review" && dataset && results ? (
          <ReviewStep
            results={results}
            alignmentLength={dataset.alignmentLength}
            sequences={dataset.sequences}
            onGetPlot={getPlot}
            onGetEventAlignment={getEventAlignment}
            onGetEventTrees={getEventTrees}
            onReviewState={setEventReviewState}
            onUpdateEvent={updateEvent}
            onUpdateEventGroup={updateEventGroup}
            onReconcileAfter={reconcileAfter}
            reconciling={reconciling}
            onSaveProject={exportProject}
            checkpointDirty={hasUnsavedCheckpoint}
            checkpointSaving={checkpointSaving}
            onBack={() => go("scan")}
            onExport={() => go("export")}
          />
        ) : null}

        {step === "export" && results ? (
          <ExportStep
            results={results}
            filename={filename}
            onCsv={exportCsv}
            onProject={exportProject}
            checkpointDirty={hasUnsavedCheckpoint}
            checkpointSaving={checkpointSaving}
            onFullAlignment={exportFullAlignment}
            onEnabledSequences={exportEnabledSequences}
            onMaskedOrDisabledSequences={exportMaskedOrDisabledSequences}
            onRecombinantSequencesRemoved={exportRecombinantSequencesRemoved}
            onRecombinantColumnsRemoved={exportRecombinantColumnsRemoved}
            onRecombinationFree={exportRecombinationFree}
            onFragmented={exportFragmented}
            onBack={() => go("review")}
          />
        ) : null}
      </main>
    </div>
  );
}
