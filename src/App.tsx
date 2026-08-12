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
  EventEdit,
  ReviewState,
  ScanOptions,
  ScanProgress,
  ScanResults,
  SignalPlot,
  WorkflowStep,
} from "./lib/types";
import { RdpWorkerClient } from "./lib/wasmClient";

const initialOptions: ScanOptions = {
  circular: true,
  pValueCutoff: 0.05,
  correction: "bonferroni",
  windowSites: 30,
  maskedSequenceIndices: [],
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
  const [options, setOptions] = useState<ScanOptions>(initialOptions);
  const [progress, setProgress] = useState<ScanProgress>(initialProgress);
  const [results, setResults] = useState<ScanResults | null>(null);
  const [loading, setLoading] = useState(false);
  const [running, setRunning] = useState(false);
  const [reconciling, setReconciling] = useState(false);
  const [error, setError] = useState("");

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

  const activeSequenceCount = Math.max(0, (dataset?.sequenceCount ?? 0) - masked.size);
  const activeTripletCount = chooseThree(activeSequenceCount);

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
    setLoading(true);
    setError("");
    try {
      if (file.name.toLowerCase().endsWith(".json")) {
        const restored = await client.current.importProject(file);
        const restoredMask = new Set(
          restored.results?.maskedSequenceIndices ??
            restored.dataset.sequences.filter((sequence) => sequence.masked).map((sequence) => sequence.index),
        );
        setDataset(restored.dataset);
        setFilename(restored.sourceFilename);
        setFileSize(file.size);
        setMasked(restoredMask);
        setResults(restored.results);
        setOptions(
          restored.results
            ? {
                circular: restored.results.circular,
                pValueCutoff: restored.results.pValueCutoff,
                correction: restored.results.correction,
                windowSites: restored.results.windowSites,
                maskedSequenceIndices: [...restoredMask],
              }
            : { ...initialOptions, maskedSequenceIndices: [...restoredMask] },
        );
        setProgress(
          restored.results
            ? {
                state: "done",
                phase: "complete",
                processedTriplets: chooseThree(restored.dataset.sequenceCount - restoredMask.size),
                totalTriplets: chooseThree(restored.dataset.sequenceCount - restoredMask.size),
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
                totalTriplets: chooseThree(restored.dataset.sequenceCount - restoredMask.size),
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
      setOptions({ ...initialOptions, maskedSequenceIndices: [...initialMask] });
      setProgress({ ...initialProgress, totalTriplets: summary.tripletCount });
      setResults(null);
      setStep("dataset");
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    } finally {
      setLoading(false);
    }
  };

  const changeMask = (index: number, shouldMask: boolean) => {
    const next = new Set(masked);
    if (shouldMask) next.add(index);
    else next.delete(index);
    const nextActiveCount = Math.max(0, (dataset?.sequenceCount ?? 0) - next.size);
    setMasked(next);
    setProgress({ ...initialProgress, totalTriplets: chooseThree(nextActiveCount) });
    setResults(null);
  };

  const changeOptions = (next: ScanOptions) => {
    setOptions(next);
    setProgress({ ...initialProgress, totalTriplets: activeTripletCount });
    setResults(null);
    setError("");
  };

  const startScan = async () => {
    if (!client.current || !dataset) return;
    setError("");
    setRunning(true);
    setResults(null);
    const scanOptions = { ...options, maskedSequenceIndices: [...masked].sort((a, b) => a - b) };
    setProgress({
      ...initialProgress,
      state: "running",
      totalTriplets: activeTripletCount,
    });
    try {
      const value = await client.current.scan(scanOptions);
      setResults(value);
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
      if (updated) setResults(updated);
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
    try {
      const project = await client.current.exportProject();
      downloadBlob(
        new Blob([project], { type: "application/json;charset=utf-8" }),
        `${safeStem(filename)}.rdpweb.json`,
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
          <span className="session-pill">Fidelity port · session 5</span>
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
            busy={loading}
            onLoad={loadAlignment}
            onMaskChange={changeMask}
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
            onReviewState={setEventReviewState}
            onUpdateEvent={updateEvent}
            onUpdateEventGroup={updateEventGroup}
            onReconcileAfter={reconcileAfter}
            reconciling={reconciling}
            onSaveProject={exportProject}
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
            onRecombinationFree={exportRecombinationFree}
            onFragmented={exportFragmented}
            onBack={() => go("review")}
          />
        ) : null}
      </main>
    </div>
  );
}
