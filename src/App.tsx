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
  analysisMode: "exploratory",
  circular: true,
  pValueCutoff: 0.05,
  correction: "bonferroni",
  windowSites: 30,
  maxChiEnabled: true,
  maxChiWindowSites: 70,
  chimaeraEnabled: true,
  chimaeraWindowSites: 60,
  geneconvEnabled: true,
  geneconvMismatchScale: 1,
  geneconvMaxOverlaps: 1,
  threeSeqEnabled: true,
  bootscanPrimaryEnabled: false,
  bootscanSecondaryEnabled: false,
  bootscanWindowSites: 200,
  bootscanStepSites: 20,
  bootscanBootstrapReplicates: 100,
  bootscanSupportCutoff: 0.7,
  bootscanRandomSeed: 3,
  polishBreakpoints: true,
  maskedSequenceIndices: [],
  disabledSequenceIndices: [],
  referenceGroupIndices: [],
};

const initialProgress: ScanProgress = {
  state: "idle",
  phase: "primary",
  processedTriplets: 0,
  totalTriplets: 0,
  correctionTests: 0,
  activeWorkingSequenceCount: 0,
  queryWorkingSequenceCount: 0,
  referenceWorkingSequenceCount: 0,
  activeReferenceGroupCount: 0,
  cumulativeTriplets: 0,
  tripletKernelEvaluations: 0,
  tripletSummariesReused: 0,
  cachedSignalsReused: 0,
  methodScansSkipped: 0,
  scanRound: 1,
  fixedEventCount: 0,
  signalCount: 0,
  eventCount: 0,
  maxChiProfilesScanned: 0,
  maxChiPeakAttempts: 0,
  maxChiCandidatesFound: 0,
  maxChiPeakLimitTriplets: 0,
  chimaeraProfilesScanned: 0,
  chimaeraPeakAttempts: 0,
  chimaeraCandidatesFound: 0,
  chimaeraPeakLimitTargets: 0,
  geneconvFragmentsScored: 0,
  geneconvQualifiedFragments: 0,
  geneconvCandidatesFound: 0,
  geneconvOverlapRejections: 0,
  geneconvNumericalFallbackTracks: 0,
  threeSeqProfilesScanned: 0,
  threeSeqExactEvaluations: 0,
  threeSeqApproximateEvaluations: 0,
  threeSeqCandidatesFound: 0,
  bootscanProfilesScanned: 0,
  bootscanCandidateRegionsScored: 0,
  bootscanCandidatesFound: 0,
  bootscanPairProfilesRequested: 0,
  bootscanPairProfileCacheHits: 0,
  bootscanPairProfileCacheMisses: 0,
  bootscanPairProfileCacheEvictions: 0,
  bootscanPairProfileCacheBytes: 0,
  bootscanPairProfileCachePeakBytes: 0,
  cycleTermination: "not-started",
  fraction: 0,
};

const NATIVE_CORRECTION_CAP = Math.floor((255 ** 4) / 2);

function chooseThree(count: number): number {
  return count < 3 ? 0 : (count * (count - 1) * (count - 2)) / 6;
}

function referenceGroupArray(sequenceCount: number, groups: Map<number, number>): number[] {
  return Array.from({ length: sequenceCount }, (_, index) =>
    normalizeReferenceGroup(groups.get(index) ?? 0)
  );
}

function normalizeReferenceGroup(value: number): number {
  return Number.isFinite(value)
    ? Math.max(0, Math.min(0xffff_ffff, Math.trunc(value)))
    : 0;
}

function inferReferenceGroups(dataset: DatasetSummary): Map<number, number> {
  const groups = new Map<number, number>();
  const groupIds = new Map<string, number>();
  dataset.sequences.forEach((sequence) => {
    // The manual documents names such as REF-A<sequence>. An unlabelled REF
    // record is treated as its own group so it can be paired with every other
    // ungrouped reference, matching the non-grouped query/reference workflow.
    const match = sequence.name.match(/(?:^|[\s|])REF(?:[-_:]?)([^<\s|:]*)/i);
    if (!match) return;
    const label = match[1]?.trim().toLocaleLowerCase();
    const key = label ? `label:${label}` : `sequence:${sequence.index}`;
    let group = groupIds.get(key);
    if (!group) {
      group = groupIds.size + 1;
      groupIds.set(key, group);
    }
    groups.set(sequence.index, group);
  });
  return groups;
}

function compactReferenceGroups(
  dataset: DatasetSummary,
  groups: Map<number, number>,
): Map<number, number> {
  const compactIds = new Map<number, number>();
  const compacted = new Map<number, number>();
  dataset.sequences.forEach((sequence) => {
    const originalGroup = normalizeReferenceGroup(groups.get(sequence.index) ?? 0);
    if (originalGroup === 0) return;
    let compactGroup = compactIds.get(originalGroup);
    if (!compactGroup) {
      compactGroup = compactIds.size + 1;
      compactIds.set(originalGroup, compactGroup);
    }
    compacted.set(sequence.index, compactGroup);
  });
  return compacted;
}

function queryReferencePlan(
  dataset: DatasetSummary | null,
  masked: Set<number>,
  disabled: Set<number>,
  groups: Map<number, number>,
  minimumValidSites: number,
) {
  let querySequenceCount = 0;
  let referenceSequenceCount = 0;
  const groupCounts = new Map<number, number>();
  dataset?.sequences.forEach((sequence) => {
    if (masked.has(sequence.index) || disabled.has(sequence.index) ||
        sequence.validSites < minimumValidSites) return;
    const group = normalizeReferenceGroup(groups.get(sequence.index) ?? 0);
    if (group === 0) {
      ++querySequenceCount;
      return;
    }
    ++referenceSequenceCount;
    groupCounts.set(group, (groupCounts.get(group) ?? 0) + 1);
  });
  const counts = [...groupCounts.values()];
  let crossGroupReferencePairs = 0;
  for (let first = 0; first < counts.length; ++first) {
    for (let second = first + 1; second < counts.length; ++second) {
      crossGroupReferencePairs += counts[first] * counts[second];
    }
  }
  return {
    querySequenceCount,
    referenceSequenceCount,
    referenceGroupCount: groupCounts.size,
    tripletCount: querySequenceCount * crossGroupReferencePairs,
    correctionTestCount: Math.min(
      NATIVE_CORRECTION_CAP,
      querySequenceCount * groupCounts.size * (groupCounts.size - 1) / 2,
    ),
  };
}

function scanEligibleSequenceCount(
  dataset: DatasetSummary | null,
  masked: Set<number>,
  disabled: Set<number>,
  minimumValidSites: number,
): number {
  return dataset?.sequences.reduce(
    (count, sequence) => count + Number(
      !masked.has(sequence.index) &&
      !disabled.has(sequence.index) &&
      sequence.validSites >= minimumValidSites,
    ),
    0,
  ) ?? 0;
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
  const [referenceGroups, setReferenceGroups] = useState<Map<number, number>>(new Map());
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

  const minimumWorkingSites = Math.max(5, options.windowSites);
  const activeSequenceCount = scanEligibleSequenceCount(
    dataset,
    masked,
    disabled,
    minimumWorkingSites,
  );
  const exploratoryTripletCount = chooseThree(activeSequenceCount);
  const referencePlan = useMemo(
    () => queryReferencePlan(
      dataset,
      masked,
      disabled,
      referenceGroups,
      minimumWorkingSites,
    ),
    [dataset, disabled, masked, minimumWorkingSites, referenceGroups],
  );
  const activeTripletCount = options.analysisMode === "query-reference"
    ? referencePlan.tripletCount
    : exploratoryTripletCount;
  const activeCorrectionTestCount = options.analysisMode === "query-reference"
    ? referencePlan.correctionTestCount
    : Math.min(NATIVE_CORRECTION_CAP, exploratoryTripletCount);

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
        const restoredReferenceGroups = new Map<number, number>();
        (restored.results?.referenceGroupIndices ?? []).forEach((group, index) => {
          if (Number.isInteger(group) && group > 0 && index < restored.dataset.sequenceCount) {
            restoredReferenceGroups.set(index, group);
          }
        });
        setDataset(restored.dataset);
        setFilename(restored.sourceFilename);
        setFileSize(file.size);
        setMasked(restoredMask);
        setDisabled(restoredDisabled);
        setReferenceGroups(restoredReferenceGroups);
        setResults(restored.results);
        setCheckpointDirty(false);
        setOptions(
          restored.results
            ? {
                analysisMode: restored.results.analysisMode,
                circular: restored.results.circular,
                pValueCutoff: restored.results.pValueCutoff,
                correction: restored.results.correction,
                windowSites: restored.results.windowSites,
                maxChiEnabled: restored.results.maxChiEnabled,
                maxChiWindowSites: restored.results.maxChiWindowSites,
                chimaeraEnabled: restored.results.chimaeraEnabled,
                chimaeraWindowSites: restored.results.chimaeraWindowSites,
                geneconvEnabled: restored.results.geneconvEnabled ?? false,
                geneconvMismatchScale: restored.results.geneconvMismatchScale ?? 1,
                geneconvMaxOverlaps: restored.results.geneconvMaxOverlaps ?? 1,
                threeSeqEnabled: restored.results.threeSeqEnabled ?? false,
                bootscanPrimaryEnabled:
                  restored.results.bootscanPrimaryEnabled ?? false,
                bootscanSecondaryEnabled:
                  restored.results.bootscanSecondaryEnabled ?? false,
                bootscanWindowSites: restored.results.bootscanWindowSites ?? 200,
                bootscanStepSites: restored.results.bootscanStepSites ?? 20,
                bootscanBootstrapReplicates:
                  restored.results.bootscanBootstrapReplicates ?? 100,
                bootscanSupportCutoff:
                  restored.results.bootscanSupportCutoff ?? 0.7,
                bootscanRandomSeed: restored.results.bootscanRandomSeed ?? 3,
                polishBreakpoints: restored.results.polishBreakpoints ?? true,
                maskedSequenceIndices: [...restoredMask],
                disabledSequenceIndices: [...restoredDisabled],
                referenceGroupIndices: referenceGroupArray(
                  restored.dataset.sequenceCount,
                  restoredReferenceGroups,
                ),
              }
            : {
                ...initialOptions,
                maskedSequenceIndices: [...restoredMask],
                disabledSequenceIndices: [...restoredDisabled],
                referenceGroupIndices: referenceGroupArray(
                  restored.dataset.sequenceCount,
                  restoredReferenceGroups,
                ),
              },
        );
        setProgress(
          restored.results
            ? {
                state: "done",
                phase: "complete",
                processedTriplets: restored.results.processedTriplets,
                totalTriplets: restored.results.totalTriplets,
                correctionTests: restored.results.correctionTests,
                activeWorkingSequenceCount: restored.results.activeWorkingSequenceCount,
                queryWorkingSequenceCount: restored.results.queryWorkingSequenceCount,
                referenceWorkingSequenceCount: restored.results.referenceWorkingSequenceCount,
                activeReferenceGroupCount: restored.results.activeReferenceGroupCount,
                cumulativeTriplets: restored.results.cumulativeTriplets,
                tripletKernelEvaluations: 0,
                tripletSummariesReused: 0,
                cachedSignalsReused: 0,
                methodScansSkipped: 0,
                scanRound: restored.results.scanRounds,
                fixedEventCount: 0,
                signalCount: restored.results.signals.length,
                eventCount: restored.results.events.length,
                maxChiProfilesScanned: restored.results.maxChiProfilesScanned,
                maxChiPeakAttempts: restored.results.maxChiPeakAttempts,
                maxChiCandidatesFound: restored.results.maxChiCandidatesFound,
                maxChiPeakLimitTriplets: restored.results.maxChiPeakLimitTriplets,
                chimaeraProfilesScanned: restored.results.chimaeraProfilesScanned,
                chimaeraPeakAttempts: restored.results.chimaeraPeakAttempts,
                chimaeraCandidatesFound: restored.results.chimaeraCandidatesFound,
                chimaeraPeakLimitTargets: restored.results.chimaeraPeakLimitTargets,
                geneconvFragmentsScored: restored.results.geneconvFragmentsScored ?? 0,
                geneconvQualifiedFragments: restored.results.geneconvQualifiedFragments ?? 0,
                geneconvCandidatesFound: restored.results.geneconvCandidatesFound ?? 0,
                geneconvOverlapRejections: restored.results.geneconvOverlapRejections ?? 0,
                geneconvNumericalFallbackTracks:
                  restored.results.geneconvNumericalFallbackTracks ?? 0,
                threeSeqProfilesScanned: restored.results.threeSeqProfilesScanned ?? 0,
                threeSeqExactEvaluations: restored.results.threeSeqExactEvaluations ?? 0,
                threeSeqApproximateEvaluations:
                  restored.results.threeSeqApproximateEvaluations ?? 0,
                threeSeqCandidatesFound: restored.results.threeSeqCandidatesFound ?? 0,
                bootscanProfilesScanned: restored.results.bootscanProfilesScanned ?? 0,
                bootscanCandidateRegionsScored:
                  restored.results.bootscanCandidateRegionsScored ?? 0,
                bootscanCandidatesFound: restored.results.bootscanCandidatesFound ?? 0,
                bootscanPairProfilesRequested:
                  restored.results.bootscanPairProfilesRequested ?? 0,
                bootscanPairProfileCacheHits:
                  restored.results.bootscanPairProfileCacheHits ?? 0,
                bootscanPairProfileCacheMisses:
                  restored.results.bootscanPairProfileCacheMisses ?? 0,
                bootscanPairProfileCacheEvictions:
                  restored.results.bootscanPairProfileCacheEvictions ?? 0,
                bootscanPairProfileCacheBytes: 0,
                bootscanPairProfileCachePeakBytes:
                  restored.results.bootscanPairProfileCachePeakBytes ?? 0,
                cycleTermination: restored.results.cycleTermination,
                fraction: 1,
              }
            : {
                ...initialProgress,
                totalTriplets: chooseThree(scanEligibleSequenceCount(
                  restored.dataset,
                  restoredMask,
                  restoredDisabled,
                  Math.max(5, initialOptions.windowSites),
                )),
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
      const inferredReferenceGroups = inferReferenceGroups(summary);
      setMasked(initialMask);
      setDisabled(new Set());
      setReferenceGroups(inferredReferenceGroups);
      setOptions({
        ...initialOptions,
        maskedSequenceIndices: [...initialMask],
        disabledSequenceIndices: [],
        referenceGroupIndices: referenceGroupArray(
          summary.sequenceCount,
          inferredReferenceGroups,
        ),
      });
      setProgress({
        ...initialProgress,
        totalTriplets: chooseThree(scanEligibleSequenceCount(
          summary,
          initialMask,
          new Set<number>(),
          Math.max(5, initialOptions.windowSites),
        )),
      });
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
    const nextActiveCount = scanEligibleSequenceCount(
      dataset,
      nextMasked,
      nextDisabled,
      minimumWorkingSites,
    );
    const nextReferencePlan = queryReferencePlan(
      dataset,
      nextMasked,
      nextDisabled,
      referenceGroups,
      minimumWorkingSites,
    );
    setMasked(nextMasked);
    setDisabled(nextDisabled);
    setOptions((current) => ({
      ...current,
      maskedSequenceIndices: [...nextMasked].sort((left, right) => left - right),
      disabledSequenceIndices: [...nextDisabled].sort((left, right) => left - right),
    }));
    setProgress({
      ...initialProgress,
      totalTriplets: options.analysisMode === "query-reference"
        ? nextReferencePlan.tripletCount
        : chooseThree(nextActiveCount),
    });
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
    const nextActiveCount = scanEligibleSequenceCount(
      dataset,
      nextMasked,
      nextDisabled,
      minimumWorkingSites,
    );
    const nextReferencePlan = queryReferencePlan(
      dataset,
      nextMasked,
      nextDisabled,
      referenceGroups,
      minimumWorkingSites,
    );
    setMasked(nextMasked);
    setDisabled(nextDisabled);
    setOptions((current) => ({
      ...current,
      maskedSequenceIndices: [...nextMasked].sort((left, right) => left - right),
      disabledSequenceIndices: [...nextDisabled].sort((left, right) => left - right),
    }));
    setProgress({
      ...initialProgress,
      totalTriplets: options.analysisMode === "query-reference"
        ? nextReferencePlan.tripletCount
        : chooseThree(nextActiveCount),
    });
    setResults(null);
    setCheckpointDirty(false);
    setError("");
  };

  const changeOptions = (next: ScanOptions) => {
    if (!confirmDiscardUnsavedAnalysis()) return;
    const nextMinimumWorkingSites = Math.max(5, next.windowSites);
    const nextEligibleCount = scanEligibleSequenceCount(
      dataset,
      masked,
      disabled,
      nextMinimumWorkingSites,
    );
    const nextReferencePlan = queryReferencePlan(
      dataset,
      masked,
      disabled,
      referenceGroups,
      nextMinimumWorkingSites,
    );
    setOptions(next);
    setProgress({
      ...initialProgress,
      totalTriplets: next.analysisMode === "query-reference"
        ? nextReferencePlan.tripletCount
        : chooseThree(nextEligibleCount),
    });
    setResults(null);
    setCheckpointDirty(false);
    setError("");
  };

  const commitReferenceGroups = (nextGroups: Map<number, number>) => {
    if (!dataset) return;
    const indices = referenceGroupArray(dataset.sequenceCount, nextGroups);
    const nextReferencePlan = queryReferencePlan(
      dataset,
      masked,
      disabled,
      nextGroups,
      minimumWorkingSites,
    );
    setReferenceGroups(nextGroups);
    setOptions((current) => ({ ...current, referenceGroupIndices: indices }));
    setProgress({
      ...initialProgress,
      totalTriplets: options.analysisMode === "query-reference"
        ? nextReferencePlan.tripletCount
        : exploratoryTripletCount,
    });
    setResults(null);
    setCheckpointDirty(false);
    setError("");
  };

  const changeReferenceGroup = (index: number, group: number) => {
    if (!dataset || !confirmDiscardUnsavedAnalysis()) return;
    if (!Number.isInteger(index) || index < 0 || index >= dataset.sequenceCount) return;
    const nextGroups = new Map(referenceGroups);
    const normalized = normalizeReferenceGroup(group);
    if (normalized === 0) nextGroups.delete(index);
    else nextGroups.set(index, normalized);
    commitReferenceGroups(nextGroups);
  };

  const changeReferenceGroups = (indices: number[], group: number) => {
    if (!dataset || !confirmDiscardUnsavedAnalysis()) return;
    const nextGroups = new Map(referenceGroups);
    const normalized = normalizeReferenceGroup(group);
    new Set(indices).forEach((index) => {
      if (!Number.isInteger(index) || index < 0 || index >= dataset.sequenceCount) return;
      if (normalized === 0) nextGroups.delete(index);
      else nextGroups.set(index, normalized);
    });
    commitReferenceGroups(nextGroups);
  };

  const changeAllReferenceGroups = (action: "detect" | "compact" | "clear") => {
    if (!dataset || !confirmDiscardUnsavedAnalysis()) return;
    const nextGroups = action === "detect"
      ? inferReferenceGroups(dataset)
      : action === "compact"
        ? compactReferenceGroups(dataset, referenceGroups)
        : new Map<number, number>();
    commitReferenceGroups(nextGroups);
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
      referenceGroupIndices: referenceGroupArray(dataset.sequenceCount, referenceGroups),
    };
    setProgress({
      ...initialProgress,
      state: "running",
      totalTriplets: activeTripletCount,
      correctionTests: activeCorrectionTestCount,
      activeWorkingSequenceCount: activeSequenceCount,
      queryWorkingSequenceCount: options.analysisMode === "query-reference"
        ? referencePlan.querySequenceCount
        : 0,
      referenceWorkingSequenceCount: options.analysisMode === "query-reference"
        ? referencePlan.referenceSequenceCount
        : 0,
      activeReferenceGroupCount: options.analysisMode === "query-reference"
        ? referencePlan.referenceGroupCount
        : 0,
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
        <div className="window-titlebar">
          <button className="mobile-menu" type="button" onClick={() => setNavOpen((open) => !open)} aria-label="Toggle workflow navigation">
            {navOpen ? <PanelLeftClose /> : <Menu />}
          </button>
          <div className="brand">
            <span className="brand-mark"><Atom size={16} strokeWidth={2.4} /></span>
            <span>
              <strong>RDP Web - Recombination Detection Program</strong>
              <small>Browser-native analysis workstation</small>
            </span>
          </div>
          <div className="window-controls" aria-hidden="true">
            <span className="window-control window-minimize">_</span>
            <span className="window-control window-maximize">□</span>
            <span className="window-control window-close">×</span>
          </div>
        </div>
        <div className="window-menubar">
          <nav className="classic-menu" aria-label="Application menu">
            <span><u>F</u>ile</span>
            <span><u>D</u>ataset</span>
            <span><u>A</u>nalysis</span>
            <span><u>V</u>iew</span>
            <span><u>H</u>elp</span>
          </nav>
          <div className="topbar-status">
            <span className={`engine-pill engine-${engine.status}`} title={engine.message}>
              <Cpu size={13} />
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
            <span className="session-pill">Win95 edition · session 20</span>
          </div>
        </div>
      </header>

      <aside className={`sidebar${navOpen ? " is-open" : ""}`}>
        <div className="sidebar-caption">Analysis workflow</div>
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
            referenceGroups={referenceGroups}
            eligibleSequenceCount={activeSequenceCount}
            exploratoryTripletCount={exploratoryTripletCount}
            busy={loading}
            onLoad={loadAlignment}
            onSequenceStateChange={changeSequenceState}
            onAllSequenceStatesChange={changeAllSequenceStates}
            onReferenceGroupChange={changeReferenceGroup}
            onReferenceGroupsChange={changeReferenceGroups}
            onAllReferenceGroupsChange={changeAllReferenceGroups}
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
            exploratoryTripletCount={exploratoryTripletCount}
            queryReferenceTripletCount={referencePlan.tripletCount}
            querySequenceCount={referencePlan.querySequenceCount}
            referenceSequenceCount={referencePlan.referenceSequenceCount}
            referenceGroupCount={referencePlan.referenceGroupCount}
            queryReferenceCorrectionTestCount={referencePlan.correctionTestCount}
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
            correctionTestCount={activeCorrectionTestCount}
            querySequenceCount={referencePlan.querySequenceCount}
            referenceSequenceCount={referencePlan.referenceSequenceCount}
            referenceGroupCount={referencePlan.referenceGroupCount}
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

      <footer className="app-statusbar" aria-live="polite">
        <span>{engine.status === "ready" ? "Ready" : engine.status === "loading" ? "Loading analysis engine…" : "Engine unavailable"}</span>
        <span>{filename || "No alignment loaded"}</span>
        <span>{results ? `${results.events.length} event${results.events.length === 1 ? "" : "s"}` : "RDP Web 0.20"}</span>
      </footer>
    </div>
  );
}
