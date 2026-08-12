/// <reference lib="webworker" />

import type {
  DatasetSummary,
  ImportedProject,
  ScanProgress,
  ScanResults,
  WorkerRequest,
  WorkerResponse,
} from "../lib/types";

interface EmscriptenModule {
  HEAPU8: Uint8Array;
  _malloc(size: number): number;
  _free(pointer: number): void;
  UTF8ToString(pointer: number): string;
  _rdp_create(): number;
  _rdp_destroy(handle: number): void;
  _rdp_version(): number;
  _rdp_load_alignment(handle: number, bytes: number, length: number): number;
  _rdp_get_summary_json(handle: number): number;
  _rdp_scan_begin(
    handle: number,
    circular: number,
    correction: number,
    pValueCutoff: number,
    windowSites: number,
    mask: number,
    maskLength: number,
  ): number;
  _rdp_scan_batch(handle: number, tripletBudget: number): number;
  _rdp_reconcile(handle: number): number;
  _rdp_cancel(handle: number): void;
  _rdp_get_progress_json(handle: number): number;
  _rdp_get_results_json(handle: number): number;
  _rdp_get_signal_plot_json(handle: number, signalId: number): number;
  _rdp_set_review_state(handle: number, signalId: number, state: number): number;
  _rdp_set_event_review_state(handle: number, eventId: number, state: number): number;
  _rdp_update_event(
    handle: number,
    eventId: number,
    recombinant: number,
    majorParent: number,
    minorParent: number,
    beginning: number,
    ending: number,
  ): number;
  _rdp_update_event_group(
    handle: number,
    eventId: number,
    sequenceIndices: number,
    sequenceCount: number,
    manualOverride: number,
  ): number;
  _rdp_reconcile_after(handle: number, eventId: number): number;
  _rdp_restore_alignment_begin(handle: number, sequenceCount: number): number;
  _rdp_restore_alignment_record(
    handle: number,
    index: number,
    name: number,
    nameLength: number,
    sequence: number,
    sequenceLength: number,
  ): number;
  _rdp_restore_alignment_finish(handle: number, format: number, formatLength: number): number;
  _rdp_restore_scan_begin(
    handle: number,
    circular: number,
    correction: number,
    pValueCutoff: number,
    windowSites: number,
    mask: number,
    maskLength: number,
  ): number;
  _rdp_restore_signal(
    handle: number,
    triplet0: number,
    triplet1: number,
    triplet2: number,
    recombinant: number,
    majorParent: number,
    minorParent: number,
    beginning: number,
    ending: number,
    wrapsOrigin: number,
    informativeBeginning: number,
    informativeEnding: number,
    localPValue: number,
    correctedPValue: number,
    correctionTests: number,
    pairSimilarity0: number,
    pairSimilarity1: number,
    pairSimilarity2: number,
    informativeSites: number,
    candidatePair: number,
    fragmentAssisted: number,
    fragmentEvent0: number,
    fragmentEvent1: number,
    fragmentEvent2: number,
    reviewState: number,
    eventId: number,
  ): number;
  _rdp_restore_scan_finish(handle: number, correctionTests: number): number;
  _rdp_restore_event_state(
    handle: number,
    eventId: number,
    anchorSignalId: number,
    recombinant: number,
    majorParent: number,
    minorParent: number,
    beginning: number,
    ending: number,
    detectionRound: number,
    tractErasedForDetection: number,
    reviewState: number,
    manualAdjusted: number,
    coRecombinantSequences: number,
    coRecombinantSequenceCount: number,
    groupManualAdjusted: number,
  ): number;
  _rdp_restore_reconciliation_required_after(handle: number, eventId: number): number;
  _rdp_export_csv(handle: number): number;
  _rdp_export_recombination_free_fasta(handle: number): number;
  _rdp_export_fragmented_fasta(handle: number): number;
  _rdp_export_project_json(handle: number): number;
  _rdp_get_error(handle: number): number;
}

type ModuleFactory = (options: {
  locateFile: (path: string) => string;
  noInitialRun: boolean;
}) => Promise<EmscriptenModule>;

let module: EmscriptenModule | null = null;
let context = 0;
let dataset: DatasetSummary | null = null;
let datasetName = "";
let threaded = false;
let scanActive = false;

const scope = self as DedicatedWorkerGlobalScope;

function respond(response: WorkerResponse): void {
  scope.postMessage(response);
}

function value(pointer: number): string {
  if (!module || pointer === 0) return "";
  return module.UTF8ToString(pointer);
}

function engineError(fallback: string): Error {
  if (!module || !context) return new Error(fallback);
  return new Error(value(module._rdp_get_error(context)) || fallback);
}

function parseJson<T>(pointer: number, fallback: string): T {
  const text = value(pointer);
  if (!text) throw engineError(fallback);
  return JSON.parse(text) as T;
}

async function importFactory(url: string): Promise<ModuleFactory> {
  const imported = (await import(/* @vite-ignore */ url)) as { default: ModuleFactory };
  return imported.default;
}

async function initialise(wasmBaseUrl: string): Promise<{ threaded: boolean; version: string }> {
  if (module) return { threaded, version: value(module._rdp_version()) };

  const base = wasmBaseUrl.endsWith("/") ? wasmBaseUrl : `${wasmBaseUrl}/`;
  const canThread = scope.crossOriginIsolated && typeof SharedArrayBuffer !== "undefined";
  const candidates = canThread
    ? [
        { name: "rdp-core-threads.mjs", threaded: true },
        { name: "rdp-core.mjs", threaded: false },
      ]
    : [{ name: "rdp-core.mjs", threaded: false }];

  let lastError: unknown = null;
  for (const candidate of candidates) {
    try {
      const factory = await importFactory(`${base}${candidate.name}`);
      module = await factory({
        noInitialRun: true,
        locateFile: (path) => `${base}${path}`,
      });
      context = module._rdp_create();
      if (!context) throw new Error("The WASM engine could not allocate an analysis context.");
      threaded = candidate.threaded;
      return { threaded, version: value(module._rdp_version()) };
    } catch (error) {
      module = null;
      context = 0;
      lastError = error;
    }
  }

  const detail = lastError instanceof Error ? lastError.message : String(lastError ?? "unknown error");
  throw new Error(
    `The RDP WASM module is not available (${detail}). Build it with npm run build:wasm before running the app.`,
  );
}

function copyBytes(bytes: Uint8Array): number {
  if (!module) throw new Error("The engine has not been initialised.");
  const pointer = module._malloc(Math.max(1, bytes.byteLength));
  if (!pointer) throw new Error("The WASM heap could not allocate the input alignment.");
  module.HEAPU8.set(bytes, pointer);
  return pointer;
}

function copyUint32(values: number[]): number {
  const encoded = new Uint32Array(values);
  return copyBytes(new Uint8Array(encoded.buffer));
}

function reviewStateCode(state: unknown): number {
  return state === "accepted" ? 1 : state === "rejected" ? 2 : 0;
}

function finiteNumber(value: unknown, fallback = 0): number {
  return typeof value === "number" && Number.isFinite(value) ? value : fallback;
}

function integer(value: unknown, fallback = 0): number {
  return Math.trunc(finiteNumber(value, fallback));
}

function requireObject(value: unknown, message: string): Record<string, unknown> {
  if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error(message);
  return value as Record<string, unknown>;
}

function requireString(value: unknown, message: string): string {
  if (typeof value !== "string") throw new Error(message);
  return value;
}

function requireArray(value: unknown, message: string): unknown[] {
  if (!Array.isArray(value)) throw new Error(message);
  return value;
}

function loadAlignment(name: string, bytes: ArrayBuffer): DatasetSummary {
  if (!module || !context) throw new Error("The engine has not been initialised.");
  const input = new Uint8Array(bytes);
  const pointer = copyBytes(input);
  try {
    if (module._rdp_load_alignment(context, pointer, input.byteLength) !== 1) {
      throw engineError("The alignment could not be loaded.");
    }
    dataset = parseJson<DatasetSummary>(
      module._rdp_get_summary_json(context),
      "The alignment summary was not returned.",
    );
    datasetName = name;
    return dataset;
  } finally {
    module._free(pointer);
  }
}

function restoreProject(name: string, bytes: ArrayBuffer): ImportedProject {
  if (!module || !context) throw new Error("The engine has not been initialised.");
  let root: Record<string, unknown>;
  try {
    root = requireObject(
      JSON.parse(new TextDecoder().decode(bytes)),
      "The selected file is not an RDP Web project.",
    );
  } catch (error) {
    if (error instanceof SyntaxError) throw new Error("The selected project contains invalid JSON.");
    throw error;
  }
  const schema = requireString(root.schema, "The project schema identifier is missing.");
  if (
    schema !== "org.rdp-web.project/v1alpha1" &&
    schema !== "org.rdp-web.project/v1alpha2" &&
    schema !== "org.rdp-web.project/v1alpha3" &&
    schema !== "org.rdp-web.project/v1alpha4" &&
    schema !== "org.rdp-web.project/v1alpha5" &&
    schema !== "org.rdp-web.project/v1alpha6"
  ) {
    throw new Error(`Unsupported RDP Web project schema: ${schema}`);
  }
  const savedDataset = requireObject(root.dataset, "The project has no saved alignment.");
  const records = requireArray(savedDataset.sequences, "The project has no saved sequences.");
  if (records.length < 3) throw new Error("The saved project contains fewer than three sequences.");

  if (module._rdp_restore_alignment_begin(context, records.length) !== 1) {
    throw engineError("The saved alignment could not be prepared for restoration.");
  }
  const encoder = new TextEncoder();
  records.forEach((record, index) => {
    const value = requireObject(record, `Saved sequence ${index + 1} is invalid.`);
    const sequenceName = requireString(value.name, `Saved sequence ${index + 1} has no name.`);
    const sequence = requireString(value.sequence, `Saved sequence ${index + 1} has no data.`);
    const nameBytes = encoder.encode(sequenceName);
    const sequenceBytes = encoder.encode(sequence);
    const namePointer = copyBytes(nameBytes);
    const sequencePointer = copyBytes(sequenceBytes);
    try {
      if (
        module!._rdp_restore_alignment_record(
          context,
          index,
          namePointer,
          nameBytes.byteLength,
          sequencePointer,
          sequenceBytes.byteLength,
        ) !== 1
      ) {
        throw engineError(`Saved sequence ${index + 1} could not be restored.`);
      }
    } finally {
      module!._free(namePointer);
      module!._free(sequencePointer);
    }
  });
  const formatBytes = encoder.encode(
    typeof savedDataset.format === "string" ? savedDataset.format : "RDP Web project",
  );
  const formatPointer = copyBytes(formatBytes);
  try {
    if (
      module._rdp_restore_alignment_finish(context, formatPointer, formatBytes.byteLength) !== 1
    ) {
      throw engineError("The saved alignment could not be restored.");
    }
  } finally {
    module._free(formatPointer);
  }

  const restoredDataset = parseJson<DatasetSummary>(
    module._rdp_get_summary_json(context),
    "The restored alignment summary was not returned.",
  );
  dataset = restoredDataset;
  datasetName = typeof root.sourceFilename === "string" ? root.sourceFilename : name;

  const rawAnalysis = root.analysis;
  if (!rawAnalysis) return { dataset: restoredDataset, results: null, sourceFilename: datasetName };
  const analysis = requireObject(rawAnalysis, "The saved analysis is invalid.");
  const savedSignals = requireArray(analysis.signals, "The saved analysis has no primary signals.");
  const savedEvents = Array.isArray(analysis.events) ? analysis.events : [];
  const pending = integer(analysis.downstreamReconciliationRequiredAfter, -1);
  if (pending >= savedEvents.length) {
    throw new Error("The saved downstream reconciliation marker refers to an unknown event.");
  }
  const eventsToRestore = (pending >= 0 ? savedEvents.slice(0, pending + 1) : savedEvents).map(
    (rawEvent, index) =>
      requireObject(rawEvent, `Saved event ${index + 1} is invalid.`),
  );

  // A role, breakpoint, group, or rejection correction invalidates every event
  // detected after the changed event. Saved projects retain that stale tail for
  // auditability, but a reload must rebuild only the valid prefix and its signal
  // evidence before the user resumes cyclic screening.
  const retainedSignalIds = new Set<number>();
  if (pending >= 0) {
    eventsToRestore.forEach((savedEvent) => {
      const anchorSignalId = integer(savedEvent.anchorSignalId, -1);
      if (anchorSignalId >= 0) retainedSignalIds.add(anchorSignalId);
      if (Array.isArray(savedEvent.supportSignalIds)) {
        savedEvent.supportSignalIds.forEach((signalId) => {
          const index = integer(signalId, -1);
          if (index >= 0) retainedSignalIds.add(index);
        });
      }
    });
    savedSignals.forEach((rawSignal, index) => {
      const signal = requireObject(rawSignal, `Saved signal ${index + 1} is invalid.`);
      const eventId = integer(signal.eventId, -1);
      if (eventId >= 0 && eventId <= pending) retainedSignalIds.add(index);
    });
  }

  const signalsToRestore = savedSignals
    .map((rawSignal, savedIndex) => ({ rawSignal, savedIndex }))
    .filter(({ savedIndex }) => pending < 0 || retainedSignalIds.has(savedIndex));
  const signalIdMap = new Map<number, number>();
  signalsToRestore.forEach(({ savedIndex }, restoredIndex) => {
    signalIdMap.set(savedIndex, restoredIndex);
  });
  const masked = new Uint8Array(restoredDataset.sequenceCount);
  if (Array.isArray(analysis.maskedSequenceIndices)) {
    analysis.maskedSequenceIndices.forEach((value) => {
      const index = integer(value, -1);
      if (index >= 0 && index < masked.length) masked[index] = 1;
    });
  }
  const maskPointer = copyBytes(masked);
  try {
    if (
      module._rdp_restore_scan_begin(
        context,
        analysis.circular === false ? 0 : 1,
        analysis.correction === "none" ? 1 : 0,
        finiteNumber(analysis.pValueCutoff, 0.05),
        integer(analysis.windowSites, 30),
        maskPointer,
        masked.length,
      ) !== 1
    ) {
      throw engineError("The saved scan settings could not be restored.");
    }
  } finally {
    module._free(maskPointer);
  }

  const inferredEventIds = new Map<number, number>();
  eventsToRestore.forEach((savedEvent, eventIndex) => {
    if (!Array.isArray(savedEvent.supportSignalIds)) return;
    savedEvent.supportSignalIds.forEach((signalId) => {
      const index = integer(signalId, -1);
      if (index >= 0) inferredEventIds.set(index, eventIndex);
    });
  });

  signalsToRestore.forEach(({ rawSignal, savedIndex }) => {
    const signal = requireObject(rawSignal, `Saved signal ${savedIndex + 1} is invalid.`);
    const triplet = requireArray(
      signal.triplet,
      `Saved signal ${savedIndex + 1} has no triplet.`,
    );
    const similarity = Array.isArray(signal.pairSimilarity) ? signal.pairSimilarity : [];
    const fragmentContext = Array.isArray(signal.fragmentEventContext)
      ? signal.fragmentEventContext
      : [];
    const explicitEventId = integer(signal.eventId, -1);
    const restoredEventId = inferredEventIds.get(savedIndex) ??
      (explicitEventId >= 0 && explicitEventId < eventsToRestore.length ? explicitEventId : -1);
    if (
      module!._rdp_restore_signal(
        context,
        integer(triplet[0], -1),
        integer(triplet[1], -1),
        integer(triplet[2], -1),
        integer(signal.recombinant, -1),
        integer(signal.majorParent, -1),
        integer(signal.minorParent, -1),
        integer(signal.beginning, -1),
        integer(signal.ending, -1),
        signal.wrapsOrigin === true ? 1 : 0,
        integer(signal.informativeBeginning),
        integer(signal.informativeEnding),
        finiteNumber(signal.localPValue, 1),
        finiteNumber(signal.correctedPValue, 1),
        integer(signal.correctionTests, integer(analysis.correctionTests)),
        finiteNumber(similarity[0]),
        finiteNumber(similarity[1]),
        finiteNumber(similarity[2]),
        integer(signal.informativeSites),
        integer(signal.candidatePair),
        signal.fragmentAssisted === true ? 1 : 0,
        integer(fragmentContext[0], -1),
        integer(fragmentContext[1], -1),
        integer(fragmentContext[2], -1),
        reviewStateCode(signal.reviewState),
        restoredEventId,
      ) !== 1
    ) {
      throw engineError(`Saved signal ${savedIndex + 1} could not be restored.`);
    }
  });
  if (
    module._rdp_restore_scan_finish(context, integer(analysis.correctionTests)) !== 1
  ) {
    throw engineError("The saved primary analysis could not be restored.");
  }

  eventsToRestore.forEach((savedEvent, index) => {
    const savedAnchorSignalId = integer(savedEvent.anchorSignalId, -1);
    const restoredAnchorSignalId = signalIdMap.get(savedAnchorSignalId);
    if (restoredAnchorSignalId === undefined) {
      throw new Error(`Saved event ${index + 1} has no restorable anchor signal.`);
    }
    const group = Array.isArray(savedEvent.coRecombinantSequenceIndices)
      ? savedEvent.coRecombinantSequenceIndices.map((value) => integer(value, -1))
      : [];
    const groupPointer = copyUint32(group);
    try {
      if (
        module!._rdp_restore_event_state(
          context,
          index,
          restoredAnchorSignalId,
          integer(savedEvent.recombinant, -1),
          integer(savedEvent.majorParent, -1),
          integer(savedEvent.minorParent, -1),
          integer(savedEvent.beginning, -1),
          integer(savedEvent.ending, -1),
          integer(savedEvent.detectionRound, index + 1),
          savedEvent.tractErasedForDetection === false ? 0 : 1,
          reviewStateCode(savedEvent.reviewState),
          savedEvent.manualAdjusted === true ? 1 : 0,
          groupPointer,
          group.length,
          savedEvent.groupManualAdjusted === true ? 1 : 0,
        ) !== 1
      ) {
        throw engineError(`Saved event ${index + 1} could not be restored.`);
      }
    } finally {
      module!._free(groupPointer);
    }
  });
  if (
    pending >= 0 &&
    module._rdp_restore_reconciliation_required_after(context, pending) !== 1
  ) {
    throw engineError("The saved downstream reconciliation marker could not be restored.");
  }
  const results = parseJson<ScanResults>(
    module._rdp_get_results_json(context),
    "The restored analysis was not returned.",
  );
  return { dataset: restoredDataset, results, sourceFilename: datasetName };
}

function exportProject(): string {
  if (!module || !context) throw new Error("The engine has not been initialised.");
  const raw = value(module._rdp_export_project_json(context));
  if (!raw) throw engineError("The project snapshot could not be exported.");
  const project = JSON.parse(raw) as Record<string, unknown>;
  project.sourceFilename = datasetName;
  return JSON.stringify(project, null, 2);
}

function emitProgress(): ScanProgress {
  if (!module || !context) throw new Error("The engine has not been initialised.");
  const progress = parseJson<ScanProgress>(
    module._rdp_get_progress_json(context),
    "Scan progress was not returned.",
  );
  scope.postMessage({ type: "progress", progress });
  return progress;
}

async function runScan(request: Extract<WorkerRequest, { type: "scan" }>): Promise<unknown> {
  if (!module || !context || !dataset) throw new Error("Load an alignment before starting a scan.");
  if (scanActive) throw new Error("A scan is already running.");

  const mask = new Uint8Array(dataset.sequenceCount);
  request.options.maskedSequenceIndices.forEach((index) => {
    if (index >= 0 && index < mask.length) mask[index] = 1;
  });
  const maskPointer = copyBytes(mask);
  try {
    const correction = request.options.correction === "bonferroni" ? 0 : 1;
    const started = module._rdp_scan_begin(
      context,
      request.options.circular ? 1 : 0,
      correction,
      request.options.pValueCutoff,
      request.options.windowSites,
      maskPointer,
      mask.length,
    );
    if (started !== 1) throw engineError("The RDP scan could not be started.");
  } finally {
    module._free(maskPointer);
  }

  scanActive = true;
  emitProgress();
  try {
    for (;;) {
      const status = module._rdp_scan_batch(context, 512);
      emitProgress();
      if (status === 1) break;
      if (status === 3) {
        await new Promise<void>((resolve) => setTimeout(resolve, 0));
        if (module._rdp_reconcile(context) !== 1) {
          throw engineError("The primary signals could not be reconciled into event hypotheses.");
        }
        emitProgress();
        break;
      }
      if (status === 2) throw new Error("The scan was cancelled.");
      if (status < 0) throw engineError("The scan failed.");
      await new Promise<void>((resolve) => setTimeout(resolve, 0));
    }
    return parseJson(module._rdp_get_results_json(context), "The scan returned no results.");
  } finally {
    scanActive = false;
  }
}

async function reidentifyLaterEvents(eventId: number): Promise<ScanResults> {
  if (!module || !context) throw new Error("The engine has not been initialised.");
  if (scanActive) throw new Error("A scan is already running.");
  if (module._rdp_reconcile_after(context, eventId) !== 1) {
    throw engineError("Later events could not be prepared for re-identification.");
  }
  scanActive = true;
  emitProgress();
  try {
    for (;;) {
      const status = module._rdp_scan_batch(context, 512);
      emitProgress();
      if (status === 3) {
        await new Promise<void>((resolve) => setTimeout(resolve, 0));
        if (module._rdp_reconcile(context) !== 1) {
          throw engineError("The rebuilt signals could not be reconciled into event hypotheses.");
        }
        emitProgress();
        break;
      }
      if (status === 2) throw new Error("Re-identification was cancelled.");
      if (status < 0) throw engineError("Re-identification failed.");
      await new Promise<void>((resolve) => setTimeout(resolve, 0));
    }
    return parseJson<ScanResults>(
      module._rdp_get_results_json(context),
      "Re-identified event results were not returned.",
    );
  } finally {
    scanActive = false;
  }
}

scope.addEventListener("message", async (event: MessageEvent<WorkerRequest>) => {
  const request = event.data;
  try {
    let result: unknown;
    switch (request.type) {
      case "init":
        result = await initialise(request.wasmBaseUrl);
        break;
      case "load":
        result = loadAlignment(request.name, request.bytes);
        break;
      case "import-project":
        result = restoreProject(request.name, request.bytes);
        break;
      case "scan":
        result = await runScan(request);
        break;
      case "cancel":
        if (module && context) module._rdp_cancel(context);
        result = undefined;
        break;
      case "plot":
        if (!module || !context) throw new Error("The engine has not been initialised.");
        result = parseJson(
          module._rdp_get_signal_plot_json(context, request.signalId),
          "Plot data was not returned.",
        );
        break;
      case "set-review-state": {
        if (!module || !context) throw new Error("The engine has not been initialised.");
        const state = request.state === "accepted" ? 1 : request.state === "rejected" ? 2 : 0;
        if (module._rdp_set_review_state(context, request.signalId, state) !== 1) {
          throw engineError("The review state could not be changed.");
        }
        result = undefined;
        break;
      }
      case "set-event-review-state": {
        if (!module || !context) throw new Error("The engine has not been initialised.");
        if (
          module._rdp_set_event_review_state(
            context,
            request.eventId,
            reviewStateCode(request.state),
          ) !== 1
        ) {
          throw engineError("The event review state could not be changed.");
        }
        result = parseJson(
          module._rdp_get_results_json(context),
          "Updated event results were not returned.",
        );
        break;
      }
      case "update-event":
        if (!module || !context) throw new Error("The engine has not been initialised.");
        if (
          module._rdp_update_event(
            context,
            request.eventId,
            request.edit.recombinant,
            request.edit.majorParent,
            request.edit.minorParent,
            request.edit.beginning,
            request.edit.ending,
          ) !== 1
        ) {
          throw engineError("The event correction could not be saved.");
        }
        result = parseJson(
          module._rdp_get_results_json(context),
          "Updated event results were not returned.",
        );
        break;
      case "update-event-group": {
        if (!module || !context) throw new Error("The engine has not been initialised.");
        const sequenceIndices = [...new Set(request.sequenceIndices.map((value) => integer(value, -1)))]
          .sort((left, right) => left - right);
        const pointer = copyUint32(sequenceIndices);
        try {
          if (
            module._rdp_update_event_group(
              context,
              request.eventId,
              pointer,
              sequenceIndices.length,
              request.manualOverride ? 1 : 0,
            ) !== 1
          ) {
            throw engineError("The co-recombinant group correction could not be saved.");
          }
        } finally {
          module._free(pointer);
        }
        result = parseJson(
          module._rdp_get_results_json(context),
          "Updated event results were not returned.",
        );
        break;
      }
      case "reconcile-after":
        result = await reidentifyLaterEvents(request.eventId);
        break;
      case "export-csv":
        if (!module || !context) throw new Error("The engine has not been initialised.");
        result = value(module._rdp_export_csv(context));
        break;
      case "export-recombination-free": {
        if (!module || !context) throw new Error("The engine has not been initialised.");
        const fasta = value(module._rdp_export_recombination_free_fasta(context));
        if (!fasta) throw engineError("The final tract-masked alignment is not ready.");
        result = fasta;
        break;
      }
      case "export-fragmented": {
        if (!module || !context) throw new Error("The engine has not been initialised.");
        const fasta = value(module._rdp_export_fragmented_fasta(context));
        if (!fasta) throw engineError("The final mosaic-fragment alignment is not ready.");
        result = fasta;
        break;
      }
      case "export-project":
        result = exportProject();
        break;
      default:
        throw new Error("Unknown worker request.");
    }
    respond({ id: request.id, ok: true, value: result });
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    respond({ id: request.id, ok: false, error: message });
  }
});
