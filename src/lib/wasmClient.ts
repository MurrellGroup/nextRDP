import type {
  DatasetSummary,
  EventEdit,
  ImportedProject,
  ReviewState,
  ScanOptions,
  ScanProgress,
  ScanResults,
  SignalPlot,
  WorkerEvent,
  WorkerRequest,
  WorkerRequestPayload,
  WorkerResponse,
} from "./types";

type Pending = {
  resolve: (value: unknown) => void;
  reject: (reason: Error) => void;
};

export class RdpWorkerClient {
  private readonly worker: Worker;
  private readonly pending = new Map<number, Pending>();
  private requestId = 0;
  private progressListeners = new Set<(progress: ScanProgress) => void>();

  constructor() {
    this.worker = new Worker(new URL("../workers/analysis.worker.ts", import.meta.url), {
      type: "module",
      name: "rdp-analysis",
    });
    this.worker.addEventListener("message", (event: MessageEvent<WorkerResponse | WorkerEvent>) => {
      const message = event.data;
      if ("type" in message) {
        if (message.type === "progress") {
          this.progressListeners.forEach((listener) => listener(message.progress));
        }
        return;
      }
      const pending = this.pending.get(message.id);
      if (!pending) return;
      this.pending.delete(message.id);
      if (message.ok) pending.resolve(message.value);
      else pending.reject(new Error(message.error));
    });
    this.worker.addEventListener("error", (event) => {
      const error = new Error(event.message || "The RDP analysis worker stopped unexpectedly.");
      this.pending.forEach(({ reject }) => reject(error));
      this.pending.clear();
    });
  }

  onProgress(listener: (progress: ScanProgress) => void): () => void {
    this.progressListeners.add(listener);
    return () => this.progressListeners.delete(listener);
  }

  init(): Promise<{ threaded: boolean; version: string }> {
    const wasmBaseUrl = new URL("wasm/", document.baseURI).href;
    return this.send({ type: "init", wasmBaseUrl }) as Promise<{
      threaded: boolean;
      version: string;
    }>;
  }

  async load(file: File): Promise<DatasetSummary> {
    const bytes = await file.arrayBuffer();
    return this.send({ type: "load", name: file.name, bytes }, [bytes]) as Promise<DatasetSummary>;
  }

  async importProject(file: File): Promise<ImportedProject> {
    const bytes = await file.arrayBuffer();
    return this.send(
      { type: "import-project", name: file.name, bytes },
      [bytes],
    ) as Promise<ImportedProject>;
  }

  scan(options: ScanOptions): Promise<ScanResults> {
    return this.send({ type: "scan", options }) as Promise<ScanResults>;
  }

  cancel(): Promise<void> {
    return this.send({ type: "cancel" }) as Promise<void>;
  }

  plot(signalId: number): Promise<SignalPlot> {
    return this.send({ type: "plot", signalId }) as Promise<SignalPlot>;
  }

  setReviewState(signalId: number, state: ReviewState): Promise<void> {
    return this.send({ type: "set-review-state", signalId, state }) as Promise<void>;
  }

  setEventReviewState(eventId: number, state: ReviewState): Promise<ScanResults> {
    return this.send({ type: "set-event-review-state", eventId, state }) as Promise<ScanResults>;
  }

  updateEvent(eventId: number, edit: EventEdit): Promise<ScanResults> {
    return this.send({ type: "update-event", eventId, edit }) as Promise<ScanResults>;
  }

  updateEventGroup(
    eventId: number,
    sequenceIndices: number[],
    manualOverride = true,
  ): Promise<ScanResults> {
    return this.send({ type: "update-event-group", eventId, sequenceIndices, manualOverride }) as Promise<ScanResults>;
  }

  reconcileAfter(eventId: number): Promise<ScanResults> {
    return this.send({ type: "reconcile-after", eventId }) as Promise<ScanResults>;
  }

  exportCsv(): Promise<string> {
    return this.send({ type: "export-csv" }) as Promise<string>;
  }

  exportRecombinationFree(): Promise<string> {
    return this.send({ type: "export-recombination-free" }) as Promise<string>;
  }

  exportFragmented(): Promise<string> {
    return this.send({ type: "export-fragmented" }) as Promise<string>;
  }

  exportProject(): Promise<string> {
    return this.send({ type: "export-project" }) as Promise<string>;
  }

  dispose(): void {
    this.worker.terminate();
    this.pending.clear();
    this.progressListeners.clear();
  }

  private send(
    request: WorkerRequestPayload,
    transfer: Transferable[] = [],
  ): Promise<unknown> {
    const id = ++this.requestId;
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      this.worker.postMessage({ ...request, id } as WorkerRequest, transfer);
    });
  }
}
