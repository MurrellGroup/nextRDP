import { lstat, readFile, readdir } from "node:fs/promises";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { createHash } from "node:crypto";

const output = resolve(process.cwd(), "dist");
const requiredFiles = [
  "index.html",
  ".nojekyll",
  "wasm/rdp-core.mjs",
  "wasm/rdp-core.wasm",
];

function fail(message) {
  throw new Error(`GitHub Pages artifact check failed: ${message}`);
}

async function requireFile(relativePath) {
  const path = resolve(output, relativePath);
  const metadata = await lstat(path).catch(() => null);
  if (!metadata?.isFile() || (relativePath !== ".nojekyll" && metadata.size === 0)) {
    fail(`${relativePath} is missing or empty`);
  }
  return path;
}

let totalBytes = 0;
async function inspectTree(directory) {
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    const path = resolve(directory, entry.name);
    const metadata = await lstat(path);
    if (metadata.isSymbolicLink()) fail(`${path} is a symbolic link`);
    if (metadata.isDirectory()) {
      await inspectTree(path);
      continue;
    }
    if (!metadata.isFile()) fail(`${path} is not a regular file`);
    if (metadata.nlink > 1) fail(`${path} is hard-linked`);
    totalBytes += metadata.size;
  }
}

const [indexPath, , modulePath, wasmPath] = await Promise.all(
  requiredFiles.map(requireFile),
);

const html = await readFile(indexPath, "utf8");
if (html.includes("/src/main.tsx")) fail("index.html still references development source");
for (const match of html.matchAll(/\b(?:href|src)=["']([^"']+)["']/g)) {
  const url = match[1];
  if (url.startsWith("/") && !url.startsWith("//")) {
    fail(`index.html contains root-relative asset URL ${url}`);
  }
}

const moduleSource = await readFile(modulePath, "utf8");
if (!moduleSource.includes("rdp-core.wasm")) {
  fail("the Emscripten module does not reference rdp-core.wasm");
}

const wasm = await readFile(wasmPath);
if (wasm.length < 8 || !wasm.subarray(0, 4).equals(Buffer.from([0x00, 0x61, 0x73, 0x6d]))) {
  fail("rdp-core.wasm does not have a WebAssembly header");
}

const { default: createRdpModule } = await import(
  `${pathToFileURL(modulePath).href}?verify=${Date.now()}`,
);
const engine = await createRdpModule({
  noInitialRun: true,
  wasmBinary: wasm,
});
if (!(engine.HEAPU8 instanceof Uint8Array)) {
  fail("the Emscripten module does not expose its HEAPU8 upload view");
}

const context = engine._rdp_create();
if (!context) fail("the WASM engine could not create a FASTA smoke-test context");
const fasta = new TextEncoder().encode(
  ">alpha\nACGTACGT\n>beta\nACGTACGA\n>gamma\nACGTTCGA\n",
);
const fastaPointer = engine._malloc(fasta.byteLength);
if (!fastaPointer) {
  engine._rdp_destroy(context);
  fail("the WASM engine could not allocate the FASTA smoke-test input");
}
try {
  engine.HEAPU8.set(fasta, fastaPointer);
  if (engine._rdp_load_alignment(context, fastaPointer, fasta.byteLength) !== 1) {
    const detail = engine.UTF8ToString(engine._rdp_get_error(context));
    fail(`FASTA upload smoke test failed${detail ? `: ${detail}` : ""}`);
  }
  const summaryPointer = engine._rdp_get_summary_json(context);
  const summary = JSON.parse(engine.UTF8ToString(summaryPointer));
  if (summary.sequenceCount !== 3 || summary.alignmentLength !== 8) {
    fail("FASTA upload smoke test returned an unexpected alignment summary");
  }
} finally {
  engine._free(fastaPointer);
  engine._rdp_destroy(context);
}

// Exercise two cyclic detections, not just parsing. The synthetic alignment
// has two independent mosaics plus four background records. A later round
// gains fragment rows, so this also proves that MakeMCCorrection remains tied
// to the initial 10-record plan while unchanged triplets use the native-style
// XOverList/BestXOList shortlist instead of re-running their RDP kernel.
let syntheticSeed = 123456789;
const syntheticRandom = () => {
  syntheticSeed = (1664525 * syntheticSeed + 1013904223) >>> 0;
  return syntheticSeed;
};
const bases = "ACGT";
const randomSequence = (length) => Array.from(
  { length },
  () => bases[syntheticRandom() % bases.length],
).join("");
const periodicMutant = (sequence, period) => Array.from(
  sequence,
  (base, index) => index % period === 0
    ? bases[(bases.indexOf(base) + 1 + (index % 2)) % bases.length]
    : base,
).join("");
const mosaic = (first, second, split) => Array.from(first, (base, index) => {
  let result = index < split ? base : second[index];
  if (index % 47 === 0 && first[index] === second[index]) {
    result = bases[(bases.indexOf(result) + 2) % bases.length];
  }
  return result;
}).join("");

const syntheticLength = 900;
const parentA = randomSequence(syntheticLength);
const parentB = periodicMutant(parentA, 3);
const recombinantA = mosaic(parentA, parentB, 450);
const parentC = randomSequence(syntheticLength);
const parentD = periodicMutant(parentC, 4);
const recombinantB = mosaic(parentC, parentD, 300);
const syntheticSequences = [
  parentA,
  parentB,
  recombinantA,
  parentC,
  parentD,
  recombinantB,
  ...Array.from({ length: 4 }, () => randomSequence(syntheticLength)),
];
const syntheticFasta = new TextEncoder().encode(
  `${syntheticSequences.map((sequence, index) =>
    `>cyclic-${index}\n${sequence}`).join("\n")}\n`,
);
const cyclicContext = engine._rdp_create();
const syntheticPointer = engine._malloc(syntheticFasta.byteLength);
const referenceGroupsPointer = engine._malloc(syntheticSequences.length * 4);
const zeroFlagsPointer = engine._malloc(syntheticSequences.length);
if (!cyclicContext || !syntheticPointer || !referenceGroupsPointer || !zeroFlagsPointer) {
  fail("the WASM engine could not allocate the cyclic-shortlist smoke test");
}
try {
  engine.HEAPU8.set(syntheticFasta, syntheticPointer);
  engine.HEAPU8.fill(
    0,
    referenceGroupsPointer,
    referenceGroupsPointer + syntheticSequences.length * 4,
  );
  engine.HEAPU8.fill(
    0,
    zeroFlagsPointer,
    zeroFlagsPointer + syntheticSequences.length,
  );
  if (engine._rdp_load_alignment(
    cyclicContext,
    syntheticPointer,
    syntheticFasta.byteLength,
  ) !== 1) {
    fail(`cyclic-shortlist alignment load failed: ${engine.UTF8ToString(
      engine._rdp_get_error(cyclicContext),
    )}`);
  }
  const scanStarted = engine._rdp_scan_begin(
    cyclicContext,
    1, // circular
    1, // no correction, while retaining the source opportunity count
    0.05,
    30,
    0, 70, // MaxChi
    0, 60, // CHIMAERA
    0, 1, 1, // GENECONV
    0, // 3SEQ
    0, 200, 20, 100, 0.7, 3, // secondary BootScan
    0, // preserve detected breakpoints
    0, // exploratory mode
    referenceGroupsPointer,
    syntheticSequences.length,
    zeroFlagsPointer,
    syntheticSequences.length,
    zeroFlagsPointer,
    syntheticSequences.length,
  );
  if (scanStarted !== 1) {
    fail(`cyclic-shortlist scan could not start: ${engine.UTF8ToString(
      engine._rdp_get_error(cyclicContext),
    )}`);
  }
  let status = 0;
  for (let batch = 0; batch < 10000 && status !== 3; ++batch) {
    status = engine._rdp_scan_batch(cyclicContext, 10000);
    if (status < 0) {
      fail(`cyclic-shortlist scan failed: ${engine.UTF8ToString(
        engine._rdp_get_error(cyclicContext),
      )}`);
    }
  }
  const progress = JSON.parse(engine.UTF8ToString(
    engine._rdp_get_progress_json(cyclicContext),
  ));
  if (status !== 3 || progress.eventCount < 2 || progress.scanRound < 3) {
    fail("cyclic-shortlist smoke test did not complete two detection rounds");
  }
  if (progress.correctionTests !== 120) {
    fail("the project correction count changed after cyclic fragment re-entry");
  }
  if (!(progress.tripletSummariesReused > 0) ||
      !(progress.methodScansSkipped > 0) ||
      !(progress.cachedSignalsReused > 0) ||
      !(progress.tripletKernelEvaluations < progress.cumulativeTriplets)) {
    fail("cyclic shortlist summaries were not reused across rounds");
  }
  if (engine._rdp_reconcile(cyclicContext) !== 1) {
    fail(`cyclic-shortlist reconciliation failed: ${engine.UTF8ToString(
      engine._rdp_get_error(cyclicContext),
    )}`);
  }
  const results = JSON.parse(engine.UTF8ToString(
    engine._rdp_get_results_json(cyclicContext),
  ));
  const selectedResult = {
    events: results.events.map((event) => [
      event.recombinant,
      event.majorParent,
      event.minorParent,
      event.beginning,
      event.ending,
      event.wrapsOrigin,
      event.detectionRound,
      event.supportSignalIds,
    ]),
    signals: results.signals.map((signal) => [
      signal.method,
      signal.triplet,
      signal.recombinant,
      signal.majorParent,
      signal.minorParent,
      signal.beginning,
      signal.ending,
      signal.correctedPValue,
      signal.eventId,
    ]),
  };
  const selectedResultDigest = createHash("sha256")
    .update(JSON.stringify(selectedResult))
    .digest("hex");
  if (selectedResultDigest !==
      "5ad90dbeeecd3ea531d52455dd3ded89498c8d0aeefc5d73c2885e451648e6fa") {
    fail(`cyclic-shortlist selected results changed (${selectedResultDigest})`);
  }

  // A user stop during a later cyclic pass is a graceful workflow boundary:
  // discard the unfinished pass, preserve its already committed event prefix,
  // and continue through reconciliation so Review receives usable results.
  const stoppedContext = engine._rdp_create();
  if (!stoppedContext) fail("the WASM engine could not create the graceful-stop context");
  try {
    if (engine._rdp_load_alignment(
      stoppedContext,
      syntheticPointer,
      syntheticFasta.byteLength,
    ) !== 1) {
      fail(`graceful-stop alignment load failed: ${engine.UTF8ToString(
        engine._rdp_get_error(stoppedContext),
      )}`);
    }
    const stopScanStarted = engine._rdp_scan_begin(
      stoppedContext,
      1, 1, 0.05, 30,
      0, 70,
      0, 60,
      0, 1, 1,
      0,
      0, 200, 20, 100, 0.7, 3,
      0,
      0,
      referenceGroupsPointer,
      syntheticSequences.length,
      zeroFlagsPointer,
      syntheticSequences.length,
      zeroFlagsPointer,
      syntheticSequences.length,
    );
    if (stopScanStarted !== 1) {
      fail(`graceful-stop scan could not start: ${engine.UTF8ToString(
        engine._rdp_get_error(stoppedContext),
      )}`);
    }
    let stopStatus = 0;
    for (let batch = 0; batch < 10000 && stopStatus !== 4; ++batch) {
      stopStatus = engine._rdp_scan_batch(stoppedContext, 10000);
      if (stopStatus < 0) {
        fail(`graceful-stop opening round failed: ${engine.UTF8ToString(
          engine._rdp_get_error(stoppedContext),
        )}`);
      }
    }
    const beforeStop = JSON.parse(engine.UTF8ToString(
      engine._rdp_get_progress_json(stoppedContext),
    ));
    if (stopStatus !== 4 || beforeStop.eventCount < 1 || beforeStop.scanRound < 2) {
      fail("graceful-stop fixture did not enter a later cyclic round");
    }
    stopStatus = engine._rdp_scan_batch(stoppedContext, 5);
    if (stopStatus !== 0) {
      fail("graceful-stop fixture did not pause inside the unfinished cyclic round");
    }
    engine._rdp_cancel(stoppedContext);
    stopStatus = engine._rdp_scan_batch(stoppedContext, 5);
    const stoppedProgress = JSON.parse(engine.UTF8ToString(
      engine._rdp_get_progress_json(stoppedContext),
    ));
    if (stopStatus !== 3 || stoppedProgress.phase !== "reconciliation" ||
        stoppedProgress.cycleTermination !== "user-stopped" ||
        stoppedProgress.eventCount !== beforeStop.eventCount) {
      fail("a cyclic stop did not preserve completed events for reconciliation");
    }
    if (engine._rdp_reconcile(stoppedContext) !== 1) {
      fail(`graceful-stop reconciliation failed: ${engine.UTF8ToString(
        engine._rdp_get_error(stoppedContext),
      )}`);
    }
    const stoppedResults = JSON.parse(engine.UTF8ToString(
      engine._rdp_get_results_json(stoppedContext),
    ));
    if (stoppedResults.cycleTermination !== "user-stopped" ||
        stoppedResults.events.length !== beforeStop.eventCount ||
        stoppedResults.signals.some((signal) => signal.eventId < 0)) {
      fail("graceful-stop results contain an unfinished-round signal or lost event");
    }
  } finally {
    engine._rdp_destroy(stoppedContext);
  }
} finally {
  engine._free(syntheticPointer);
  engine._free(referenceGroupsPointer);
  engine._free(zeroFlagsPointer);
  engine._rdp_destroy(cyclicContext);
}

await inspectTree(output);
console.log(
  `GitHub Pages artifact verified: ${requiredFiles.length} required files, FASTA upload, cyclic-shortlist, and graceful-stop smoke tests passed, ${totalBytes.toLocaleString()} total bytes.`,
);
