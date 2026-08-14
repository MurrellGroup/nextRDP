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
    0, 0, 200, 20, 100, 0.7, 3, // primary/secondary BootScan
    0, 0, 200, 20, 100, 1000, 3, // primary/secondary SISCAN
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
  if (status !== 3 || progress.eventCount < 2 || progress.scanRound < 3 ||
      progress.processedTriplets !== progress.totalTriplets) {
    fail("cyclic-shortlist smoke test did not complete two detection rounds");
  }
  if (progress.correctionTests !== 120) {
    fail("the project correction count changed after cyclic fragment re-entry");
  }
  if (!(progress.tripletSummariesReused > 0) ||
      !(progress.cleanTripletsPruned > 0) ||
      !(progress.methodScansSkipped > 0) ||
      !(progress.cachedSignalsReused > 0) ||
      !(progress.invalidScheduleTripletsSkipped > 0) ||
      !(progress.fragmentSequencesPruned > 0) ||
      !(progress.tripletKernelEvaluations < progress.cumulativeTriplets)) {
    fail("cyclic shortlist or event-free fragment pruning was not applied across rounds");
  }
  if (engine._rdp_reconcile(cyclicContext) !== 1) {
    fail(`cyclic-shortlist reconciliation failed: ${engine.UTF8ToString(
      engine._rdp_get_error(cyclicContext),
    )}`);
  }
  const resultsText = engine.UTF8ToString(
    engine._rdp_get_results_json(cyclicContext),
  );
  const results = JSON.parse(resultsText);
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

  // Reconstruct the supplied PHYLPRO diagnostic lazily from the reconciled
  // event. It is deliberately a review-only profile: the supplied RDP5 route
  // has no active significance test and therefore must not invent a p-value
  // or feed a new detection back into the cyclic scheduler.
  const phylproPointer = engine._rdp_get_event_phylpro_json(
    cyclicContext,
    results.events[0].id,
    60,
    0, // ignore missing observations pairwise
    0, // exclude the zero-distance self observation
  );
  if (!phylproPointer) {
    fail(`PHYLPRO lazy review failed: ${engine.UTF8ToString(
      engine._rdp_get_error(cyclicContext),
    )}`);
  }
  const phylpro = JSON.parse(engine.UTF8ToString(phylproPointer));
  if (phylpro.status !== "source-shaped-active-unvalidated" ||
      phylpro.kernel !== "FindSubSeqPP-MakePDstMat-UpdatePDstMat-PPRegression" ||
      phylpro.significanceTest !== "not-implemented-in-supplied-rdp5" ||
      phylpro.optimization !== "three-target-rows-linear-in-context" ||
      phylpro.circular !== true ||
      phylpro.contextSequences !== syntheticSequences.length ||
      phylpro.sequenceIndices.length !== 3 ||
      phylpro.minimumBySequence.length !== 3 ||
      phylpro.breakpoints.length !== 2 ||
      !(phylpro.eligibleColumns > 0) ||
      !(phylpro.evaluatedPoints >= phylpro.returnedPoints) ||
      !Array.isArray(phylpro.points) || phylpro.points.length < 2 ||
      !phylpro.points.every((point) =>
        Number.isFinite(point.recombinant) &&
        Number.isFinite(point.majorParent) &&
        Number.isFinite(point.minorParent))) {
    fail("PHYLPRO lazy review profile/provenance is incomplete");
  }
  const resultsAfterPhylpro = engine.UTF8ToString(
    engine._rdp_get_results_json(cyclicContext),
  );
  if (resultsAfterPhylpro !== resultsText) {
    fail("PHYLPRO review mutated the reconciled discovery result");
  }

  // Exercise primary distance-mode BootScan through the public WASM ABI. The
  // ten-record fixture contains two mosaics and 45 unique sequence pairs, so
  // a complete opening round must discover support regions while reusing the
  // same pair/window/bootstrap summaries across different triplets.
  const bootscanContext = engine._rdp_create();
  if (!bootscanContext) fail("the WASM engine could not create the primary-BootScan context");
  try {
    if (engine._rdp_load_alignment(
      bootscanContext,
      syntheticPointer,
      syntheticFasta.byteLength,
    ) !== 1) {
      fail(`primary-BootScan alignment load failed: ${engine.UTF8ToString(
        engine._rdp_get_error(bootscanContext),
      )}`);
    }
    const bootscanStarted = engine._rdp_scan_begin(
      bootscanContext,
      1, 1, 0.05, 30,
      0, 70,
      0, 60,
      0, 1, 1,
      0,
      1, 0, 100, 20, 30, 0.7, 3,
      0, 0, 200, 20, 100, 1000, 3,
      0,
      0,
      referenceGroupsPointer,
      syntheticSequences.length,
      zeroFlagsPointer,
      syntheticSequences.length,
      zeroFlagsPointer,
      syntheticSequences.length,
    );
    if (bootscanStarted !== 1) {
      fail(`primary-BootScan scan could not start: ${engine.UTF8ToString(
        engine._rdp_get_error(bootscanContext),
      )}`);
    }
    let bootscanStatus = 0;
    for (let batch = 0;
      batch < 10000 && bootscanStatus !== 4 && bootscanStatus !== 3;
      ++batch) {
      bootscanStatus = engine._rdp_scan_batch(bootscanContext, 10000);
      if (bootscanStatus < 0) {
        fail(`primary-BootScan opening round failed: ${engine.UTF8ToString(
          engine._rdp_get_error(bootscanContext),
        )}`);
      }
    }
    const bootscanProgress = JSON.parse(engine.UTF8ToString(
      engine._rdp_get_progress_json(bootscanContext),
    ));
    if (bootscanStatus !== 4 || bootscanProgress.eventCount < 1) {
      fail("primary-BootScan fixture did not commit an opening-round event");
    }
    if (!(bootscanProgress.bootscanProfilesScanned > 0) ||
        !(bootscanProgress.bootscanCandidateRegionsScored > 0) ||
        !(bootscanProgress.bootscanCandidatesFound > 0) ||
        !(bootscanProgress.bootscanPairProfilesRequested > 0) ||
        !(bootscanProgress.bootscanPairProfileCacheHits > 0) ||
        !(bootscanProgress.bootscanPairProfileCacheMisses > 0) ||
        bootscanProgress.bootscanPairProfileCacheHits +
          bootscanProgress.bootscanPairProfileCacheMisses !==
          bootscanProgress.bootscanPairProfilesRequested ||
        !(bootscanProgress.bootscanPairProfileCachePeakBytes > 0) ||
        bootscanProgress.bootscanPairProfileCachePeakBytes > 64 * 1024 * 1024) {
      fail("primary-BootScan discovery/cache telemetry is inconsistent");
    }

    // Stop at the new round boundary so the test retains the completed event
    // but avoids spending CI time rediscovering the second mosaic.
    engine._rdp_cancel(bootscanContext);
    bootscanStatus = engine._rdp_scan_batch(bootscanContext, 1);
    if (bootscanStatus !== 3 || engine._rdp_reconcile(bootscanContext) !== 1) {
      fail(`primary-BootScan result finalization failed: ${engine.UTF8ToString(
        engine._rdp_get_error(bootscanContext),
      )}`);
    }
    const bootscanResults = JSON.parse(engine.UTF8ToString(
      engine._rdp_get_results_json(bootscanContext),
    ));
    const bootscanSignal = bootscanResults.signals.find(
      (signal) => signal.method === "BOOTSCAN" && signal.bootscanDiscovery,
    );
    if (!bootscanResults.bootscanPrimaryEnabled || !bootscanSignal ||
        bootscanSignal.bootscanDiscovery.strictClosestPairVoting !== true ||
        bootscanSignal.bootscanDiscovery.probabilityModel !== "MakeScoresBS-binomial" ||
        !(bootscanSignal.bootscanDiscovery.rawPValue > 0) ||
        !(bootscanSignal.bootscanDiscovery.correctedPValue > 0)) {
      fail("primary-BootScan evidence did not survive event reconciliation");
    }
    const bootscanPlot = JSON.parse(engine.UTF8ToString(
      engine._rdp_get_signal_plot_json(bootscanContext, bootscanSignal.id),
    ));
    if (bootscanPlot.method !== "BOOTSCAN" ||
        bootscanPlot.metric !== "bootstrap-support" ||
        !Array.isArray(bootscanPlot.points) || bootscanPlot.points.length < 2) {
      fail("primary-BootScan review plot was not reconstructed");
    }
  } finally {
    engine._rdp_destroy(bootscanContext);
  }

  // Exercise the supplied SISCAN path through the production Emscripten ABI.
  // The fourth sequence is the deterministic nearest outlier; the recombinant
  // switches from parent one to parent two for a planted 60-site tract.
  const siscanLength = 240;
  const siscanAlphabet = "ACGT";
  const siscanParentOne = Array.from(
    { length: siscanLength },
    (_, index) => siscanAlphabet[index % 4],
  ).join("");
  const siscanParentTwo = Array.from(
    { length: siscanLength },
    (_, index) => siscanAlphabet[(index + 1) % 4],
  ).join("");
  const siscanOutlier = Array.from(
    { length: siscanLength },
    (_, index) => siscanAlphabet[(index + 2) % 4],
  ).join("");
  const siscanRecombinant = Array.from(
    siscanParentOne,
    (base, index) => index >= 80 && index <= 139 ? siscanParentTwo[index] : base,
  ).join("");
  const siscanFasta = new TextEncoder().encode(
    `>siscan-recombinant\n${siscanRecombinant}\n` +
    `>siscan-parent-one\n${siscanParentOne}\n` +
    `>siscan-parent-two\n${siscanParentTwo}\n` +
    `>siscan-outlier\n${siscanOutlier}\n`,
  );
  const siscanPointer = engine._malloc(siscanFasta.byteLength);
  const siscanContext = engine._rdp_create();
  if (!siscanPointer || !siscanContext) {
    fail("the WASM engine could not allocate the SISCAN smoke test");
  }
  try {
    engine.HEAPU8.set(siscanFasta, siscanPointer);
    if (engine._rdp_load_alignment(
      siscanContext,
      siscanPointer,
      siscanFasta.byteLength,
    ) !== 1) {
      fail(`SISCAN alignment load failed: ${engine.UTF8ToString(
        engine._rdp_get_error(siscanContext),
      )}`);
    }
    const siscanStarted = engine._rdp_scan_begin(
      siscanContext,
      0, 1, 0.05, 30,
      0, 70,
      0, 60,
      0, 1, 1,
      0,
      0, 0, 200, 20, 100, 0.7, 3,
      1, 0, 40, 5, 20, 100, 3,
      0,
      0,
      referenceGroupsPointer,
      4,
      zeroFlagsPointer,
      4,
      zeroFlagsPointer,
      4,
    );
    if (siscanStarted !== 1) {
      fail(`SISCAN scan could not start: ${engine.UTF8ToString(
        engine._rdp_get_error(siscanContext),
      )}`);
    }
    let siscanStatus = 0;
    for (let batch = 0;
      batch < 1000 && siscanStatus !== 4 && siscanStatus !== 3;
      ++batch) {
      siscanStatus = engine._rdp_scan_batch(siscanContext, 64);
      if (siscanStatus < 0) {
        fail(`SISCAN scan failed: ${engine.UTF8ToString(
          engine._rdp_get_error(siscanContext),
        )}`);
      }
    }
    const siscanProgress = JSON.parse(engine.UTF8ToString(
      engine._rdp_get_progress_json(siscanContext),
    ));
    if (siscanStatus !== 4 || siscanProgress.eventCount < 1 ||
        !(siscanProgress.siscanProfilesScanned > 0) ||
        !(siscanProgress.siscanWindowsScored > 0) ||
        !(siscanProgress.siscanCandidateRegionsScored > 0) ||
        !(siscanProgress.siscanCandidatesFound > 0) ||
        !(siscanProgress.siscanPermutationDraws > 0) ||
        siscanProgress.siscanContextBuilds !== 1 ||
        siscanProgress.siscanContextTreeMerges !== 3 ||
        !(siscanProgress.siscanRandomValuesGenerated > 0)) {
      fail("SISCAN discovery/context/random-prefix telemetry is inconsistent");
    }
    engine._rdp_cancel(siscanContext);
    siscanStatus = engine._rdp_scan_batch(siscanContext, 1);
    if (siscanStatus !== 3 || engine._rdp_reconcile(siscanContext) !== 1) {
      fail(`SISCAN result finalization failed: ${engine.UTF8ToString(
        engine._rdp_get_error(siscanContext),
      )}`);
    }
    const siscanResults = JSON.parse(engine.UTF8ToString(
      engine._rdp_get_results_json(siscanContext),
    ));
    const siscanSignal = siscanResults.signals.find(
      (signal) => signal.method === "SISCAN" && signal.siscanDiscovery,
    );
    if (!siscanResults.siscanPrimaryEnabled || !siscanSignal ||
        siscanSignal.siscanDiscovery.outlierMode !== "nearest-source-wpgma" ||
        siscanSignal.siscanDiscovery.permutationGenerator !== "microsoft-crt-flat-prefix" ||
        siscanSignal.siscanDiscovery.sourceFastWindowQuirk !== true ||
        !(siscanSignal.siscanDiscovery.maximumZ > 0) ||
        !(siscanSignal.siscanDiscovery.normalTailPValue > 0) ||
        !(siscanSignal.siscanDiscovery.windowAdjustedPValue > 0) ||
        !(siscanSignal.siscanDiscovery.correctedPValue > 0)) {
      fail("SISCAN discovery evidence did not survive reconciliation");
    }
    const siscanPlot = JSON.parse(engine.UTF8ToString(
      engine._rdp_get_signal_plot_json(siscanContext, siscanSignal.id),
    ));
    if (siscanPlot.method !== "SISCAN" ||
        siscanPlot.metric !== "sister-scan-z-score" ||
        !(siscanPlot.minimumValue < 0) ||
        !(siscanPlot.maximumValue > 0) ||
        !Array.isArray(siscanPlot.points) || siscanPlot.points.length < 2) {
      fail("SISCAN review plot was not reconstructed");
    }
  } finally {
    engine._free(siscanPointer);
    engine._rdp_destroy(siscanContext);
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
      0, 0, 200, 20, 100, 0.7, 3,
      0, 0, 200, 20, 100, 1000, 3,
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
  `GitHub Pages artifact verified: ${requiredFiles.length} required files, FASTA upload, primary-BootScan/cache, SISCAN/context/random-prefix, cyclic-shortlist, PHYLPRO review, and graceful-stop smoke tests passed, ${totalBytes.toLocaleString()} total bytes.`,
);
