import { readFile } from "node:fs/promises";
import { resolve } from "node:path";

const root = process.cwd();
const read = (path) => readFile(resolve(root, path), "utf8");
const [
  packageSource,
  lockSource,
  cmake,
  header,
  implementation,
  method,
  burt,
  maxchi,
  chimaeraHeader,
  geneconvHeader,
  geneconv,
  threeSeqHeader,
  threeSeq,
  worker,
  types,
  client,
  app,
  styles,
  datasetStep,
  settings,
  scan,
  review,
  signalPlot,
  eventAlignmentInspector,
  eventTreeInspector,
  exportStep,
  pagesWorkflow,
  viteConfig,
  pagesVerifier,
  sessionHandoff,
  session17Handoff,
  session18Handoff,
  chimaeraTrace,
  geneconvTrace,
  threeSeqTrace,
  lateTrace,
  shortlistTrace,
  readme,
  status,
] = await Promise.all([
  read("package.json"),
  read("package-lock.json"),
  read("wasm/CMakeLists.txt"),
  read("wasm/include/rdp_api.h"),
  read("wasm/src/rdp_api.cpp"),
  read("wasm/src/rdp_method.cpp"),
  read("wasm/src/burt_confidence.cpp"),
  read("wasm/src/maxchi.cpp"),
  read("wasm/src/chimaera.hpp"),
  read("wasm/src/geneconv.hpp"),
  read("wasm/src/geneconv.cpp"),
  read("wasm/src/threeseq.hpp"),
  read("wasm/src/threeseq.cpp"),
  read("src/workers/analysis.worker.ts"),
  read("src/lib/types.ts"),
  read("src/lib/wasmClient.ts"),
  read("src/App.tsx"),
  read("src/styles.css"),
  read("src/components/DatasetStep.tsx"),
  read("src/components/SettingsStep.tsx"),
  read("src/components/ScanStep.tsx"),
  read("src/components/ReviewStep.tsx"),
  read("src/components/SignalPlot.tsx"),
  read("src/components/EventAlignmentInspector.tsx"),
  read("src/components/EventTreeInspector.tsx"),
  read("src/components/ExportStep.tsx"),
  read(".github/workflows/deploy-pages.yml"),
  read("vite.config.ts"),
  read("scripts/verify-pages-output.mjs"),
  read("docs/session-15-handoff.md"),
  read("docs/session-17-handoff.md"),
  read("docs/session-18-handoff.md"),
  read("docs/native-chimaera-discovery-trace.md"),
  read("docs/native-geneconv-discovery-trace.md"),
  read("docs/native-threeseq-discovery-trace.md"),
  read("docs/native-late-consensus-trace.md"),
  read("docs/native-cyclic-shortlist-trace.md"),
  read("README.md"),
  read("STATUS.md"),
]);

const packageMetadata = JSON.parse(packageSource);
const lockMetadata = JSON.parse(lockSource);

function fail(message) {
  throw new Error(`Source contract check failed: ${message}`);
}

function captures(source, expression) {
  return [...source.matchAll(expression)].map((match) => match[1]);
}

function uniqueSorted(values) {
  return [...new Set(values)].sort();
}

function compareSets(label, expected, actual) {
  const expectedSet = new Set(expected);
  const actualSet = new Set(actual);
  const missing = [...expectedSet].filter((value) => !actualSet.has(value));
  const extra = [...actualSet].filter((value) => !expectedSet.has(value));
  if (missing.length || extra.length) {
    fail(
      `${label}; missing [${missing.join(", ")}], extra [${extra.join(", ")}]`,
    );
  }
}

const headerFunctions = uniqueSorted(
  captures(
    header,
    /\b(?:const\s+char\s*\*|void|int|uint32_t)\s+(rdp_[a-z0-9_]+)\s*\(/g,
  ),
);
const implementationFunctions = uniqueSorted(
  captures(
    implementation,
    /RDP_KEEPALIVE\s+(?:const\s+char\s*\*|[\w:]+)\s+(rdp_[a-z0-9_]+)\s*\(/g,
  ),
);
compareSets("rdp_api.h and rdp_api.cpp exports differ", headerFunctions, implementationFunctions);

function parameterCount(source, functionName) {
  const match = new RegExp(`\\b${functionName}\\s*\\(`).exec(source);
  if (!match) fail(`could not find ${functionName} while checking ABI arity`);
  const opening = match.index + match[0].lastIndexOf("(");
  let depth = 1;
  let commas = 0;
  let nonWhitespace = false;
  for (let index = opening + 1; index < source.length; index += 1) {
    const character = source[index];
    if (character === "(") depth += 1;
    else if (character === ")") {
      depth -= 1;
      if (depth === 0) {
        const parameters = source.slice(opening + 1, index).trim();
        if (parameters === "void") return 0;
        return nonWhitespace ? commas + (/,$/.test(parameters) ? 0 : 1) : 0;
      }
    } else if (character === "," && depth === 1) {
      commas += 1;
    } else if (depth === 1 && !/\s/.test(character)) {
      nonWhitespace = true;
    }
  }
  fail(`unterminated parameter list for ${functionName}`);
}

for (const functionName of headerFunctions) {
  const declared = parameterCount(header, functionName);
  const implemented = parameterCount(implementation, functionName);
  if (declared !== implemented) {
    fail(`${functionName} ABI arity differs: header ${declared}, implementation ${implemented}`);
  }
  if (worker.includes(`_${functionName}(`)) {
    const bridged = parameterCount(worker, `_${functionName}`);
    if (declared !== bridged) {
      fail(`${functionName} worker ABI arity differs: C ${declared}, worker ${bridged}`);
    }
  }
}

const cmakeExports = uniqueSorted(
  captures(cmake, /'(_(?:rdp_[a-z0-9_]+|malloc|free))'/g),
);
const expectedCmakeExports = uniqueSorted([
  "_malloc",
  "_free",
  ...headerFunctions.map((name) => `_${name}`),
]);
compareSets("CMake exported function list differs from rdp_api.h", expectedCmakeExports, cmakeExports);

if (!cmake.includes("'UTF8ToString','HEAPU8'") ||
    !worker.includes("module.HEAPU8.set(bytes, pointer)")) {
  fail("the worker FASTA upload path does not have an exported WASM heap view");
}

const workerFunctions = uniqueSorted(
  captures(worker, /\b_(rdp_[a-z0-9_]+)\b/g).map((name) => `_${name}`),
);
const missingWorkerExports = workerFunctions.filter((name) => !cmakeExports.includes(name));
if (missingWorkerExports.length) {
  fail(`worker calls unexported ABI functions [${missingWorkerExports.join(", ")}]`);
}

const version = packageMetadata.version;
if (typeof version !== "string" || !version) fail("package version is missing");
if (lockMetadata.version !== version || lockMetadata.packages?.[""]?.version !== version) {
  fail("package.json and package-lock.json versions differ");
}
if (!implementation.includes(`\\\"engineVersion\\\":\\\"${version}\\\"`) ||
    !implementation.includes(`return \"${version}\";`)) {
  fail("rdp_api.cpp engine version differs from package.json");
}
if (!method.includes(`\\\"engineVersion\\\":\\\"${version}\\\"`)) {
  fail("rdp_method.cpp result version differs from package.json");
}

const schema = "org.rdp-web.project/v1alpha16";
if (!implementation.includes(schema) || !worker.includes(schema)) {
  fail(`emitted/imported project schema ${schema} is not aligned`);
}
for (let generation = 1; generation <= 16; generation += 1) {
  if (!worker.includes(`schema !== "org.rdp-web.project/v1alpha${generation}"`)) {
    fail(`v1alpha${generation} project import is missing`);
  }
}
if (!worker.includes('schema === "org.rdp-web.project/v1alpha10" || supportsReferenceGroups') ||
    !worker.includes("supportsMaxChiDiscovery && analysis.maxChiEnabled !== false")) {
  fail("older project imports do not preserve their RDP-only discovery semantics");
}
if (!worker.includes('schema === "org.rdp-web.project/v1alpha12" ||') ||
    !worker.includes("supportsChimaeraDiscovery && analysis.chimaeraEnabled !== false")) {
  fail("pre-v12 project imports do not preserve their pre-CHIMAERA discovery semantics");
}
if (!worker.includes("supportsGeneconvDiscovery") ||
    !worker.includes('schema === "org.rdp-web.project/v1alpha13" ||') ||
    !worker.includes("supportsGeneconvDiscovery && analysis.geneconvEnabled !== false")) {
  fail("pre-v13 project imports do not preserve their pre-GENECONV discovery semantics");
}
if (!worker.includes("supportsThreeSeqDiscovery") ||
    !worker.includes('schema === "org.rdp-web.project/v1alpha14" ||') ||
    !worker.includes("supportsThreeSeqDiscovery && analysis.threeSeqEnabled !== false")) {
  fail("pre-v14 project imports do not preserve their pre-3SEQ discovery semantics");
}
if (!worker.includes("supportsThreeSeqSplit") ||
    !worker.includes('schema === "org.rdp-web.project/v1alpha15" ||') ||
    !worker.includes('schema === "org.rdp-web.project/v1alpha16"') ||
    !worker.includes("claims post-erasure split evidence in a pre-v15 project")) {
  fail("pre-v15 project imports do not reject unsupported 3SEQ split evidence");
}
if (!worker.includes("assetVersion") || !worker.includes("loadedVersion !== assetVersion")) {
  fail("worker asset/engine version guard is missing");
}
if (!method.includes("triplet_contains(signal.recombinant)") ||
    !method.includes("triplet_contains(signal.major_parent)") ||
    !method.includes("triplet_contains(signal.minor_parent)") ||
    !method.includes("three distinct triplet members")) {
  fail("saved signal restore does not enforce distinct role membership in its triplet");
}

for (const workflowContract of [
  "actions/checkout@v7",
  "actions/setup-node@v7",
  "emscripten-core/setup-emsdk@v16",
  "actions/configure-pages@v6",
  "actions/upload-pages-artifact@v5.0.0",
  "include-hidden-files: true",
  "actions/deploy-pages@v5",
  "pages: write",
  "id-token: write",
  "node-version: 20",
  "version: 5.0.1",
  "npm ci",
  "npm run check:source",
  "npm run check:types",
  "npm run build",
  "path: dist",
]) {
  if (!pagesWorkflow.includes(workflowContract)) {
    fail(`GitHub Pages workflow is missing ${workflowContract}`);
  }
}
if (!pagesWorkflow.includes("github.event.repository.default_branch")) {
  fail("GitHub Pages workflow does not gate automatic deployment to the default branch");
}
if (!viteConfig.includes('base: "./"')) {
  fail("Vite base is not relative for repository-subpath Pages deployment");
}
for (const artifactContract of [
  '".nojekyll"',
  '"wasm/rdp-core.mjs"',
  '"wasm/rdp-core.wasm"',
  "WebAssembly header",
  "FASTA upload smoke test",
  "cyclic-shortlist",
  "engine.HEAPU8.set(fasta, fastaPointer)",
  "is a symbolic link",
  "is hard-linked",
]) {
  if (!pagesVerifier.includes(artifactContract)) {
    fail(`Pages artifact verifier is missing ${artifactContract}`);
  }
}

const activeLateConsensusStatus =
  "active-rdp-maxchi-chimaera-geneconv-threeseq-plus-optional-bootscan-post-group-recheck";
if (!method.includes(activeLateConsensusStatus) || !types.includes(activeLateConsensusStatus)) {
  fail("active late-consensus status differs between core and web contract");
}
if (!method.includes("nativeGroupMembershipComplete") ||
    !types.includes("nativeGroupMembershipComplete: true")) {
  fail("native group-membership completion flag differs between core and web contract");
}
if (!method.includes("nativePrimaryRdpRecheckComplete") ||
    !types.includes("nativePrimaryRdpRecheckComplete: true")) {
  fail("primary-RDP post-group recheck flag differs between core and web contract");
}
for (const statusContract of [
  "maxChiTripletRecheckApplied",
  "maxChiPostGroupRecheckApplied",
  "source-shaped-strongest-peak-unvalidated",
  "maxChiEventDiscoveryApplied",
  "maxChiDiscoveryFeedsCyclicScheduler",
  "source-shaped-multi-peak-destroy-retry-unvalidated",
  "nativeMaxChiFullRecheckComplete",
]) {
  if (!method.includes(statusContract) || !types.includes(statusContract)) {
    fail(`MaxChi late-consensus status is missing ${statusContract}`);
  }
}
for (const statusContract of [
  "bootscanSecondaryEnabled",
  "bootscanTripletRecheckApplied",
  "bootscanPostGroupRecheckApplied",
  "source-shaped-distance-bootstrap-binomial-unvalidated",
  "nativeBootscanFullRecheckComplete",
]) {
  if (!method.includes(statusContract) || !types.includes(statusContract)) {
    fail(`BootScan secondary-recheck status is missing ${statusContract}`);
  }
}
for (const statusContract of [
  "threeSeqKernelStatus",
  "source-shaped-hypergeometric-random-walk-unvalidated",
  "threeSeqEventDiscoveryApplied",
  "threeSeqDiscoveryFeedsCyclicScheduler",
  "threeSeqTripletRecheckApplied",
  "threeSeqPostGroupRecheckApplied",
  "threeSeqRecheckKernelStatus",
  "source-shaped-findall-two-orientation-unvalidated",
  "nativeThreeSeqFullRecheckComplete",
]) {
  if (!method.includes(statusContract) || !types.includes(statusContract)) {
    fail(`3SEQ late-consensus status is missing ${statusContract}`);
  }
}
for (const statusContract of [
  "geneconvKernelStatus",
  "source-shaped-six-track-ka-fragments-unvalidated",
  "geneconvEventDiscoveryApplied",
  "geneconvDiscoveryFeedsCyclicScheduler",
  "geneconvTripletRecheckApplied",
  "geneconvPostGroupRecheckApplied",
  "geneconvRecheckKernelStatus",
  "source-shaped-six-track-best-fragment-unvalidated",
  "nativeGeneconvFullRecheckComplete",
]) {
  if (!method.includes(statusContract) || !types.includes(statusContract)) {
    fail(`GENECONV late-consensus status is missing ${statusContract}`);
  }
}
for (const statusContract of [
  "chimaeraKernelStatus",
  "source-shaped-target-profile-multi-peak-destroy-retry-unvalidated",
  "chimaeraEventDiscoveryApplied",
  "chimaeraDiscoveryFeedsCyclicScheduler",
  "chimaeraTripletRecheckApplied",
  "chimaeraPostGroupRecheckApplied",
  "source-shaped-three-target-strongest-peak-unvalidated",
  "nativeChimaeraFullRecheckComplete",
]) {
  if (!method.includes(statusContract) || !types.includes(statusContract)) {
    fail(`CHIMAERA late-consensus status is missing ${statusContract}`);
  }
}
for (const field of [
  "finalTrimFirstExpansionAdded",
  "finalTrimSecondExpansionAdded",
  "selectedRolePrunedOut",
  "finalTrimMember",
  "consensusPrimaryMember",
  "consensusEquivalentMember",
  "consensusStragglerMember",
  "consensusRebuiltMember",
  "consensusFallbackRestored",
  "selectedTreeCleanupPrunedOut",
  "selectedTreeCleanupAdded",
  "finalDistanceMember",
  "postGroupRdpRecheck",
  "localPValueCutoff",
  "emittedSignalCount",
  "candidateSignalCount",
  "overlappingSignalCount",
  "eventRedetected",
  "bestCorrectedPValue",
  "postGroupMaxChiRecheck",
  "postGroupChimaeraRecheck",
  "postGroupGeneconvRecheck",
  "postGroupThreeSeqRecheck",
]) {
  if (!method.includes(`\\\"${field}\\\"`) || !types.includes(field)) {
    fail(`late-consensus field ${field} differs between core and web contract`);
  }
}

for (const field of [
  "erasureWithinRdpWindow",
  "uncertainDueToErasure",
  "rdpWindowInformativeSites",
  "nearestErasureInformativeSites",
  "uncertainPriorEventIds",
  "uncertainErasureEventIds",
  "nativeCheckEndsApplied",
  "nativeCheckEndsWarning",
  "informationProfileAvailable",
  "inputMissingDataInCheckRange",
  "linearEdgeWithinRdpWindow",
  "nativeCheckRange",
  "nativeCheckEndsStatus",
  "nativeCheckEndsAfterFirstEvent",
  "inputMissingRunLength",
  "uncertaintyReasons",
  "statisticalConfidence",
  "complete-active-unvalidated",
  "BURT/BenHMM",
  "hmmCyclesArgument",
  "serialTrainingStarts",
  "posteriorThresholds",
  "randomAdapter",
  "optionAvailable",
  "enabledByDefault",
  "canRepositionDetectedEvents",
  "breakpointConfidence",
  "candidateIntervalCount",
  "bestLogLikelihood",
  "confidence99",
  "hmmCoordinate",
  "confidence95",
  "insufficientInsideOrOutsideReverted",
]) {
  if (!method.includes(`\\\"${field}\\\"`) || !types.includes(field)) {
    fail(`breakpoint-inspection field ${field} differs between core and web contract`);
  }
}

if (!header.includes("int polish_breakpoints") ||
    !implementation.includes("int polish_breakpoints") ||
    !method.includes("options_.polish_breakpoints") ||
    !worker.includes("polishBreakpoints: number") ||
    !types.includes("polishBreakpoints: boolean") ||
    !method.includes('\\"polishBreakpoints\\"')) {
  fail("BURT polish-breakpoint option differs across ABI, worker, core, and result contract");
}
if (!app.includes("polishBreakpoints: true") ||
    !settings.includes('checked={options.polishBreakpoints}') ||
    !method.includes('unavailable_reason = "disabled"') ||
    !types.includes('| "disabled"')) {
  fail("BURT default-enabled/disabled-state workflow contract is incomplete");
}

for (const curationContract of [
  "disabledSequenceIndices",
  "disabled_sequences",
  "options.disabled",
]) {
  if (!header.includes(curationContract) &&
      !implementation.includes(curationContract) &&
      !method.includes(curationContract) &&
      !worker.includes(curationContract) &&
      !types.includes(curationContract) &&
      !app.includes(curationContract)) {
    fail(`masked/disabled sequence contract is missing ${curationContract}`);
  }
}
if (!datasetStep.includes('value="enabled"') ||
    !datasetStep.includes('value="masked"') ||
    !datasetStep.includes('value="disabled"') ||
    !method.includes("sequence_masked") ||
    !method.includes("sequence_disabled") ||
    !method.includes("tree_candidates.push_back(static_cast<std::uint32_t>(sequence))") ||
    !method.includes("Disabled sequences cannot be assigned an analysed event role.")) {
  fail("manual enabled/masked/disabled workflow contract is incomplete");
}
for (const bulkStateAction of ["auto-mask", "enable-all", "mask-all", "disable-all"]) {
  if (!app.includes(`\"${bulkStateAction}\"`) ||
      !datasetStep.includes(`\"${bulkStateAction}\"`)) {
    fail(`bulk sequence-curation workflow is missing ${bulkStateAction}`);
  }
}
if (!method.includes("signal_candidates_scratch_") ||
    !method.includes("append_candidate_signals") ||
    method.includes("auto found = candidate_signals")) {
  fail("primary triplet candidate buffer reuse contract is incomplete");
}

if (!cmake.includes("src/burt_confidence.cpp")) {
  fail("BURT/BenHMM implementation is not part of the WASM target");
}
if (!cmake.includes("src/maxchi.cpp")) {
  fail("MaxChi discovery/recheck implementation is not part of the WASM target");
}
if (!cmake.includes("src/geneconv.cpp")) {
  fail("GENECONV discovery implementation is not part of the WASM target");
}
if (!cmake.includes("src/threeseq.cpp")) {
  fail("3SEQ discovery implementation is not part of the WASM target");
}
if (!cmake.includes("_rdp_restore_chimaera_discovery")) {
  fail("CHIMAERA discovery restore is not exported by the WASM target");
}
if (!cmake.includes("_rdp_restore_geneconv_discovery")) {
  fail("GENECONV discovery restore is not exported by the WASM target");
}
if (!cmake.includes("_rdp_restore_threeseq_discovery")) {
  fail("3SEQ discovery restore is not exported by the WASM target");
}
for (const sourceContract of [
  "source_normal_z",
  "source_chi_p_value",
  "source_critical_difference",
  "build_variable_profile",
  "make_banned_windows",
  "strongest_peak",
  "grow_peak",
  "options.p_value_cutoff / 6.0",
  "evidence.within_triplet_p_value",
  "evidence.corrected_p_value",
  "maxchi_discover",
  "maxchi_discover_prepared",
  "source_smooth_chi",
  "constexpr std::ptrdiff_t divisor = q_window * 2 + 1",
  "position <= 1 + q_window",
  "grow_discovery_peak",
  "find_side_chi",
  "optimize_left_breakpoint",
  "optimize_right_breakpoint",
  "destroy_completed_peak_region",
  "destroy_rejected_peak",
  "include_source_peak_in_destroy_region",
  "assign_discovery_roles",
  "std::priority_queue",
  "maximum_peak_attempts",
  "wasted_attempts >= 3",
  "maxchi_plot_profile",
]) {
  if (!maxchi.includes(sourceContract)) {
    fail(`MaxChi source contract is missing ${sourceContract}`);
  }
}
if (!method.includes("kScanGeneconv | kScanMaxchi | kScanChimaera") ||
    !method.includes("kScanThreeseq") ||
    !method.includes("maxchi_discover_prepared(")) {
  fail("combined discovery scan does not reuse its one-pass prepared profile");
}
for (const uiContract of [
  "maxChiTripletRecheck",
  "postGroupMaxChiRecheck",
  "MaxChi recheck",
  "MaxChi peak-to-tract construction",
  "active MCXoverF MaxChi discovery path",
  "Raw chi-square pair profiles",
]) {
  if (!review.includes(uiContract)) {
    fail(`MaxChi review workflow is missing ${uiContract}`);
  }
}
for (const evidenceContract of [
  "maxChiTripletRecheck",
  "missingDataWindowFilterApplied",
  "linearEdgeWindowFilterApplied",
  "bonferroniApplied",
  "correctionTests",
  "variableSites",
  "criticalDifference",
  "grownHalfWindow",
  "peakAlignmentPosition",
  "maximumChiSquare",
  "withinTripletPValue",
  "correctedPValue",
  "sourceRecheckHit",
  "eventDiscoveryApplied",
  "maxChiDiscovery",
  "peakAttempt",
  "peakPair",
  "tractSide",
  "rawPValue",
  "leftFlankChiSquare",
  "rightFlankChiSquare",
  "smoothingUse",
]) {
  if (!method.includes(`\\\"${evidenceContract}\\\"`) || !types.includes(evidenceContract)) {
    fail(`MaxChi evidence contract differs between core and web for ${evidenceContract}`);
  }
}
if (!method.includes("twelve-term-eleven-divisor-source-basin-destruction-only") ||
    !types.includes("twelve-term-eleven-divisor-source-basin-destruction-only")) {
  fail("MaxChi supplied smoothing quirk differs between core and web contract");
}
if (!settings.includes("RDP + GENECONV + MaxChi + CHIMAERA + 3SEQ discovery") ||
    !settings.includes("Discover events with MaxChi") ||
    !settings.includes("Included in event discovery")) {
  fail("method settings do not expose combined RDP/MaxChi discovery");
}
for (const sourceContract of [
  "ChimaeraDiscoveryOptions",
  "ChimaeraRecheckOptions",
  "ChimaeraDiscoveryCandidate",
  "ChimaeraDiscoverySummary",
  "ChimaeraRecheckEvidence",
  "ChimaeraPlotProfile",
  "chimaera_discover_prepared",
  "chimaera_recheck_prepared",
  "chimaera_plot_profile",
]) {
  if (!chimaeraHeader.includes(sourceContract)) {
    fail(`CHIMAERA header contract is missing ${sourceContract}`);
  }
}
for (const sourceContract of [
  "kChimaeraParentOne",
  "kChimaeraParentTwo",
  "kChimaeraScorePair",
  "kChimaeraOtherPair",
  "build_chimaera_target_profile",
  "discover_chimaera_target",
  "options.p_value_cutoff / 6.0",
  "source_smooth_chi",
  "grow_discovery_peak",
  "find_side_chi",
  "optimize_left_breakpoint",
  "optimize_right_breakpoint",
  "destroy_completed_peak_region",
  "destroy_rejected_peak",
  "wasted_attempts >= 3",
  "std::priority_queue",
  "target_profiles_scanned",
  "peak_limit_targets",
  "ChimaeraRecheckEvidence chimaera_recheck_prepared",
  "std::make_tuple",
]) {
  if (!maxchi.includes(sourceContract)) {
    fail(`CHIMAERA source contract is missing ${sourceContract}`);
  }
}
if (!method.includes("chimaera_discover_prepared(") ||
    !method.includes("maxchi_workspace_.triplet_missing_data") ||
    !method.includes("profile.similarities") ||
    !method.includes("SignalMethod::chimaera")) {
  fail("CHIMAERA does not reuse the fused triplet preparation or enter the cyclic signal path");
}
for (const evidenceContract of [
  "chimaeraDiscovery",
  "AlistChi-FastRecCheckChim-CXoverA",
  "target-specific-information-rich-binary-string",
  "raw-chi-square-lazy-heap-per-target",
  "targetLocal",
  "informationRichSites",
  "insideParentOneMatchRate",
  "outsideParentOneMatchRate",
]) {
  if (!method.includes(`\\\"${evidenceContract}\\\"`) || !types.includes(evidenceContract)) {
    fail(`CHIMAERA evidence contract differs between core and web for ${evidenceContract}`);
  }
}
for (const optionContract of [
  "int chimaera_enabled",
  "chimaera_window_sites",
]) {
  if (!header.includes(optionContract) || !implementation.includes(optionContract)) {
    fail(`CHIMAERA scan ABI option is missing ${optionContract}`);
  }
}
for (const workerContract of [
  "chimaeraEnabled: number",
  "chimaeraWindowSites: number",
  "_rdp_restore_chimaera_discovery",
  'savedMethod === "CHIMAERA"',
  "analysis.chimaeraProfilesScanned",
  "analysis.chimaeraPeakAttempts",
  "analysis.chimaeraCandidatesFound",
  "analysis.chimaeraPeakLimitTargets",
]) {
  if (!worker.includes(workerContract)) {
    fail(`CHIMAERA worker/restore contract is missing ${workerContract}`);
  }
}
for (const resultContract of [
  "chimaeraProfilesScanned",
  "chimaeraPeakAttempts",
  "chimaeraCandidatesFound",
  "chimaeraPeakLimitTargets",
  "chimaeraEnabled",
  "chimaeraWindowSites",
]) {
  if (!method.includes(resultContract) || !types.includes(`${resultContract}:`)) {
    fail(`CHIMAERA result contract differs between core and web for ${resultContract}`);
  }
  if (!app.includes(resultContract) || !scan.includes(resultContract)) {
    fail(`CHIMAERA application progress/options surface is missing ${resultContract}`);
  }
}
if (!app.includes("chimaeraEnabled: true") ||
    !app.includes("chimaeraWindowSites: 60") ||
    !settings.includes("Discover events with CHIMAERA") ||
    !settings.includes("target-specific parent-match string")) {
  fail("CHIMAERA discovery defaults/settings are incomplete");
}
for (const uiContract of [
  "CHIMAERA target profile and tract construction",
  "insideParentOneMatchRate",
  "AlistChi → FastRecCheckChim → CXoverA",
  "Target-specific CHIMAERA chi-square profile",
]) {
  if (!review.includes(uiContract)) {
    fail(`CHIMAERA review workflow is missing ${uiContract}`);
  }
}
if (!signalPlot.includes('signal.method === "CHIMAERA"') ||
    !signalPlot.includes("CHIMAERA target χ² profile") ||
    !signalPlot.includes("target-to-parent-one")) {
  fail("CHIMAERA one-trace plot contract is incomplete");
}
if (!method.includes('\\\"maxChiChimaeraOnlySupport\\\"') ||
    !method.includes("MaxChi/CHIMAERA-only support caution") ||
    !types.includes("maxChiChimaeraOnlySupport: boolean") ||
    !review.includes("closely related methods, not independent confirmation")) {
  fail("MaxChi/CHIMAERA related-method support caution is not preserved end to end");
}
for (const recheckContract of [
  "fast_method_triplet_rechecks",
  "chimaera_recheck_prepared(",
  "chimaeraTripletRecheck",
  "postGroupChimaeraRecheck",
  "FastRecCheckChim-three-target-strongest-peak",
  "CHIMAERA triplet recheck status",
]) {
  if (!method.includes(recheckContract)) {
    fail(`CHIMAERA representative/finalized-list recheck is missing ${recheckContract}`);
  }
}
for (const recheckContract of [
  "interface ChimaeraRecheckEvidence",
  "chimaeraTripletRecheck: ChimaeraRecheckEvidence",
  "postGroupChimaeraRecheck: ChimaeraRecheckEvidence",
  "targetProfilesScanned",
  "bestTarget",
]) {
  if (!types.includes(recheckContract)) {
    fail(`CHIMAERA web recheck contract is missing ${recheckContract}`);
  }
}
for (const recheckContract of [
  "CHIMAERA recheck",
  "FastRecCheckChim three-target strongest-peak statistic",
  "postGroupChimaeraRecheck",
  "five late rechecks",
  "related-method evidence rather than two independent confirmations",
]) {
  if (!review.includes(recheckContract)) {
    fail(`CHIMAERA review recheck surface is missing ${recheckContract}`);
  }
}
for (const sourceContract of [
  "GeneconvDiscoveryOptions",
  "GeneconvDiscoveryCandidate",
  "GeneconvDiscoverySummary",
  "GeneconvRecheckEvidence",
  "GeneconvCategoryRun",
  "GeneconvWorkspace",
  "GeneconvPlotProfile",
  "geneconv_discover_prepared",
  "geneconv_recheck_prepared",
  "geneconv_plot_profile",
]) {
  if (!geneconvHeader.includes(sourceContract)) {
    fail(`GENECONV header contract is missing ${sourceContract}`);
  }
}
for (const sourceContract of [
  "track_positive",
  "build_categories",
  "append_source_run",
  "build_source_runs",
  "category_runs",
  "solve_source_lambda_k",
  "source_critical_score",
  "source_ka_probability",
  "-std::log(1.0 - p_value_cutoff)",
  "probability = 1.0 - std::exp(-tail)",
  "std::clamp(probability, 0.0, 1.0)",
  "fragment.raw_p_value < raw_cutoff",
  "static_cast<float>(mismatch_penalty)",
  "prepare_prefix_queries",
  "next_lower_prefix",
  "range_maximum",
  "std::stable_sort",
  "fragment_maximum_coverage",
  "overlap_range_add",
  "geneconv_discover_prepared",
  "geneconv_recheck_prepared",
  "geneconv_plot_profile",
  "std::priority_queue",
]) {
  if (!geneconv.includes(sourceContract)) {
    fail(`GENECONV source contract is missing ${sourceContract}`);
  }
}
if (geneconv.includes("const std::size_t repeat = circular ? 2 : 1")) {
  fail("GENECONV must retain GetFragsP's terminal circular run instead of repeating every signed run");
}
if (geneconv.includes("std::clamp(probability, kMinimumProbability") ||
    geneconv.includes("std::max(kMinimumProbability, probability)")) {
  fail("GENECONV discovery must retain a source-underflowed zero probability");
}
if (!method.includes("geneconv_discover_prepared(") ||
    !method.includes("SignalMethod::geneconv") ||
    !method.includes("profile.similarities") ||
    !method.includes("maxchi_workspace_")) {
  fail("GENECONV does not reuse fused triplet preparation or enter the cyclic signal path");
}
for (const recheckContract of [
  "geneconv_recheck_prepared(",
  "geneconvTripletRecheck",
  "postGroupGeneconvRecheck",
  "GCXoverD-six-track-best-fragment",
  "GENECONV triplet recheck status",
]) {
  if (!method.includes(recheckContract)) {
    fail(`GENECONV representative/finalized-list recheck is missing ${recheckContract}`);
  }
}
for (const recheckContract of [
  "interface GeneconvRecheckEvidence",
  "geneconvTripletRecheck: GeneconvRecheckEvidence",
  "postGroupGeneconvRecheck: GeneconvRecheckEvidence",
  "sourceSkewFilterRejected",
  "overlapRejectedFragments",
]) {
  if (!types.includes(recheckContract)) {
    fail(`GENECONV web recheck contract is missing ${recheckContract}`);
  }
}
for (const recheckContract of [
  "GENECONV recheck",
  "GCXoverD six-track ordinary-kernel recheck",
  "postGroupGeneconvRecheck",
  "five late rechecks",
]) {
  if (!review.includes(recheckContract)) {
    fail(`GENECONV review recheck surface is missing ${recheckContract}`);
  }
}
for (const tieContract of [
  "source_scan_method_priority",
  "case SignalMethod::rdp: return 0",
  "case SignalMethod::geneconv: return 1",
  "case SignalMethod::maxchi: return 2",
  "case SignalMethod::chimaera: return 3",
  "case SignalMethod::threeseq: return 4",
]) {
  if (!method.includes(tieContract)) {
    fail(`Source method-major tie ordering is missing ${tieContract}`);
  }
}
for (const evidenceContract of [
  "geneconvDiscovery",
  "FindSubSeqGCAP6-GetFragsP-GetMaxFragScoreP-CalcKMaxP-GCCalcPValP2-GCXoverD",
  "karlin-altschul",
  "stable-lowest-p-configured-coverage",
  "minimumFragmentFiltersApplied",
  "polymorphicSites",
  "positiveSites",
  "discordantSites",
  "mismatchPenalty",
  "fragmentScore",
  "criticalScore",
  "karlinAltschulK",
  "rawPValue",
  "correctedPValue",
]) {
  if (!method.includes(`\\\"${evidenceContract}\\\"`) || !types.includes(evidenceContract)) {
    fail(`GENECONV evidence contract differs between core and web for ${evidenceContract}`);
  }
}
for (const optionContract of [
  "int geneconv_enabled",
  "geneconv_mismatch_scale",
  "geneconv_max_overlaps",
]) {
  if (!header.includes(optionContract) || !implementation.includes(optionContract)) {
    fail(`GENECONV scan ABI option is missing ${optionContract}`);
  }
}
for (const workerContract of [
  "geneconvEnabled: number",
  "geneconvMismatchScale: number",
  "geneconvMaxOverlaps: number",
  "_rdp_restore_geneconv_discovery",
  'savedMethod === "GENECONV"',
  "analysis.geneconvFragmentsScored",
  "analysis.geneconvQualifiedFragments",
  "analysis.geneconvCandidatesFound",
  "analysis.geneconvOverlapRejections",
  "analysis.geneconvNumericalFallbackTracks",
]) {
  if (!worker.includes(workerContract)) {
    fail(`GENECONV worker/restore contract is missing ${workerContract}`);
  }
}
for (const resultContract of [
  "geneconvFragmentsScored",
  "geneconvQualifiedFragments",
  "geneconvCandidatesFound",
  "geneconvOverlapRejections",
  "geneconvNumericalFallbackTracks",
  "geneconvEnabled",
  "geneconvMismatchScale",
  "geneconvMaxOverlaps",
]) {
  if (!method.includes(resultContract) || !types.includes(`${resultContract}:`)) {
    fail(`GENECONV result contract differs between core and web for ${resultContract}`);
  }
  if (!app.includes(resultContract) || !scan.includes(resultContract)) {
    fail(`GENECONV application progress/options surface is missing ${resultContract}`);
  }
}
if (!app.includes("geneconvEnabled: true") ||
    !app.includes("geneconvMismatchScale: 1") ||
    !app.includes("geneconvMaxOverlaps: 1") ||
    !settings.includes("Discover events with GENECONV") ||
    !settings.includes("stable") ||
    !settings.includes("lowest-P")) {
  fail("GENECONV discovery defaults/settings are incomplete");
}
for (const uiContract of [
  "GENECONV six-track fragment score",
  "Corrected KA hit",
  "Karlin–Altschul",
  "GENECONV Karlin–Altschul fragment envelope",
]) {
  if (!review.includes(uiContract)) {
    fail(`GENECONV review workflow is missing ${uiContract}`);
  }
}
if (!signalPlot.includes('signal.method === "GENECONV"') ||
    !signalPlot.includes("GENECONV negative log10 KA P fragment envelope") ||
    !signalPlot.includes("inner and outer GENECONV fragment envelopes")) {
  fail("GENECONV three-colour fragment plot contract is incomplete");
}
for (const sourceContract of [
  "ThreeSeqDiscoveryOptions",
  "ThreeSeqWalkDirection",
  "ThreeSeqDiscoveryCandidate",
  "ThreeSeqDiscoverySummary",
  "ThreeSeqRecheckEvidence",
  "ThreeSeqPlotProfile",
  "ThreeSeqProbabilityKey",
  "ThreeSeqWorkspace",
  "threeseq_discover_prepared",
  "threeseq_recheck_prepared",
  "threeseq_plot_profile",
]) {
  if (!threeSeqHeader.includes(sourceContract)) {
    fail(`3SEQ header contract is missing ${sourceContract}`);
  }
}
for (const sourceContract of [
  "kSourceTargetOrder{{2, 0, 1}}",
  "kTargetParentOnePair{{0, 2, 1}}",
  "kTargetParentTwoPair{{1, 0, 2}}",
  "build_target_walk",
  "workspace.coordinates.size() >= 4",
  "source_siegmund_discrete",
  "exact_transition_budget_ok",
  "exact_excursion_probability",
  "maximum_exact_state_transitions",
  "exact_probability_cache.size() >= 8192",
  "source_scale",
  "result.probability = 1.0e-300",
  "source_excursion",
  "Literal CheckwrapC prefix extension",
  "result.probability_magnitude = result.magnitude",
  "source_sub_probability",
  "std::upper_bound",
  "source_check_split",
  "post_erasure_split_enabled",
  "alternative_split",
  "missing_data_split_applied = split.missing_data_found",
  "source_corrected_probability",
  "source_findall_corrected_probability",
  "raw > 1.0e-15",
  "source_threshold_passes",
  "std::expm1",
  "std::log1p",
  "selected_minus - selected_plus == selected.magnitude",
  "probability.probability > 0.0",
  "threeseq_discover_prepared",
  "threeseq_recheck_prepared",
  "threeseq_plot_profile",
]) {
  if (!threeSeq.includes(sourceContract)) {
    fail(`3SEQ source contract is missing ${sourceContract}`);
  }
}
if (!threeSeqHeader.includes("std::vector<float> probability_state") ||
    !threeSeqHeader.includes("std::vector<float> probability_next")) {
  fail("3SEQ exact DP does not retain supplied Single state precision");
}
if (!method.includes("threeseq_discover_prepared(") ||
    !method.includes("threeseq_recheck_prepared(") ||
    !method.includes("SignalMethod::threeseq") ||
    !method.includes("maxchi_workspace_") ||
    !method.includes("maxchi_workspace_.triplet_missing_data") ||
    !method.includes("post_erasure_split_enabled = !events_.empty()") ||
    !method.includes("profile.similarities")) {
  fail("3SEQ does not reuse fused triplet preparation or enter the cyclic signal path");
}
const geneconvDispatch = method.indexOf("if (options_.geneconv_enabled)", method.indexOf("void RdpScanner::scan_triplet"));
const maxChiDispatch = method.indexOf("if (options_.maxchi_enabled)", method.indexOf("void RdpScanner::scan_triplet"));
const chimaeraDispatch = method.indexOf("if (options_.chimaera_enabled)", method.indexOf("void RdpScanner::scan_triplet"));
const threeSeqDispatch = method.indexOf("if (options_.threeseq_enabled)", method.indexOf("void RdpScanner::scan_triplet"));
if (!(geneconvDispatch >= 0 && geneconvDispatch < maxChiDispatch &&
      maxChiDispatch < chimaeraDispatch && chimaeraDispatch < threeSeqDispatch)) {
  fail("combined scan does not retain the supplied GENECONV/MaxChi/CHIMAERA/3SEQ dispatch order");
}
for (const evidenceContract of [
  "threeSeqDiscovery",
  "FindSubSeqTS-Seq3PVals-CheckwrapC-TSXOver",
  "target-specific-information-rich-random-walk",
  "exact-hypergeometric-walk-with-siegmund-fallback",
  "dunn-sidak-when-project-correction-enabled",
  "targetLocal",
  "walkDirection",
  "informationRichSites",
  "parentOneMatches",
  "parentTwoMatches",
  "probabilityExcursion",
  "maximumExcursion",
  "exactProbability",
  "siegmundFallback",
  "missingDataSplitApplied",
]) {
  if (!method.includes(`\\\"${evidenceContract}\\\"`) || !types.includes(evidenceContract)) {
    fail(`3SEQ evidence contract differs between core and web for ${evidenceContract}`);
  }
}
if (!header.includes("int threeseq_enabled") ||
    !implementation.includes("int threeseq_enabled")) {
  fail("3SEQ scan ABI option is missing");
}
for (const restoreContract of [
  "information_rich_sites < 4",
  "maximum_excursion > parent_two_matches",
  "probability_excursion > information_rich_sites",
  "missing_data_split_applied == 0",
  "probability_excursion > maximum_excursion",
  "signal.informative_sites != information_rich_sites",
  "signal.local_p_value != raw_p_value",
  "missing_data_split_applied != 0",
  "discovery.major_parent_local == expected_major",
  "discovery.minor_parent_local == expected_minor",
  "discovery.candidate_pair == expected_pair",
]) {
  if (!implementation.includes(restoreContract)) {
    fail(`3SEQ restore validation is missing ${restoreContract}`);
  }
}
for (const workerContract of [
  "threeSeqEnabled: number",
  "_rdp_restore_threeseq_discovery",
  'savedMethod === "3SEQ"',
  "supportsThreeSeqDiscovery",
  "supportsThreeSeqSplit",
  "analysis.threeSeqProfilesScanned",
  "analysis.threeSeqExactEvaluations",
  "analysis.threeSeqApproximateEvaluations",
  "analysis.threeSeqCandidatesFound",
]) {
  if (!worker.includes(workerContract)) {
    fail(`3SEQ worker/restore contract is missing ${workerContract}`);
  }
}
for (const resultContract of [
  "threeSeqProfilesScanned",
  "threeSeqExactEvaluations",
  "threeSeqApproximateEvaluations",
  "threeSeqCandidatesFound",
  "threeSeqEnabled",
]) {
  if (!method.includes(resultContract) || !types.includes(`${resultContract}:`)) {
    fail(`3SEQ result contract differs between core and web for ${resultContract}`);
  }
  if (!app.includes(resultContract) || !scan.includes(resultContract)) {
    fail(`3SEQ application progress/options surface is missing ${resultContract}`);
  }
}
if (!app.includes("threeSeqEnabled: true") ||
    !settings.includes("Discover events with 3SEQ") ||
    !settings.includes("Siegmund fallback") ||
    !settings.includes('name: "3SEQ"') ||
    !settings.includes('state: "ready"')) {
  fail("3SEQ discovery defaults/settings are incomplete");
}
for (const uiContract of [
  "3SEQ hypergeometric random-walk excursion",
  "Corrected 3SEQ hit",
  "Dunn–Šidák corrected",
  "FindSubSeqTS → Seq3PVals/GetTSPVal → CheckwrapC → TSXOver",
  "Target-specific hypergeometric random walks",
]) {
  if (!review.includes(uiContract)) {
    fail(`3SEQ review workflow is missing ${uiContract}`);
  }
}
for (const recheckContract of [
  "threeSeqTripletRecheck",
  "postGroupThreeSeqRecheck",
  "TSXOver-Findall-two-orientations",
  "3SEQ triplet recheck status",
  "/3SEQ-recheck=",
  "/3SEQ-source-list-entries=",
]) {
  if (!method.includes(recheckContract)) {
    fail(`3SEQ representative/finalized-list recheck is missing ${recheckContract}`);
  }
}
for (const recheckContract of [
  "interface ThreeSeqRecheckEvidence",
  "threeSeqTripletRecheck: ThreeSeqRecheckEvidence",
  "postGroupThreeSeqRecheck: ThreeSeqRecheckEvidence",
  "qualifyingOrientations",
  "sourceListEntries",
  "bestDirection",
]) {
  if (!types.includes(recheckContract)) {
    fail(`3SEQ web recheck contract is missing ${recheckContract}`);
  }
}
for (const recheckContract of [
  "3SEQ Findall recheck",
  "TSXOver(1) Findall recheck",
  "postGroupThreeSeqRecheck",
  "five late rechecks",
  "inverse-interval list copy",
]) {
  if (!review.includes(recheckContract)) {
    fail(`3SEQ review recheck surface is missing ${recheckContract}`);
  }
}
if (!signalPlot.includes('signal.method === "3SEQ"') ||
    !signalPlot.includes("3SEQ target-specific hypergeometric random walks") ||
    !signalPlot.includes("Candidate target") ||
    !signalPlot.includes("plot.minimumValue")) {
  fail("3SEQ three-target random-walk plot contract is incomplete");
}
for (const traceContract of [
  "FindSubSeqDP3/6",
  "FastRecCheckChim",
  "AlistChi",
  "CXoverA",
  "No alternate RDP implementation was consulted",
]) {
  if (!chimaeraTrace.includes(traceContract)) {
    fail(`CHIMAERA supplied-source trace is missing ${traceContract}`);
  }
}
for (const traceContract of [
  "TSXOver(1)",
  "both Findall orientations",
  "inverse-parent/inverse-interval",
  "five",
]) {
  if (!lateTrace.includes(traceContract)) {
    fail(`late-consensus supplied-source trace is missing ${traceContract}`);
  }
}
for (const traceContract of [
  "FindSubSeqTS",
  "FindSubSeqTS2",
  "CheckwrapC",
  "Seq3PVals",
  "Get3SeqPvalC",
  "GetTSPVal",
  "SiegmundDiscrete",
  "TSXOver",
  "Dunn–Šidák",
  "No alternate RDP or 3SEQ implementation was consulted",
  "source-shaped and active",
  "CheckSplit3Seq",
  "native saved-output fixtures",
]) {
  if (!threeSeqTrace.includes(traceContract)) {
    fail(`3SEQ supplied-source trace is missing ${traceContract}`);
  }
}
for (const traceContract of [
  "FindSubSeqGCAP6",
  "GCXoverDP2",
  "GetFragsP",
  "GetMaxFragScoreP",
  "CalcKMaxP",
  "GCCalcPValP2",
  "Module31.bas::GCXoverD",
  "No alternate RDP or GENECONV implementation was\nconsulted",
  "source-shaped and active",
  "fragment-count-bound indexing quirk",
  "bounded bisection fallback",
]) {
  if (!geneconvTrace.includes(traceContract)) {
    fail(`GENECONV supplied-source trace is missing ${traceContract}`);
  }
}
for (const queryReferenceContract of [
  "AnalysisMode::query_reference",
  "query_sequences_",
  "reference_sequences_",
  "reference_pair_is_valid",
  "advance_reference_pair",
  "one-query-two-different-reference-groups",
  "reference-group-pairs-times-query-origins",
  "queryReferenceInputRole",
  "referenceGroup",
  "saved query-vs-reference signal violates its one-query",
  "referenceGroupIndices",
  "query-reference",
]) {
  if (!method.includes(queryReferenceContract) &&
      !types.includes(queryReferenceContract) &&
      !worker.includes(queryReferenceContract) &&
      !settings.includes(queryReferenceContract)) {
    fail(`query-vs-reference workflow is missing ${queryReferenceContract}`);
  }
}
for (const [label, headerContract, implementationContract] of [
  ["mode", "int query_reference_mode", "int query_reference_mode"],
  ["groups", "const uint32_t* reference_groups", "const std::uint32_t* reference_groups"],
  ["group count", "size_t reference_group_count", "std::size_t reference_group_count"],
]) {
  if (!header.includes(headerContract) ||
      !implementation.includes(implementationContract)) {
    fail(`query-vs-reference ABI is missing ${label}`);
  }
}
if (!datasetStep.includes("Detect REF names") ||
    !settings.includes("Query vs reference") ||
    !app.includes("inferReferenceGroups") ||
    !worker.includes("referenceGroupsPointer")) {
  fail("query-vs-reference configuration is not wired across dataset, settings, and worker");
}
for (const bulkReferenceContract of [
  "selectedSequenceIndices",
  "Select every sequence matching the current filter",
  "Assign group",
  "Make queries",
  "matchingSequences.length > 500",
]) {
  if (!datasetStep.includes(bulkReferenceContract)) {
    fail(`bulk query/reference assignment is missing ${bulkReferenceContract}`);
  }
}
if (!app.includes("compactReferenceGroups") ||
    !datasetStep.includes("Compact groups") ||
    !datasetStep.includes("eligibleSequenceCount")) {
  fail("reference-group compaction or exact dataset eligibility feedback is missing");
}
if (!types.includes("correctionTests: number") ||
    !app.includes("queryReferenceCorrectionTestCount") ||
    !app.includes("Math.floor((255 ** 4) / 2)") ||
    !settings.includes("group-pair × query opportunities after the native cap")) {
  fail("query-vs-reference correction count is not exposed in live workflow state");
}
for (const progressRoleContract of [
  "activeWorkingSequenceCount",
  "queryWorkingSequenceCount",
  "referenceWorkingSequenceCount",
  "activeReferenceGroupCount",
]) {
  if (!method.includes(`\\\"${progressRoleContract}\\\"`) ||
      !types.includes(`${progressRoleContract}: number`) ||
      !app.includes(progressRoleContract) ||
      !scan.includes(progressRoleContract)) {
    fail(`live cyclic role progress is missing ${progressRoleContract}`);
  }
}
for (const shortlistProgressContract of [
  "tripletKernelEvaluations",
  "tripletSummariesReused",
  "cachedSignalsReused",
  "methodScansSkipped",
]) {
  if (!method.includes(`\\\"${shortlistProgressContract}\\\"`) ||
      !types.includes(`${shortlistProgressContract}: number`) ||
      !app.includes(shortlistProgressContract) ||
      !scan.includes(shortlistProgressContract)) {
    fail(`cyclic-shortlist progress is missing ${shortlistProgressContract}`);
  }
}
for (const shortlistCoreContract of [
  "working_triplet_available",
  "round_triplet_signal_summaries_",
  "carried_triplet_signal_summaries_",
  "dirty_working_sequences_",
  "reuse_carried_triplet_signals",
  "triplet_touches_dirty_sequence",
  "refresh_threeseq_on_unchanged_triplets_",
  "correction_tests_frozen_",
  "first_post_erasure_threeseq_refresh",
]) {
  if (!method.includes(shortlistCoreContract)) {
    fail(`cyclic-shortlist core is missing ${shortlistCoreContract}`);
  }
}
if (!pagesVerifier.includes("progress.correctionTests !== 120") ||
    !pagesVerifier.includes("progress.tripletKernelEvaluations < progress.cumulativeTriplets") ||
    !pagesVerifier.includes("progress.cachedSignalsReused > 0")) {
  fail("Pages verification does not exercise fixed correction and shortlist reuse");
}
if (!review.includes("RDP5 XOverList equivalent") ||
    !review.includes("BURT can move the displayed event") ||
    !review.includes("initial scan-plan opportunities")) {
  fail("GENECONV probability scope and BURT ordering are not explicit in review");
}
if (!worker.includes("supportsReferenceGroups && analysis.analysisMode") ||
    !app.includes("restored.results.analysisMode")) {
  fail("query-vs-reference projects do not preserve their schema-gated scan plan");
}
if (!review.includes("Reference recombinant · group") ||
    !review.includes("allowed and explicitly documented query-vs-reference outcome")) {
  fail("review does not distinguish a reference sequence called as recombinant");
}
if (!method.includes("signal_reference_group") ||
    !types.includes("interface RdpSignal") ||
    !types.includes("queryReferenceInputRole: QueryReferenceInputRole")) {
  fail("per-signal query/reference input-role metadata is incomplete");
}
if (!eventAlignmentInspector.includes("reference group") ||
    !eventTreeInspector.includes("leaf.queryReferenceInputRole")) {
  fail("query-vs-reference input roles are missing from graphical event inspectors");
}
if (!method.includes("Active reference assignments") ||
    !method.includes("Recombinant input role") ||
    !exportStep.includes("reference-recombinant call")) {
  fail("query-vs-reference event exports do not retain visible input-role context");
}
for (const handoffContract of [
  "0.15.0-session-15",
  "org.rdp-web.project/v1alpha15",
  "FindSubSeqTS",
  "CheckwrapC",
  "exact hypergeometric random-walk",
  "SiegmundDiscrete",
  "Dunn–Šidák",
  "CheckSplit3Seq",
  "TSXOver(1)",
  "Findall",
  "inverse-parent/inverse-interval",
  "no supplied/native source, C++, WebAssembly, TypeScript compiler",
]) {
  if (!sessionHandoff.includes(handoffContract)) {
    fail(`Session 15 handoff is missing ${handoffContract}`);
  }
}
for (const handoffContract of [
  "0.17.0-session-17",
  "org.rdp-web.project/v1alpha16",
  "Raw KA P",
  "Project corrected",
  "GCCalcPValP2",
  "MakeMCCorrection",
  "BestXOList",
  "Worthwhilescan",
  "triplet-kernel evaluations",
  "No alternate RDP implementation was consulted",
]) {
  if (!session17Handoff.includes(handoffContract)) {
    fail(`Session 17 handoff is missing ${handoffContract}`);
  }
}
for (const shortlistContract of [
  "XOverList",
  "XOverDefine",
  "BestXOList",
  "Worthwhilescan",
  "StoreLPV",
  "FindBetterRecSignal",
  "correctionTests",
  "tripletKernelEvaluations",
  "CheckSplit3Seq",
  "No alternate RDP implementation",
]) {
  if (!shortlistTrace.includes(shortlistContract)) {
    fail(`cyclic-shortlist supplied-source trace is missing ${shortlistContract}`);
  }
}
for (const handoffContract of [
  "0.18.0-session-18",
  "org.rdp-web.project/v1alpha16",
  "Windows 95",
  "user-stopped",
  "round_signal_begin_",
  "graceful-stop regression",
]) {
  if (!session18Handoff.includes(handoffContract)) {
    fail(`Session 18 handoff is missing ${handoffContract}`);
  }
}
if (!readme.includes("Session 18 source checkpoint") ||
    !readme.includes("v1alpha16") ||
    !readme.includes("cyclic-shortlist trace") ||
    !readme.includes("session-18-handoff") ||
    !status.includes("Port status — session 18") ||
    !status.includes("3SEQ exploratory discovery") ||
    !status.includes("3SEQ Findall recheck") ||
    !status.includes("XOverList/BestXOList-style shortlist")) {
  fail("Session 18 README/status documentation is stale");
}
if (!exportStep.includes("Session 18 snapshot") ||
    !exportStep.includes("target-rotated 3SEQ") ||
    !exportStep.includes("CheckSplit3Seq") ||
    !exportStep.includes("TSXOver(1)")) {
  fail("Session 18 export fidelity boundary is stale");
}
for (const win95Contract of [
  "Windows 95 visual skin",
  "--win-face: #c0c0c0",
  "--win-blue: #000080",
  "window-titlebar",
  "window-menubar",
  "app-statusbar",
  "repeating-linear-gradient",
]) {
  if (!styles.includes(win95Contract) && !app.includes(win95Contract)) {
    fail(`Windows 95 visual contract is missing ${win95Contract}`);
  }
}
for (const stopContract of [
  'cycle_termination_ = "user-stopped"',
  "signals_.resize(round_signal_begin_)",
  "cancelled_.store(false)",
]) {
  if (!method.includes(stopContract)) {
    fail(`graceful cyclic stop is missing ${stopContract}`);
  }
}
if (!scan.includes("Stop and review completed events") ||
    !scan.includes("The unfinished cyclic round was discarded") ||
    !pagesVerifier.includes("graceful-stop results") ||
    !pagesVerifier.includes('cycleTermination !== "user-stopped"')) {
  fail("graceful cyclic stop UI/runtime regression coverage is incomplete");
}
if (!method.includes("std::numeric_limits<std::uint64_t>::max() / reference_group_pairs")) {
  fail("query-vs-reference correction multiplication is not overflow guarded");
}
for (const maxChiOptionContract of [
  "int maxchi_enabled",
  "maxchi_window_sites",
]) {
  if (!header.includes(maxChiOptionContract) ||
      !implementation.includes(maxChiOptionContract)) {
    fail(`MaxChi scan ABI option is missing ${maxChiOptionContract}`);
  }
}
for (const maxChiWorkerContract of [
  "maxChiEnabled: number",
  "maxChiWindowSites: number",
  "_rdp_restore_maxchi_discovery",
  'savedMethod === "MAXCHI"',
  "analysis.maxChiProfilesScanned",
  "analysis.maxChiPeakAttempts",
  "analysis.maxChiCandidatesFound",
  "analysis.maxChiPeakLimitTriplets",
]) {
  if (!worker.includes(maxChiWorkerContract)) {
    fail(`MaxChi worker/restore contract is missing ${maxChiWorkerContract}`);
  }
}
for (const restoreCounterContract of [
  "double cumulative_triplets",
  "double maxchi_profiles_scanned",
  "double maxchi_peak_attempts",
  "double maxchi_candidates_found",
  "double maxchi_peak_limit_triplets",
  "double chimaera_profiles_scanned",
  "double chimaera_peak_attempts",
  "double chimaera_candidates_found",
  "double chimaera_peak_limit_targets",
  "double geneconv_fragments_scored",
  "double geneconv_qualified_fragments",
  "double geneconv_candidates_found",
  "double geneconv_overlap_rejections",
  "double geneconv_numerical_fallback_tracks",
  "double threeseq_profiles_scanned",
  "double threeseq_exact_evaluations",
  "double threeseq_approximate_evaluations",
  "double threeseq_candidates_found",
]) {
  if (!header.includes(restoreCounterContract) || !implementation.includes(restoreCounterContract)) {
    fail(`completed-project scan counter restore is missing ${restoreCounterContract}`);
  }
}
for (const restoreSummaryContract of [
  "cycle_termination_length",
  "cumulative_triplets_authoritative_",
  "processedTriplets",
  "totalTriplets",
]) {
  if (!method.includes(restoreSummaryContract) &&
      !header.includes(restoreSummaryContract) &&
      !implementation.includes(restoreSummaryContract)) {
    fail(`completed-project summary restore is missing ${restoreSummaryContract}`);
  }
}
if (!worker.includes("cycleTerminationBytes") ||
    !app.includes("restored.results.processedTriplets") ||
    !app.includes("restored.results.totalTriplets")) {
  fail("completed-project stop reason or final-round progress is not restored exactly");
}
for (const plotContract of [
  "plot_sample_indices",
  "peak_alignment_position",
  "detectionProfileExact",
  "original-alignment-reconstruction",
  "random-walk-height",
  "minimumValue",
]) {
  if (!method.includes(plotContract) &&
      !types.includes(plotContract) &&
      !signalPlot.includes(plotContract)) {
    fail(`method-aware plot contract is missing ${plotContract}`);
  }
}
if (!signalPlot.includes("Original-alignment reconstruction")) {
  fail("later-round plot reconstruction is not disclosed in review UI");
}
for (const maxChiResultContract of [
  "discoveryMethods",
  "maxChiProfilesScanned",
  "maxChiPeakAttempts",
  "maxChiCandidatesFound",
  "maxChiPeakLimitTriplets",
  "anchorMethod",
  "detectionMethods",
  "chi-square",
]) {
  if (!method.includes(maxChiResultContract) || !types.includes(maxChiResultContract)) {
    fail(`MaxChi result contract differs between core and web for ${maxChiResultContract}`);
  }
}
if (!app.includes("maxChiEnabled: true") ||
    !app.includes("maxChiWindowSites: 70")) {
  fail("MaxChi discovery defaults are not enabled in the application state");
}
for (const exportContract of [
  ["rdp_export_enabled_sequences_fasta", implementation, worker],
  ["rdp_export_masked_or_disabled_sequences_fasta", implementation, worker],
  ["export-enabled-sequences", worker, types, client],
  ["export-masked-or-disabled-sequences", worker, types, client],
  ["onEnabledSequences", app, exportStep],
  ["onMaskedOrDisabledSequences", app, exportStep],
  ["onExportFullAlignment", app, datasetStep],
  ["onFullAlignment", app, exportStep],
  ["rdp_export_recombinant_sequences_removed_fasta", implementation, worker],
  ["rdp_export_recombinant_columns_removed_fasta", implementation, worker],
  ["export-recombinant-sequences-removed", worker, types, client],
  ["export-recombinant-columns-removed", worker, types, client],
  ["onRecombinantSequencesRemoved", app, exportStep],
  ["onRecombinantColumnsRemoved", app, exportStep],
]) {
  const [name, ...sources] = exportContract;
  if (sources.some((source) => !source.includes(name))) {
    fail(`accepted-event FASTA export contract is missing ${name}`);
  }
}
for (const curationExportContract of [
  "curated_sequences_fasta",
  "exportCuratedSequences",
  "maskedSequenceIndices",
  "disabledSequenceIndices",
]) {
  if (!method.includes(curationExportContract) &&
      !implementation.includes(curationExportContract) &&
      !worker.includes(curationExportContract)) {
    fail(`pre-scan curation export contract is missing ${curationExportContract}`);
  }
}
for (const erasureContract of [
  "if (beginning == ending) return true",
  "coordinate >= beginning || coordinate <= ending",
  "coordinate >= beginning && coordinate <= ending",
  "breakpoint_polish_missing_scratch_",
  "prior_index < event.id && prior_index < events_.size()",
  "triplet_missing[coordinate] = 1",
]) {
  if (!method.includes(erasureContract)) {
    fail(`inclusive erasure/BURT missing-data contract is missing ${erasureContract}`);
  }
}
if (method.includes("coordinate > beginning || coordinate < ending") ||
    method.includes("coordinate > beginning && coordinate < ending")) {
  fail("the shared active tract helper regressed to endpoint-exclusive erasure");
}
if (!method.includes("beginning != ending && region_length > 2")) {
  fail("primary RDP detection no longer preserves FastRecCheckP's EN != BE gate");
}
for (const sourceContract of [
  "kHmmCyclesArgument = 20",
  "kHmmTrainingStarts = kHmmCyclesArgument + 1",
  "kMaximumTrainingIterations = 100",
  "kPseudocount = 0.01",
  "kSourceDefaultSeed = 3",
  "class SourceMsvcRandom",
  "0.995",
  "0.999",
  "source_match_breakpoint",
  "source_missing_data_reposition",
  "triplet_missing_data",
  "information_count == 0",
  "workspace.lattice.assign(training_stride * kStates, 0.0)",
  "insufficient_inside_or_outside_reverted",
]) {
  if (!burt.includes(sourceContract)) {
    fail(`BURT/BenHMM source contract is missing ${sourceContract}`);
  }
}

console.log(
  `Source contracts verified: ${headerFunctions.length} RDP ABI functions, engine ${version}, schema ${schema}.`,
);
