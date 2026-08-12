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
  worker,
  types,
  client,
  app,
  datasetStep,
  settings,
  review,
  exportStep,
  pagesWorkflow,
  viteConfig,
  pagesVerifier,
] = await Promise.all([
  read("package.json"),
  read("package-lock.json"),
  read("wasm/CMakeLists.txt"),
  read("wasm/include/rdp_api.h"),
  read("wasm/src/rdp_api.cpp"),
  read("wasm/src/rdp_method.cpp"),
  read("wasm/src/burt_confidence.cpp"),
  read("wasm/src/maxchi.cpp"),
  read("src/workers/analysis.worker.ts"),
  read("src/lib/types.ts"),
  read("src/lib/wasmClient.ts"),
  read("src/App.tsx"),
  read("src/components/DatasetStep.tsx"),
  read("src/components/SettingsStep.tsx"),
  read("src/components/ReviewStep.tsx"),
  read("src/components/ExportStep.tsx"),
  read(".github/workflows/deploy-pages.yml"),
  read("vite.config.ts"),
  read("scripts/verify-pages-output.mjs"),
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

const cmakeExports = uniqueSorted(
  captures(cmake, /'(_(?:rdp_[a-z0-9_]+|malloc|free))'/g),
);
const expectedCmakeExports = uniqueSorted([
  "_malloc",
  "_free",
  ...headerFunctions.map((name) => `_${name}`),
]);
compareSets("CMake exported function list differs from rdp_api.h", expectedCmakeExports, cmakeExports);

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

const schema = "org.rdp-web.project/v1alpha9";
if (!implementation.includes(schema) || !worker.includes(schema)) {
  fail(`emitted/imported project schema ${schema} is not aligned`);
}
if (!worker.includes("assetVersion") || !worker.includes("loadedVersion !== assetVersion")) {
  fail("worker asset/engine version guard is missing");
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
  "npm run check:source",
  "npm run check:types",
  "npm run build",
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
  "is a symbolic link",
  "is hard-linked",
]) {
  if (!pagesVerifier.includes(artifactContract)) {
    fail(`Pages artifact verifier is missing ${artifactContract}`);
  }
}

const activeLateConsensusStatus = "active-rdp-maxchi-post-group-recheck";
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
  "nativeMaxChiFullRecheckComplete",
]) {
  if (!method.includes(statusContract) || !types.includes(statusContract)) {
    fail(`MaxChi late-consensus status is missing ${statusContract}`);
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
  fail("MaxChi recheck implementation is not part of the WASM target");
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
]) {
  if (!maxchi.includes(sourceContract)) {
    fail(`MaxChi source contract is missing ${sourceContract}`);
  }
}
for (const uiContract of [
  "maxChiTripletRecheck",
  "postGroupMaxChiRecheck",
  "MaxChi recheck",
  "This is triplet corroboration, not MaxChi event discovery.",
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
]) {
  if (!method.includes(`\\\"${evidenceContract}\\\"`) || !types.includes(evidenceContract)) {
    fail(`MaxChi evidence contract differs between core and web for ${evidenceContract}`);
  }
}
if (!settings.includes("RDP full · MaxChi recheck") ||
    !settings.includes("Included in confirmation")) {
  fail("method settings do not distinguish MaxChi confirmation from RDP discovery");
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
