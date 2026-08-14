export type WorkflowStep = "dataset" | "settings" | "scan" | "review" | "export";

export type CorrectionMode = "bonferroni" | "none";
export type AnalysisMode = "exploratory" | "query-reference";
export type QueryReferenceInputRole = "query" | "reference" | "not-applied";
export type SequenceAnalysisState = "enabled" | "masked" | "disabled";
export type PhylproGapMode = "ignore-missing-pairwise" | "strip-any-missing-column";
export type DiscoveryMethod = "RDP" | "MAXCHI" | "CHIMAERA" | "GENECONV" | "3SEQ" | "BOOTSCAN" | "SISCAN";

export interface SequenceSummary {
  index: number;
  name: string;
  validSites: number;
  missingSites: number;
  missingFraction: number;
  masked: boolean;
}

export interface DatasetSummary {
  format: string;
  sequenceCount: number;
  alignmentLength: number;
  activeSequenceCount: number;
  tripletCount: number;
  variableSiteCount: number;
  informativeSiteCount: number;
  minimumPairIdentity: number | null;
  meanPairIdentity: number | null;
  recommendedMinimumDistance: number;
  partitionBoundaries: number[];
  sequences: SequenceSummary[];
  warnings: string[];
}

export interface ScanOptions {
  cpuThreads: number;
  analysisMode: AnalysisMode;
  circular: boolean;
  pValueCutoff: number;
  correction: CorrectionMode;
  windowSites: number;
  maxChiEnabled: boolean;
  maxChiWindowSites: number;
  chimaeraEnabled: boolean;
  chimaeraWindowSites: number;
  geneconvEnabled: boolean;
  geneconvMismatchScale: number;
  geneconvMaxOverlaps: number;
  threeSeqEnabled: boolean;
  bootscanPrimaryEnabled: boolean;
  bootscanSecondaryEnabled: boolean;
  bootscanWindowSites: number;
  bootscanStepSites: number;
  bootscanBootstrapReplicates: number;
  bootscanSupportCutoff: number;
  bootscanRandomSeed: number;
  siscanPrimaryEnabled: boolean;
  siscanSecondaryEnabled: boolean;
  siscanWindowSites: number;
  siscanStepSites: number;
  siscanScanPermutations: number;
  siscanPValuePermutations: number;
  siscanRandomSeed: number;
  polishBreakpoints: boolean;
  maskedSequenceIndices: number[];
  disabledSequenceIndices: number[];
  referenceGroupIndices: number[];
}

export interface ScanRoundTiming {
  round: number;
  elapsedMs: number;
  completed: boolean;
}

export interface ScanTiming {
  startedAt: string;
  totalMs: number;
  setupMs: number;
  primaryMs: number;
  cyclicRescanMs: number;
  reconciliationMs: number;
  currentRoundMs: number;
  rounds: ScanRoundTiming[];
}

export interface ScanExecution {
  mode: "single-worker" | "wasm-pthreads";
  requestedThreads: number;
  activeThreads: number;
  hardwareConcurrency: number;
}

export interface EngineRuntimeInfo {
  threaded: boolean;
  version: string;
  hardwareConcurrency: number;
  maximumThreads: number;
  recommendedThreads: number;
}

export interface ScanProgress {
  state: "idle" | "running" | "done" | "cancelled" | "error";
  phase: "primary" | "cyclic-rescan" | "reconciliation" | "complete";
  processedTriplets: number;
  totalTriplets: number;
  correctionTests: number;
  activeWorkingSequenceCount: number;
  queryWorkingSequenceCount: number;
  referenceWorkingSequenceCount: number;
  activeReferenceGroupCount: number;
  cumulativeTriplets: number;
  tripletKernelEvaluations: number;
  tripletSummariesReused: number;
  cleanTripletsPruned: number;
  cachedSignalsReused: number;
  methodScansSkipped: number;
  invalidScheduleTripletsSkipped: number;
  fragmentSequencesPruned: number;
  scanRound: number;
  fixedEventCount: number;
  signalCount: number;
  eventCount: number;
  maxChiProfilesScanned: number;
  maxChiPeakAttempts: number;
  maxChiCandidatesFound: number;
  maxChiPeakLimitTriplets: number;
  chimaeraProfilesScanned: number;
  chimaeraPeakAttempts: number;
  chimaeraCandidatesFound: number;
  chimaeraPeakLimitTargets: number;
  geneconvFragmentsScored: number;
  geneconvQualifiedFragments: number;
  geneconvCandidatesFound: number;
  geneconvOverlapRejections: number;
  geneconvNumericalFallbackTracks: number;
  threeSeqProfilesScanned: number;
  threeSeqExactEvaluations: number;
  threeSeqApproximateEvaluations: number;
  threeSeqCandidatesFound: number;
  bootscanProfilesScanned: number;
  bootscanCandidateRegionsScored: number;
  bootscanCandidatesFound: number;
  bootscanPairProfilesRequested: number;
  bootscanPairProfileCacheHits: number;
  bootscanPairProfileCacheMisses: number;
  bootscanPairProfileCacheEvictions: number;
  bootscanPairProfileCacheBytes: number;
  bootscanPairProfileCachePeakBytes: number;
  siscanProfilesScanned: number;
  siscanWindowsScored: number;
  siscanCandidateRegionsScored: number;
  siscanCandidatesFound: number;
  siscanPermutationDraws: number;
  siscanContextBuilds: number;
  siscanContextPairComparisons: number;
  siscanContextTreeMerges: number;
  siscanRandomValuesGenerated: number;
  cycleTermination: string;
  fraction: number;
  timing: ScanTiming | null;
  execution: ScanExecution;
}

export type ReviewState = "unreviewed" | "accepted" | "rejected";

export interface MaxChiDiscoveryEvidence {
  status: "source-shaped-active-unvalidated";
  kernel: "MCXoverF-multi-peak-destroy-retry";
  peakOrdering: "raw-chi-square-lazy-heap";
  smoothingUse: "twelve-term-eleven-divisor-source-basin-destruction-only";
  peakAttempt: number;
  peakPair: 0 | 1 | 2 | null;
  tractSide: "left" | "right" | "unavailable";
  peakAlignmentPosition: number;
  variableSites: number;
  initialHalfWindow: number;
  grownHalfWindow: number;
  criticalDifference: number;
  maximumChiSquare: number;
  rawPValue: number;
  withinTripletPValue: number;
  correctedPValue: number;
  leftFlankChiSquare: number;
  rightFlankChiSquare: number;
  missingDataWindowFilterApplied: boolean;
  linearEdgeWindowFilterApplied: boolean;
}

export interface ChimaeraDiscoveryEvidence {
  status: "source-shaped-active-unvalidated";
  kernel: "AlistChi-FastRecCheckChim-CXoverA";
  profile: "target-specific-information-rich-binary-string";
  peakOrdering: "raw-chi-square-lazy-heap-per-target";
  smoothingUse: "twelve-term-eleven-divisor-source-basin-destruction-only";
  targetLocal: 0 | 1 | 2;
  peakAttempt: number;
  tractSide: "left" | "right" | "unavailable";
  peakAlignmentPosition: number;
  informationRichSites: number;
  initialHalfWindow: number;
  grownHalfWindow: number;
  criticalDifference: number;
  maximumChiSquare: number;
  rawPValue: number;
  withinTripletPValue: number;
  correctedPValue: number;
  leftFlankChiSquare: number;
  rightFlankChiSquare: number;
  insideParentOneMatchRate: number;
  outsideParentOneMatchRate: number;
  missingDataWindowFilterApplied: boolean;
  linearEdgeWindowFilterApplied: boolean;
}

export interface GeneconvDiscoveryEvidence {
  status: "source-shaped-active-unvalidated";
  kernel: "FindSubSeqGCAP6-GetFragsP-GetMaxFragScoreP-CalcKMaxP-GCCalcPValP2-GCXoverD";
  probabilityModel: "karlin-altschul";
  indelMode: "ignored";
  overlapPolicy: "stable-lowest-p-configured-coverage";
  minimumFragmentFiltersApplied: false;
  track: 0 | 1 | 2 | 3 | 4 | 5;
  polymorphicSites: number;
  positiveSites: number;
  discordantSites: number;
  mismatchPenalty: number;
  fragmentScore: number;
  criticalScore: number;
  lambda: number;
  karlinAltschulK: number;
  rawPValue: number;
  correctedPValue: number;
  karlinAltschulProbability: true;
  ignoredIndels: true;
  overlapFilterApplied: true;
}

export interface ThreeSeqDiscoveryEvidence {
  status: "source-shaped-active-unvalidated";
  kernel: "FindSubSeqTS-Seq3PVals-CheckwrapC-TSXOver";
  profile: "target-specific-information-rich-random-walk";
  probabilityModel: "exact-hypergeometric-walk-with-siegmund-fallback";
  correctionModel: "dunn-sidak-when-project-correction-enabled";
  targetLocal: 0 | 1 | 2;
  walkDirection: "ascent" | "descent";
  informationRichSites: number;
  parentOneMatches: number;
  parentTwoMatches: number;
  probabilityExcursion: number;
  maximumExcursion: number;
  rawPValue: number;
  correctedPValue: number;
  exactProbability: boolean;
  siegmundFallback: boolean;
  missingDataSplitApplied: boolean;
}

export interface BootscanDiscoveryEvidence {
  status: "source-shaped-active-unvalidated";
  kernel: "BSXoverR-SEQBOOT2-FastBootDist-GetPltVal-ScanBSPlots-MakeBSEvent";
  mode: "jukes-cantor-distance";
  probabilityModel: "MakeScoresBS-binomial";
  strictClosestPairVoting: true;
  supportedPair: 0 | 1 | 2;
  windowsScored: number;
  usableWindows: number;
  informativeSites: number;
  tractInformativeSites: number;
  tractPairMatches: number;
  outsidePairMatches: number;
  maximumPairSupport: number;
  meanPairSupport: number;
  bootstrapPValue: number;
  rawPValue: number;
  correctedPValue: number;
  erasedWindowFilterApplied: boolean;
}

export interface SiscanDiscoveryEvidence {
  status: "source-shaped-active-unvalidated";
  kernel: "SSXoverC-GetSSOL-Get3Score-GetPScores2-DoPerms3-MakeZValue2-DoSums-FindMaxZ-ShrinkRegionC";
  outlierMode: "nearest-source-wpgma";
  permutationGenerator: "microsoft-crt-flat-prefix";
  gapMode: "strip";
  variablePatternMode: "one-two-three-variable";
  sourceFastWindowQuirk: true;
  globalPair: 0 | 1 | 2;
  candidatePair: 0 | 1 | 2;
  outlierSequence: number;
  windowsInRegion: number;
  informativeSites: number;
  permutationDraws: number;
  selectedScore: number;
  selectedScoreFamily: "partition" | "summed";
  maximumZ: number;
  normalTailPValue: number;
  regionLengthAdjustedPValue: number;
  windowAdjustedPValue: number;
  correctedPValue: number;
}

export interface RdpSignal {
  id: number;
  method: DiscoveryMethod;
  triplet: [number, number, number];
  tripletNames: [string, string, string];
  recombinant: number;
  recombinantName: string;
  queryReferenceInputRole: QueryReferenceInputRole;
  referenceGroup: number | null;
  majorParent: number;
  majorParentName: string;
  minorParent: number;
  minorParentName: string;
  beginning: number;
  ending: number;
  wrapsOrigin: boolean;
  informativeBeginning: number;
  informativeEnding: number;
  localPValue: number;
  correctedPValue: number;
  correctionTests: number;
  pairSimilarity: [number, number, number];
  informativeSites: number;
  candidatePair: number;
  maxChiDiscovery: MaxChiDiscoveryEvidence | null;
  chimaeraDiscovery: ChimaeraDiscoveryEvidence | null;
  geneconvDiscovery: GeneconvDiscoveryEvidence | null;
  threeSeqDiscovery: ThreeSeqDiscoveryEvidence | null;
  bootscanDiscovery: BootscanDiscoveryEvidence | null;
  siscanDiscovery: SiscanDiscoveryEvidence | null;
  fragmentAssisted: boolean;
  fragmentEventContext: [number | null, number | null, number | null];
  eventId: number | null;
  reviewState: ReviewState;
  provisionalRoles: true;
}

export interface TraceEvidence {
  sequenceIndex: number;
  sequenceName: string;
  beginning: number;
  ending: number;
  wrapsOrigin: boolean;
  localPValue: number;
  correctedPValue: number;
  significant: boolean;
}

export interface RoleCandidateIndices {
  recombinant: number[];
  majorParent: number[];
  minorParent: number[];
}

export interface BreakpointErasureBoundary {
  erasureAdjacent: boolean;
  erasureWithinRdpWindow: boolean;
  uncertainDueToErasure: boolean;
  nativeCheckEndsApplied: boolean;
  nativeCheckEndsWarning: boolean;
  informationProfileAvailable: boolean;
  inputMissingDataInCheckRange: boolean;
  linearEdgeWithinRdpWindow: boolean;
  nativeCheckRange: {
    beginning: number;
    ending: number;
    wrapsOrigin: boolean;
    coordinateCount: number;
  };
  rdpWindowInformativeSites: number;
  nearestErasureInformativeSites: number | null;
  uncertainPriorEventIds: number[];
  priorEventIds: number[];
}

export interface BreakpointContext {
  source: "cyclic-erasure-history";
  beginning: BreakpointErasureBoundary;
  ending: BreakpointErasureBoundary;
}

export interface BreakpointConfidenceBoundary {
  name: "beginning" | "ending";
  inputCoordinate: number;
  polishedCoordinate: number;
  intervalAvailable: boolean;
  sourceIntervalContainsInput: boolean;
  confidence99: {
    beginning: number;
    ending: number;
    wrapsOrigin: boolean;
  };
  hmmCoordinate: number;
  confidence95: {
    beginning: number;
    ending: number;
    wrapsOrigin: boolean;
  };
  repositioned: boolean;
  missingDataAdjusted: boolean;
  finalGapAdjusted: boolean;
}

export interface BreakpointConfidenceEvidence {
  status: "complete-active-unvalidated" | "unavailable";
  method: "BURT/BenHMM";
  attempted: boolean;
  available: boolean;
  appliedToEvent: boolean;
  informationRichSites: number;
  candidateIntervalCount: number;
  bestLogLikelihood: number;
  randomSeed: 3;
  randomAdapter: "msvc-rand-15-bit";
  sourceRandomAdapter: true;
  hmmCyclesArgument: 20;
  serialTrainingStarts: 21;
  posteriorThresholds: [0.995, 0.999];
  polishedBeginning: number;
  polishedEnding: number;
  singleTransitionAssignment: boolean;
  insufficientInsideOrOutsideReverted: boolean;
  unavailableReason:
    | "not-run"
    | "disabled"
    | "invalid-triplet-or-breakpoint"
    | "no-information-rich-sites"
    | "hmm-training-unavailable"
    | "no-hmm-state-transition"
    | "no-matched-breakpoint-interval"
    | null;
  boundaries: [BreakpointConfidenceBoundary, BreakpointConfidenceBoundary];
}

export interface FinalTrimMatrixEvidence {
  status: "complete-active-rff0";
  implementedOKSeqElements: [7, 8, 9, 10, 11, 12, 13, 14];
  inactiveZeroOKSeqElements: [10, 11];
  collapsedTreePositionScore: number;
  rawTreePositionScore: number;
  relativeDistanceScore: number;
  breakpointDistanceScores: [number, number];
  breakpointScoreAvailable: [boolean, boolean];
  detectedRegionDistanceScore: number;
  detectedRegionMatchDistances: [number, number, number, number, number, number, number];
  detectedEventMatch: boolean;
  detectedEventOverlap: number;
  detectedEventBeginning: number | null;
  detectedEventEnding: number | null;
  detectedEventSignalId: number | null;
  detectedRegionSaturated: boolean;
  sourceSequenceIndexQuirkApplied: true;
  activeConsensusMatrixScore: number;
  treeDistanceFallback: boolean;
  appliesToNonrepresentative: boolean;
}

export interface CalcMatchEvidence {
  status: "complete-active-rff0" | "unavailable";
  implementedOKSeqElements: [17, 18];
  standardGroupingThresholds: true;
  fragmentVariableSites: number;
  targetHalfWindow: number;
  smoothingHalfWindow: number;
  regionalMatchScore: number;
  rawBreakpointMatchClass: -1 | 0 | 1 | 2;
  breakpointMatchClass: -1 | 0 | 1 | 2;
  checkpointMatches: [number, number, number, number, number, number];
  breakpointsExist: [boolean, boolean];
  topologyFiltered: boolean;
  topologyDistanceFallback: boolean;
  unavailableReason?: "insufficient-variable-sites-or-source-fragment-bound";
}

export interface ConsensusScoreEvidence {
  status: "complete-active-rff0";
  implementedOKSeqElements: [0, 1, 2, 3, 4, 5, 6, 15];
  correctedCorrelationPValue: number;
  detectedEventOverlap: number;
  patternScore: number;
  detectableSetMember: boolean;
  initialRListMember: boolean;
  duplicateCleanedMember: boolean;
  nearestNonrecombinantMember: boolean;
  finalTrimFirstExpansionAdded: boolean;
  finalTrimSecondExpansionAdded: boolean;
  selectedRolePrunedOut: boolean;
  finalTrimMember: boolean;
  maximumDirectCorrelation: number;
  baseScoreBeforeFinalMembership: number;
  scoreAfterFinalMembership: number;
  scoreAfterRcorrx: number;
  sourceLongMatrixMultiplier: number;
  sourceLongNsSemantics: true;
  finalScore: number;
  consensusPrimaryMember: boolean;
  consensusEquivalentMember: boolean;
  consensusStragglerMember: boolean;
  consensusRebuiltMember: boolean;
  consensusFallbackRestored: boolean;
  selectedTreeCleanupPrunedOut: boolean;
  selectedTreeCleanupAdded: boolean;
  finalDistanceMember: boolean;
  representativeSentinel: boolean;
  otherRepresentativeZero: boolean;
  complete: true;
}

export interface PostGroupRdpRecheckEvidence {
  status:
    | "complete"
    | "representative-skipped"
    | "not-in-final-distance-list"
    | "profile-unavailable";
  requested: boolean;
  representativeSkipped: boolean;
  profileAvailable: boolean;
  thresholdMode: "native-lowp-times-100000-lift";
  localPValueCutoff: number;
  emittedSignalCount: number;
  candidateSignalCount: number;
  overlappingSignalCount: number;
  eventRedetected: boolean;
  significant: boolean;
  bestBeginning: number | null;
  bestEnding: number | null;
  bestWrapsOrigin: boolean;
  bestOverlap: number;
  bestLocalPValue: number | null;
  bestCorrectedPValue: number | null;
}

export interface MaxChiRecheckEvidence {
  status:
    | "complete-active-unvalidated"
    | "representative-skipped"
    | "not-in-final-distance-list"
    | "profile-unavailable";
  kernel: "FastRecCheckMC2-strongest-peak";
  eventDiscoveryApplied: false;
  requested: boolean;
  representativeSkipped: boolean;
  profileAvailable: boolean;
  missingDataWindowFilterApplied: boolean;
  linearEdgeWindowFilterApplied: boolean;
  bonferroniApplied: boolean;
  correctionTests: number;
  variableSites: number;
  fixedWindowSites: number;
  halfWindow: number;
  criticalDifference: number;
  grownHalfWindow: number;
  bestPair: 0 | 1 | 2 | null;
  peakAlignmentPosition: number | null;
  maximumChiSquare: number | null;
  localPValue: number | null;
  withinTripletPValue: number | null;
  correctedPValue: number | null;
  sourceRecheckHit: boolean;
}

export interface ChimaeraRecheckEvidence {
  status:
    | "complete-active-unvalidated"
    | "representative-skipped"
    | "not-in-final-distance-list"
    | "profile-unavailable";
  kernel: "FastRecCheckChim-three-target-strongest-peak";
  eventDiscoveryApplied: false;
  requested: boolean;
  representativeSkipped: boolean;
  profileAvailable: boolean;
  missingDataWindowFilterApplied: boolean;
  linearEdgeWindowFilterApplied: boolean;
  bonferroniApplied: boolean;
  correctionTests: number;
  fixedWindowSites: number;
  targetProfilesScanned: number;
  bestTarget: 0 | 1 | 2 | null;
  informationRichSites: number;
  halfWindow: number;
  criticalDifference: number;
  grownHalfWindow: number;
  peakAlignmentPosition: number | null;
  maximumChiSquare: number | null;
  localPValue: number | null;
  withinTripletPValue: number | null;
  correctedPValue: number | null;
  sourceRecheckHit: boolean;
}

export interface GeneconvRecheckEvidence {
  status:
    | "complete-active-unvalidated"
    | "representative-skipped"
    | "not-in-final-distance-list"
    | "profile-unavailable";
  kernel: "GCXoverD-six-track-best-fragment";
  eventDiscoveryApplied: false;
  requested: boolean;
  representativeSkipped: boolean;
  profileAvailable: boolean;
  bonferroniApplied: boolean;
  ignoredIndels: true;
  overlapFilterApplied: true;
  minimumFragmentFiltersApplied: false;
  sourceSkewFilterRejected: boolean;
  correctionTests: number;
  polymorphicSites: number;
  tracksScreened: number;
  fragmentsScored: number;
  qualifiedFragments: number;
  overlapRejectedFragments: number;
  numericalFallbackTracks: number;
  bestTrack: 0 | 1 | 2 | 3 | 4 | 5 | null;
  recombinantLocal: 0 | 1 | 2 | null;
  beginning: number | null;
  ending: number | null;
  wrapsOrigin: boolean;
  fragmentScore: number | null;
  criticalScore: number | null;
  rawPValue: number | null;
  correctedPValue: number | null;
  sourceRecheckHit: boolean;
}

export interface ThreeSeqRecheckEvidence {
  status:
    | "complete-active-unvalidated"
    | "representative-skipped"
    | "not-in-final-distance-list"
    | "profile-unavailable";
  kernel: "TSXOver-Findall-two-orientations";
  eventDiscoveryApplied: false;
  requested: boolean;
  representativeSkipped: boolean;
  profileAvailable: boolean;
  correctionApplied: boolean;
  correctionTests: number;
  targetProfilesScanned: number;
  exactProbabilityEvaluations: number;
  approximateProbabilityEvaluations: number;
  qualifyingOrientations: number;
  sourceListEntries: number;
  bestTarget: 0 | 1 | 2 | null;
  bestDirection: "ascent" | "descent" | null;
  informationRichSites: number;
  parentOneMatches: number;
  parentTwoMatches: number;
  probabilityExcursion: number;
  maximumExcursion: number;
  beginning: number | null;
  ending: number | null;
  wrapsOrigin: boolean;
  rawPValue: number | null;
  correctedPValue: number | null;
  exactProbability: boolean;
  siegmundFallback: boolean;
  missingDataSplitApplied: boolean;
  sourceRecheckHit: boolean;
}

export interface BootscanRecheckEvidence {
  status:
    | "complete-active-unvalidated"
    | "representative-skipped"
    | "not-requested"
    | "profile-unavailable";
  kernel: "BSXoverM-SEQBOOT2-FastBootDistIP-DrawBSPlotsIII";
  eventDiscoveryApplied: false;
  coordinateChanging: false;
  requested: boolean;
  representativeSkipped: boolean;
  profileAvailable: boolean;
  sourceDistanceMode: true;
  sourceBinomialProbability: true;
  sourceCircularWindows: true;
  erasedWindowFilterApplied: boolean;
  bonferroniApplied: boolean;
  correctionTests: number;
  windowSites: number;
  stepSites: number;
  bootstrapReplicates: number;
  randomSeed: number;
  supportCutoff: number;
  windowsScanned: number;
  eventWindowsScored: number;
  usableEventWindows: number;
  informativeSites: number;
  tractInformativeSites: number;
  tractPairMatches: number;
  outsidePairMatches: number;
  scoredPair: 0 | 1 | 2 | null;
  maximumPairSupport: number;
  meanScoredPairSupport: number;
  localPValue: number | null;
  correctedPValue: number | null;
  supportGatePassed: boolean;
  sourceRecheckHit: boolean;
}

export interface SiscanRecheckEvidence {
  status:
    | "complete-active-unvalidated"
    | "representative-skipped"
    | "not-requested"
    | "outlier-unavailable"
    | "profile-unavailable";
  kernel: "GetSSOL-Get3Score-GetPScores2-DoPerms3P-MakeZValue2-DoSums";
  eventDiscoveryApplied: false;
  coordinateChanging: false;
  requested: boolean;
  representativeSkipped: boolean;
  profileAvailable: boolean;
  outlierAvailable: boolean;
  sourceNearestOutlier: true;
  sourceGapStripping: true;
  sourceVariablePatterns: true;
  bonferroniApplied: boolean;
  correctionTests: number;
  outlierSequence: number | null;
  informativeSites: number;
  permutationDraws: number;
  globalPair: 0 | 1 | 2;
  scoredPair: 0 | 1 | 2 | null;
  selectedScore: number;
  selectedScoreFamily: "partition" | "summed" | "unavailable";
  maximumZ: number;
  normalTailPValue: number;
  regionLengthAdjustedPValue: number;
  windowAdjustedPValue: number;
  correctedPValue: number;
  sourceRecheckHit: boolean;
}

export interface DistanceCorrelationEvidence {
  sequenceIndex: number;
  sequenceName: string;
  correlations: [number, number, number];
  directCorrelations: [number, number, number];
  pValues: [number, number, number];
  inversionCodes: [number, number, number];
  warningFiltered: [boolean, boolean, boolean];
  duplicateFiltered: [boolean, boolean, boolean];
  minimumComparableSites: [number, number, number];
  breakpointOverlapSites: [number, number];
  aggregateScore: number;
  aggregateTarget: number;
  overlapEligible: boolean;
  acceptableAffinity: boolean;
  strongCorrelationOverride: boolean;
  bestMatrixPair: number | null;
  significant: boolean;
  detectableSupport: boolean;
  positiveSupport: boolean;
  inverseSupport: boolean;
  strippedInverseOnly: boolean;
  duplicateCleanedSupport: boolean;
  finalTrimMatrix: FinalTrimMatrixEvidence;
  calcMatch: CalcMatchEvidence;
  consensusScore: ConsensusScoreEvidence;
  postGroupRdpRecheck: PostGroupRdpRecheckEvidence;
  postGroupMaxChiRecheck: MaxChiRecheckEvidence;
  postGroupChimaeraRecheck: ChimaeraRecheckEvidence;
  postGroupGeneconvRecheck: GeneconvRecheckEvidence;
  postGroupThreeSeqRecheck: ThreeSeqRecheckEvidence;
  postGroupBootscanRecheck: BootscanRecheckEvidence;
  postGroupSiscanRecheck: SiscanRecheckEvidence;
}

export interface PhylogeneticCorrelationEvidence {
  sequenceIndex: number;
  sequenceName: string;
  collapsedAffinityMargins: [number, number, number];
  rawAffinityMargins: [number, number, number];
  collapsedPairSupport: [boolean, boolean, boolean];
  rawPairSupport: [boolean, boolean, boolean];
  bestTreePair: number | null;
  supportingTreePairs: number;
  included: boolean;
  distanceFallback: boolean;
  disabledExcluded: boolean;
}

export interface TreeRegionSummary {
  name:
    | "5-prime-outside"
    | "5-prime-inside"
    | "3-prime-outside"
    | "3-prime-inside"
    | "outside-tract"
    | "inside-tract";
  sites: number;
  sequences: number;
  bootstrapReplicates: number;
  supportedInternalBranches: number;
  internalBranches: number;
  rawDistanceRankLevels: number;
  collapsedDistanceRankLevels: number;
  negativeBranchesNormalized: number;
  bootstrapRandomSeed: number;
  usable: boolean;
}

export interface TreePanelSummary {
  sequenceCount: number;
  subsampled: boolean;
  sequenceCap: number;
  njKernel: "supplied-clearcut-float";
  distanceEncoding: "source-midpoint-ultrametric-ranks";
  bootstrapGenerator: "microsoft-crt-seqboot2";
  bootstrapSupport: "base-tree-pseudocount";
  negativeBranchPolicy: "absolute-five-decimal-serialization";
  analyticalBranchParsing: "four-decimal-clamped-complete-edge-repair";
  treeRooting: "source-midpoint-ultrametric";
  collapseEncoding: "parent-rank-promotion-no-recompression";
  randomSeed: number;
  flankVariableSiteTarget: number;
  regions: TreeRegionSummary[];
}

export interface RoleMetricEvidence {
  method:
    | "PhPr"
    | "TreePhPr"
    | "CollapsedTreePhPr"
    | "SubPhPr"
    | "TreeSubPhPr"
    | "SubDist"
    | "TreeSubDist"
    | "TrpScore"
    | "ThreeSetSupport";
  scores: [number, number, number];
  contributions: [number, number, number];
  weight: number;
  winningRole: number | null;
  higherIsRecombinant: boolean;
  informative: boolean;
}

export interface RoleConsensusEvidence {
  method: "source-decision-tree-subset";
  nativeWeightParity: false;
  informative: boolean;
  recommendedRole: number | null;
  recommendedRecombinant: number;
  recommendedRecombinantName: string;
  recommendedMajorParent: number;
  recommendedMajorParentName: string;
  recommendedMinorParent: number;
  recommendedMinorParentName: string;
  confidence: number;
  votes: [number, number, number];
  metrics: RoleMetricEvidence[];
}

export interface RoleHypothesisEvidence {
  presumedRecombinant: number;
  presumedRecombinantName: string;
  parentOne: number;
  parentOneName: string;
  parentTwo: number;
  parentTwoName: string;
  testedSequences: number;
  validSequences: number;
  detectableSignalSetIndices: number[];
  detectableSignalSetNames: string[];
  distanceCorrelationSetIndices: number[];
  distanceCorrelationSetNames: string[];
  phylogeneticCorrelationSetIndices: number[];
  phylogeneticCorrelationSetNames: string[];
  completeTwoOfThreeSetIndices: number[];
  completeTwoOfThreeSetNames: string[];
  correlationWarnings: [boolean, boolean, boolean];
  distanceCorrelationEvidence: DistanceCorrelationEvidence[];
  phylogeneticCorrelationEvidence: PhylogeneticCorrelationEvidence[];
  phylogeneticCorrelationStatus: "complete";
  evidenceSetConsensusComplete: true;
  finalTrimDuplicateCorrelationStatus: "complete";
  finalTrimMatrixStatus: "complete-active-rff0";
  finalTrimMembershipStatus: "complete-active-rff0";
  calcMatchStatus: "complete-active-rff0-with-bounded-unavailable-cases";
  consensusScoreStatus: "complete-active-rff0";
  selectedTreeCleanupStatus: "complete-active";
  nativeGroupMembershipComplete: true;
  primaryRdpPostGroupRecheckStatus: "complete-active";
  nativePrimaryRdpRecheckComplete: true;
  maxChiPostGroupRecheckStatus: "source-shaped-strongest-peak-unvalidated";
  nativeMaxChiFullRecheckComplete: false;
  geneconvPostGroupRecheckStatus: "source-shaped-six-track-best-fragment-unvalidated";
  nativeGeneconvFullRecheckComplete: false;
  threeSeqPostGroupRecheckStatus: "source-shaped-findall-two-orientation-unvalidated";
  nativeThreeSeqFullRecheckComplete: false;
  bootscanPostGroupRecheckStatus: "source-shaped-distance-bootstrap-binomial-unvalidated";
  nativeBootscanFullRecheckComplete: false;
  siscanPostGroupRecheckStatus: "source-shaped-fixed-region-vertical-permutation-unvalidated";
  nativeSiscanFullRecheckComplete: false;
  lateNativeConsensusComplete: false;
}

export interface ReconciledEvent {
  id: number;
  anchorSignalId: number;
  anchorMethod: DiscoveryMethod;
  detectionMethods: DiscoveryMethod[];
  maxChiChimaeraOnlySupport: boolean;
  detectionRound: number;
  erasedNucleotideSites: number;
  erasedWorkingSites: number;
  fragmentSequencesAdded: number;
  fragmentAssistedDetection: boolean;
  tractErasedForDetection: boolean;
  reconciliationBasis: "two-shared-sequences-and-30-percent-overlap";
  recombinant: number;
  recombinantName: string;
  queryReferenceInputRole: QueryReferenceInputRole;
  referenceGroup: number | null;
  majorParent: number;
  majorParentName: string;
  minorParent: number;
  minorParentName: string;
  beginning: number;
  ending: number;
  wrapsOrigin: boolean;
  breakpointConfidence: BreakpointConfidenceEvidence;
  breakpointContext: BreakpointContext;
  maxChiTripletRecheck: MaxChiRecheckEvidence;
  chimaeraTripletRecheck: ChimaeraRecheckEvidence;
  geneconvTripletRecheck: GeneconvRecheckEvidence;
  threeSeqTripletRecheck: ThreeSeqRecheckEvidence;
  bootscanTripletRecheck: BootscanRecheckEvidence;
  siscanTripletRecheck: SiscanRecheckEvidence;
  bestLocalPValue: number;
  bestCorrectedPValue: number;
  supportSignalIds: number[];
  detectableSequenceIndices: number[];
  detectableSequenceNames: string[];
  roleCandidateIndices: RoleCandidateIndices;
  automaticCoRecombinantSequenceIndices: number[];
  automaticCoRecombinantSequenceNames: string[];
  coRecombinantSequenceIndices: number[];
  coRecombinantSequenceNames: string[];
  treePanel: TreePanelSummary;
  roleConsensus: RoleConsensusEvidence;
  roleHypotheses: [RoleHypothesisEvidence, RoleHypothesisEvidence, RoleHypothesisEvidence];
  traceEvidence: TraceEvidence[];
  reviewState: ReviewState;
  manualAdjusted: boolean;
  groupManualAdjusted: boolean;
  rolesProvisional: true;
}

export interface EventEdit {
  recombinant: number;
  majorParent: number;
  minorParent: number;
  beginning: number;
  ending: number;
}

export interface SignalPlotPoint {
  alignmentPosition: number;
  pair12: number;
  pair13: number;
  pair23: number;
}

export interface SignalPlot {
  signalId: number;
  windowSites: number;
  method: DiscoveryMethod;
  metric: "pair-identity" | "chi-square" | "negative-log10-p-value" | "random-walk-height" | "bootstrap-support" | "sister-scan-z-score";
  profileContext: "detection-alignment" | "original-alignment-reconstruction";
  detectionProfileExact: boolean;
  minimumValue: number;
  maximumValue: number;
  points: SignalPlotPoint[];
}

export type EventAlignmentRole =
  | "recombinant"
  | "major-parent"
  | "minor-parent"
  | "co-recombinant"
  | "evidence";

export interface EventAlignmentPanel {
  name: "beginning" | "ending";
  center: number;
  centerIndex: number;
  statisticalConfidence: BreakpointConfidenceBoundary;
  parentTransition: {
    expectedLeft: "major" | "minor";
    expectedRight: "major" | "minor";
    leftInformativeCoordinate: number | null;
    rightInformativeCoordinate: number | null;
    spanSites: number | null;
    supported: boolean;
  };
  adjacentErasureEventIds: number[];
  erasureAdjacent: boolean;
  uncertainErasureEventIds: number[];
  erasureWithinRdpWindow: boolean;
  uncertainDueToErasure: boolean;
  nativeCheckEndsApplied: boolean;
  nativeCheckEndsWarning: boolean;
  informationProfileAvailable: boolean;
  inputMissingDataInCheckRange: boolean;
  linearEdgeWithinRdpWindow: boolean;
  nativeCheckRange: {
    beginning: number;
    ending: number;
    wrapsOrigin: boolean;
    coordinateCount: number;
  };
  rdpWindowInformativeSites: number;
  nearestErasureInformativeSites: number | null;
  coordinates: number[];
}

export interface EventAlignmentRow {
  sequenceIndex: number;
  sequenceName: string;
  role: EventAlignmentRole;
  queryReferenceInputRole: QueryReferenceInputRole;
  referenceGroup: number | null;
  masked: boolean;
  disabled: boolean;
  currentGroupMember: boolean;
  automaticGroupMember: boolean;
  trace: boolean;
  panels: [string, string];
}

export interface EventAlignmentView {
  eventId: number;
  alignmentLength: number;
  circular: boolean;
  fragmentAssisted: boolean;
  requestedFlankSites: number;
  candidateRowCount: number;
  omittedRowCount: number;
  rows: EventAlignmentRow[];
  panels: [EventAlignmentPanel, EventAlignmentPanel];
}

export interface EventTreeLeaf {
  node: number;
  workingSequenceIndex: number;
  sequenceIndex: number;
  sequenceName: string;
  fragmentEventId: number | null;
  role: EventAlignmentRole;
  queryReferenceInputRole: QueryReferenceInputRole;
  referenceGroup: number | null;
  masked: boolean;
  disabled: boolean;
  currentGroupMember: boolean;
  automaticGroupMember: boolean;
  trace: boolean;
}

export interface EventTreeEdge {
  from: number;
  to: number;
  length: number;
  bootstrapSupport: number | null;
  internal: boolean;
  collapsed: boolean;
}

export interface EventTreeRegion {
  name: TreeRegionSummary["name"];
  sites: number;
  sequences: number;
  usable: boolean;
  nodeCount: number;
  root: number;
  bootstrapReplicates: number;
  supportedInternalBranches: number;
  internalBranches: number;
  rawDistanceRankLevels: number;
  collapsedDistanceRankLevels: number;
  negativeBranchesNormalized: number;
  bootstrapRandomSeed: number;
  edges: EventTreeEdge[];
}

export interface EventTreeView {
  eventId: number;
  method: "neighbour-joining";
  distance: "Jukes-Cantor";
  njKernel: "supplied-clearcut-float";
  distanceEncoding: "source-midpoint-ultrametric-ranks";
  bootstrapGenerator: "microsoft-crt-seqboot2";
  bootstrapSupport: "base-tree-pseudocount";
  negativeBranchPolicy: "absolute-five-decimal-serialization";
  analyticalBranchParsing: "four-decimal-clamped-complete-edge-repair";
  treeRooting: "source-midpoint-ultrametric";
  collapseEncoding: "parent-rank-promotion-no-recompression";
  displayRooting: "arbitrary-internal-node";
  bootstrapCollapseCutoff: number;
  bootstrapReplicates: number;
  randomSeed: number;
  flankVariableSiteTarget: number;
  subsampled: boolean;
  sequenceCap: number;
  fragmentAssisted: boolean;
  leaves: EventTreeLeaf[];
  regions: EventTreeRegion[];
}

export interface PhylproPlotPoint {
  alignmentPosition: number;
  recombinant: number;
  majorParent: number;
  minorParent: number;
}

export interface PhylproMinimum {
  sequenceIndex: number;
  position: number;
  correlation: number;
}

export interface PhylproBreakpoint {
  name: "beginning" | "ending";
  eventPosition: number;
  profilePosition: number;
  correlations: [number, number, number];
}

export interface EventPhylproView {
  eventId: number;
  status: "source-shaped-active-unvalidated";
  kernel: "FindSubSeqPP-MakePDstMat-UpdatePDstMat-PPRegression";
  roleOrder: "recombinant-major-parent-minor-parent";
  columnSelection: "polymorphic-after-gap-policy";
  distance: "pairwise-Hamming-count";
  correlation: "Pearson-source-single-output";
  significanceTest: "not-implemented-in-supplied-rdp5";
  optimization: "three-target-rows-linear-in-context";
  maskedContextIncluded: true;
  disabledContextExcluded: true;
  fragmentContextIncluded: false;
  circular: boolean;
  windowSites: number;
  halfWindowSites: number;
  windowCapped: boolean;
  gapMode: PhylproGapMode;
  includeSelf: boolean;
  eligibleColumns: number;
  contextSequences: number;
  targetContextComparisons: number;
  rollingUpdates: number;
  evaluatedPoints: number;
  returnedPoints: number;
  minimumValue: number;
  maximumValue: number;
  sequenceIndices: [number, number, number];
  sequenceNames: [string, string, string];
  minimumBySequence: [PhylproMinimum, PhylproMinimum, PhylproMinimum];
  breakpoints: [PhylproBreakpoint, PhylproBreakpoint];
  points: PhylproPlotPoint[];
}

export interface LateConsensusStatus {
  status: "active-rdp-geneconv-bootscan-maxchi-chimaera-siscan-threeseq-plus-optional-bootscan-siscan-post-group-recheck";
  groupPruningApplied: true;
  nativeGroupMembershipComplete: true;
  primaryRdpPostGroupRecheckApplied: true;
  nativePrimaryRdpRecheckComplete: true;
  maxChiTripletRecheckApplied: true;
  maxChiPostGroupRecheckApplied: true;
  maxChiKernelStatus: "source-shaped-multi-peak-destroy-retry-unvalidated";
  maxChiEventDiscoveryApplied: boolean;
  maxChiDiscoveryFeedsCyclicScheduler: boolean;
  chimaeraKernelStatus: "source-shaped-target-profile-multi-peak-destroy-retry-unvalidated";
  chimaeraEventDiscoveryApplied: boolean;
  chimaeraDiscoveryFeedsCyclicScheduler: boolean;
  chimaeraTripletRecheckApplied: true;
  chimaeraPostGroupRecheckApplied: true;
  chimaeraRecheckKernelStatus: "source-shaped-three-target-strongest-peak-unvalidated";
  geneconvKernelStatus: "source-shaped-six-track-ka-fragments-unvalidated";
  geneconvEventDiscoveryApplied: boolean;
  geneconvDiscoveryFeedsCyclicScheduler: boolean;
  geneconvTripletRecheckApplied: true;
  geneconvPostGroupRecheckApplied: true;
  geneconvRecheckKernelStatus: "source-shaped-six-track-best-fragment-unvalidated";
  threeSeqKernelStatus: "source-shaped-hypergeometric-random-walk-unvalidated";
  threeSeqEventDiscoveryApplied: boolean;
  threeSeqDiscoveryFeedsCyclicScheduler: boolean;
  threeSeqTripletRecheckApplied: true;
  threeSeqPostGroupRecheckApplied: true;
  threeSeqRecheckKernelStatus: "source-shaped-findall-two-orientation-unvalidated";
  bootscanPrimaryEnabled: boolean;
  bootscanEventDiscoveryApplied: boolean;
  bootscanDiscoveryFeedsCyclicScheduler: boolean;
  bootscanDiscoveryKernelStatus: "source-shaped-bsxoverr-distance-bootstrap-binomial-unvalidated";
  bootscanSecondaryEnabled: boolean;
  bootscanTripletRecheckApplied: boolean;
  bootscanPostGroupRecheckApplied: boolean;
  bootscanRecheckKernelStatus: "source-shaped-distance-bootstrap-binomial-unvalidated";
  nativeBootscanFullRecheckComplete: false;
  siscanPrimaryEnabled: boolean;
  siscanEventDiscoveryApplied: boolean;
  siscanDiscoveryFeedsCyclicScheduler: boolean;
  siscanDiscoveryKernelStatus: "source-shaped-ssxoverc-wpgma-vertical-permutation-unvalidated";
  siscanSecondaryEnabled: boolean;
  siscanTripletRecheckApplied: boolean;
  siscanPostGroupRecheckApplied: boolean;
  siscanRecheckKernelStatus: "source-shaped-fixed-region-vertical-permutation-unvalidated";
  nativeSiscanFullRecheckComplete: false;
  nativeMaxChiFullRecheckComplete: false;
  nativeChimaeraFullRecheckComplete: false;
  nativeGeneconvFullRecheckComplete: false;
  nativeThreeSeqFullRecheckComplete: false;
  implementedStages: string[];
  pendingStages: string[];
}

export interface BreakpointInspectionStatus {
  available: true;
  source: "original-alignment";
  maxFlankSites: number;
  maxRows: number;
  nativeCheckEndsStatus: "complete-active-unvalidated";
  nativeCheckEndsAfterFirstEvent: true;
  inputMissingRunLength: 10;
  uncertaintyReasons: [
    "prior-erasure",
    "input-missing-data",
    "linear-edge",
    "profile-unavailable",
  ];
  statisticalConfidence: {
    available: true;
    status: "complete-active-unvalidated";
    method: "BURT/BenHMM";
    hmmCyclesArgument: 20;
    serialTrainingStarts: 21;
    posteriorThresholds: [0.995, 0.999];
    randomSeed: 3;
    randomAdapter: "msvc-rand-15-bit";
    optionAvailable: true;
    enabledByDefault: true;
    canRepositionDetectedEvents: true;
  };
}

export interface TreeInspectionStatus {
  available: true;
  source: "reconciliation-tree-panel";
  regionCount: 6;
  bootstrapCollapseCutoff: number;
  payload: "on-demand-edge-lists";
}

export interface PhylproInspectionStatus {
  available: true;
  status: "source-shaped-active-unvalidated";
  source: "original-alignment-current-event-roles";
  payload: "on-demand-three-target-correlation-profile";
  defaultWindowSites: 60;
  defaultGapMode: "ignore-missing-pairwise";
  defaultIncludeSelf: false;
  significanceTest: "not-implemented-in-supplied-rdp5";
  canChangeEvents: false;
}

export interface ScanResults {
  engineVersion: string;
  status: "cyclic-three-set-reconciled";
  method: string;
  analysisMode: AnalysisMode;
  queryReference: {
    active: boolean;
    querySequenceCount: number;
    referenceSequenceCount: number;
    referenceGroupCount: number;
    tripletConstraint: "one-query-two-different-reference-groups";
    sourceCorrectionRule: "reference-group-pairs-times-query-origins";
  };
  discoveryMethods: DiscoveryMethod[];
  reconciliationTier: "detectable-distance-phylogenetic";
  cycleMode: "strongest-first-tract-erasure-with-bounded-fragment-reentry";
  lateConsensus: LateConsensusStatus;
  breakpointInspection: BreakpointInspectionStatus;
  treeInspection: TreeInspectionStatus;
  phylproInspection: PhylproInspectionStatus;
  finalAlignmentReady: boolean;
  fragmentReentry: boolean;
  fragmentReentryAlignmentLengthLimit: number;
  fragmentSequenceCap: number;
  fragmentReentryCapped: boolean;
  workingSequenceCount: number;
  workingFragmentSequenceCount: number;
  activeWorkingSequenceCount: number;
  queryWorkingSequenceCount: number;
  referenceWorkingSequenceCount: number;
  activeReferenceGroupCount: number;
  processedTriplets: number;
  totalTriplets: number;
  scanRounds: number;
  cumulativeTriplets: number;
  maxChiProfilesScanned: number;
  maxChiPeakAttempts: number;
  maxChiCandidatesFound: number;
  maxChiPeakLimitTriplets: number;
  chimaeraProfilesScanned: number;
  chimaeraPeakAttempts: number;
  chimaeraCandidatesFound: number;
  chimaeraPeakLimitTargets: number;
  geneconvFragmentsScored: number;
  geneconvQualifiedFragments: number;
  geneconvCandidatesFound: number;
  geneconvOverlapRejections: number;
  geneconvNumericalFallbackTracks: number;
  threeSeqProfilesScanned: number;
  threeSeqExactEvaluations: number;
  threeSeqApproximateEvaluations: number;
  threeSeqCandidatesFound: number;
  bootscanProfilesScanned: number;
  bootscanCandidateRegionsScored: number;
  bootscanCandidatesFound: number;
  bootscanPairProfilesRequested: number;
  bootscanPairProfileCacheHits: number;
  bootscanPairProfileCacheMisses: number;
  bootscanPairProfileCacheEvictions: number;
  bootscanPairProfileCachePeakBytes: number;
  siscanProfilesScanned: number;
  siscanWindowsScored: number;
  siscanCandidateRegionsScored: number;
  siscanCandidatesFound: number;
  siscanPermutationDraws: number;
  siscanContextBuilds: number;
  siscanContextPairComparisons: number;
  siscanContextTreeMerges: number;
  siscanRandomValuesGenerated: number;
  cycleTermination: string;
  correction: CorrectionMode;
  correctionTests: number;
  circular: boolean;
  maskedSequenceIndices: number[];
  disabledSequenceIndices: number[];
  referenceGroupIndices: number[];
  downstreamReconciliationRequiredAfter: number | null;
  pValueCutoff: number;
  windowSites: number;
  maxChiEnabled: boolean;
  maxChiWindowSites: number;
  chimaeraEnabled: boolean;
  chimaeraWindowSites: number;
  geneconvEnabled: boolean;
  geneconvMismatchScale: number;
  geneconvMaxOverlaps: number;
  threeSeqEnabled: boolean;
  bootscanPrimaryEnabled: boolean;
  bootscanSecondaryEnabled: boolean;
  bootscanWindowSites: number;
  bootscanStepSites: number;
  bootscanBootstrapReplicates: number;
  bootscanSupportCutoff: number;
  bootscanRandomSeed: number;
  siscanPrimaryEnabled: boolean;
  siscanSecondaryEnabled: boolean;
  siscanWindowSites: number;
  siscanStepSites: number;
  siscanScanPermutations: number;
  siscanPValuePermutations: number;
  siscanRandomSeed: number;
  polishBreakpoints: boolean;
  timing: ScanTiming | null;
  execution: ScanExecution;
  signals: RdpSignal[];
  events: ReconciledEvent[];
  notes: string[];
}

export interface ImportedProject {
  dataset: DatasetSummary;
  results: ScanResults | null;
  sourceFilename: string;
}

export type WorkerRequest =
  | { id: number; type: "init"; wasmBaseUrl: string; assetVersion: string }
  | { id: number; type: "load"; name: string; bytes: ArrayBuffer }
  | { id: number; type: "import-project"; name: string; bytes: ArrayBuffer }
  | { id: number; type: "scan"; options: ScanOptions }
  | { id: number; type: "cancel" }
  | { id: number; type: "plot"; signalId: number }
  | { id: number; type: "event-alignment"; eventId: number; flankSites: number; rowLimit: number }
  | { id: number; type: "event-trees"; eventId: number }
  | {
      id: number;
      type: "event-phylpro";
      eventId: number;
      windowSites: number;
      gapMode: PhylproGapMode;
      includeSelf: boolean;
    }
  | { id: number; type: "set-review-state"; signalId: number; state: ReviewState }
  | { id: number; type: "set-event-review-state"; eventId: number; state: ReviewState }
  | { id: number; type: "update-event"; eventId: number; edit: EventEdit }
  | { id: number; type: "update-event-group"; eventId: number; sequenceIndices: number[]; manualOverride: boolean }
  | { id: number; type: "reconcile-after"; eventId: number; cpuThreads: number }
  | { id: number; type: "export-csv" }
  | {
      id: number;
      type: "export-enabled-sequences";
      maskedSequenceIndices: number[];
      disabledSequenceIndices: number[];
    }
  | {
      id: number;
      type: "export-masked-or-disabled-sequences";
      maskedSequenceIndices: number[];
      disabledSequenceIndices: number[];
    }
  | { id: number; type: "export-recombinant-sequences-removed" }
  | { id: number; type: "export-recombinant-columns-removed" }
  | { id: number; type: "export-recombination-free" }
  | { id: number; type: "export-fragmented" }
  | { id: number; type: "export-project" };

export type WorkerEvent =
  | { type: "progress"; progress: ScanProgress }
  | ({ type: "engine" } & EngineRuntimeInfo);

export type WorkerResponse =
  | { id: number; ok: true; value: unknown }
  | { id: number; ok: false; error: string };

type WithoutRequestId<Request> = Request extends { id: number }
  ? Omit<Request, "id">
  : never;

export type WorkerRequestPayload = WithoutRequestId<WorkerRequest>;
