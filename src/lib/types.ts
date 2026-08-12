export type WorkflowStep = "dataset" | "settings" | "scan" | "review" | "export";

export type CorrectionMode = "bonferroni" | "none";
export type SequenceAnalysisState = "enabled" | "masked" | "disabled";

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
  circular: boolean;
  pValueCutoff: number;
  correction: CorrectionMode;
  windowSites: number;
  polishBreakpoints: boolean;
  maskedSequenceIndices: number[];
  disabledSequenceIndices: number[];
}

export interface ScanProgress {
  state: "idle" | "running" | "done" | "cancelled" | "error";
  phase: "primary" | "cyclic-rescan" | "reconciliation" | "complete";
  processedTriplets: number;
  totalTriplets: number;
  cumulativeTriplets: number;
  scanRound: number;
  fixedEventCount: number;
  signalCount: number;
  eventCount: number;
  cycleTermination: string;
  fraction: number;
}

export type ReviewState = "unreviewed" | "accepted" | "rejected";

export interface RdpSignal {
  id: number;
  method: "RDP";
  triplet: [number, number, number];
  tripletNames: [string, string, string];
  recombinant: number;
  recombinantName: string;
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
  fixedWindowSites: 70;
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
  usable: boolean;
}

export interface TreePanelSummary {
  sequenceCount: number;
  subsampled: boolean;
  sequenceCap: number;
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
  lateNativeConsensusComplete: false;
}

export interface ReconciledEvent {
  id: number;
  anchorSignalId: number;
  detectionRound: number;
  erasedNucleotideSites: number;
  erasedWorkingSites: number;
  fragmentSequencesAdded: number;
  fragmentAssistedDetection: boolean;
  tractErasedForDetection: boolean;
  reconciliationBasis: "two-shared-sequences-and-30-percent-overlap";
  recombinant: number;
  recombinantName: string;
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
  edges: EventTreeEdge[];
}

export interface EventTreeView {
  eventId: number;
  method: "neighbour-joining";
  distance: "Jukes-Cantor";
  displayRooting: "arbitrary-internal-node";
  bootstrapCollapseCutoff: number;
  bootstrapReplicates: number;
  subsampled: boolean;
  sequenceCap: number;
  fragmentAssisted: boolean;
  leaves: EventTreeLeaf[];
  regions: EventTreeRegion[];
}

export interface LateConsensusStatus {
  status: "active-rdp-maxchi-post-group-recheck";
  groupPruningApplied: true;
  nativeGroupMembershipComplete: true;
  primaryRdpPostGroupRecheckApplied: true;
  nativePrimaryRdpRecheckComplete: true;
  maxChiTripletRecheckApplied: true;
  maxChiPostGroupRecheckApplied: true;
  maxChiKernelStatus: "source-shaped-strongest-peak-unvalidated";
  maxChiEventDiscoveryApplied: false;
  nativeMaxChiFullRecheckComplete: false;
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

export interface ScanResults {
  engineVersion: string;
  status: "cyclic-three-set-reconciled";
  method: "RDP";
  reconciliationTier: "detectable-distance-phylogenetic";
  cycleMode: "strongest-first-tract-erasure-with-bounded-fragment-reentry";
  lateConsensus: LateConsensusStatus;
  breakpointInspection: BreakpointInspectionStatus;
  treeInspection: TreeInspectionStatus;
  finalAlignmentReady: boolean;
  fragmentReentry: boolean;
  fragmentReentryAlignmentLengthLimit: number;
  fragmentSequenceCap: number;
  fragmentReentryCapped: boolean;
  workingSequenceCount: number;
  workingFragmentSequenceCount: number;
  scanRounds: number;
  cumulativeTriplets: number;
  cycleTermination: string;
  correction: CorrectionMode;
  correctionTests: number;
  circular: boolean;
  maskedSequenceIndices: number[];
  disabledSequenceIndices: number[];
  downstreamReconciliationRequiredAfter: number | null;
  pValueCutoff: number;
  windowSites: number;
  polishBreakpoints: boolean;
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
  | { id: number; type: "set-review-state"; signalId: number; state: ReviewState }
  | { id: number; type: "set-event-review-state"; eventId: number; state: ReviewState }
  | { id: number; type: "update-event"; eventId: number; edit: EventEdit }
  | { id: number; type: "update-event-group"; eventId: number; sequenceIndices: number[]; manualOverride: boolean }
  | { id: number; type: "reconcile-after"; eventId: number }
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
  | { type: "engine"; threaded: boolean; version: string };

export type WorkerResponse =
  | { id: number; ok: true; value: unknown }
  | { id: number; ok: false; error: string };

type WithoutRequestId<Request> = Request extends { id: number }
  ? Omit<Request, "id">
  : never;

export type WorkerRequestPayload = WithoutRequestId<WorkerRequest>;
