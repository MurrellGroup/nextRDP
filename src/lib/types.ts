export type WorkflowStep = "dataset" | "settings" | "scan" | "review" | "export";

export type CorrectionMode = "bonferroni" | "none";

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
  maskedSequenceIndices: number[];
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
  maskedExcluded: boolean;
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

export interface ScanResults {
  engineVersion: string;
  status: "cyclic-three-set-reconciled";
  method: "RDP";
  reconciliationTier: "detectable-distance-phylogenetic";
  cycleMode: "strongest-first-tract-erasure-with-bounded-fragment-reentry";
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
  downstreamReconciliationRequiredAfter: number | null;
  pValueCutoff: number;
  windowSites: number;
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
  | { id: number; type: "init"; wasmBaseUrl: string }
  | { id: number; type: "load"; name: string; bytes: ArrayBuffer }
  | { id: number; type: "import-project"; name: string; bytes: ArrayBuffer }
  | { id: number; type: "scan"; options: ScanOptions }
  | { id: number; type: "cancel" }
  | { id: number; type: "plot"; signalId: number }
  | { id: number; type: "set-review-state"; signalId: number; state: ReviewState }
  | { id: number; type: "set-event-review-state"; eventId: number; state: ReviewState }
  | { id: number; type: "update-event"; eventId: number; edit: EventEdit }
  | { id: number; type: "update-event-group"; eventId: number; sequenceIndices: number[]; manualOverride: boolean }
  | { id: number; type: "reconcile-after"; eventId: number }
  | { id: number; type: "export-csv" }
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
