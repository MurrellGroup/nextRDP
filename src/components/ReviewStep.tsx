import {
  AlertTriangle,
  ArrowLeft,
  ArrowRight,
  Check,
  ChevronLeft,
  ChevronRight,
  GitCompareArrows,
  GitBranch,
  Layers3,
  Pencil,
  RefreshCw,
  RotateCcw,
  Save,
  Search,
  X,
} from "lucide-react";
import { useEffect, useMemo, useState } from "react";

import type {
  BreakpointErasureBoundary,
  EventAlignmentView,
  EventEdit,
  EventTreeView,
  EventPhylproView,
  PhylproGapMode,
  ReconciledEvent,
  ReviewState,
  ScanResults,
  SequenceSummary,
  SignalPlot as SignalPlotData,
} from "../lib/types";
import { EventAlignmentInspector } from "./EventAlignmentInspector";
import { EventTreeInspector } from "./EventTreeInspector";
import { EventPhylproInspector } from "./EventPhylproInspector";
import { SignalPlot } from "./SignalPlot";

interface ReviewStepProps {
  results: ScanResults;
  alignmentLength: number;
  sequences: SequenceSummary[];
  onGetPlot: (signalId: number) => Promise<SignalPlotData>;
  onGetEventAlignment: (
    eventId: number,
    flankSites: number,
    rowLimit: number,
  ) => Promise<EventAlignmentView>;
  onGetEventTrees: (eventId: number) => Promise<EventTreeView>;
  onGetEventPhylpro: (
    eventId: number,
    windowSites: number,
    gapMode: PhylproGapMode,
    includeSelf: boolean,
  ) => Promise<EventPhylproView>;
  onReviewState: (eventId: number, state: ReviewState) => void;
  onUpdateEvent: (eventId: number, edit: EventEdit) => Promise<void>;
  onUpdateEventGroup: (
    eventId: number,
    sequenceIndices: number[],
    manualOverride?: boolean,
  ) => Promise<void>;
  onReconcileAfter: (eventId: number) => void;
  reconciling: boolean;
  onSaveProject: () => void;
  checkpointDirty: boolean;
  checkpointSaving: boolean;
  onBack: () => void;
  onExport: () => void;
}

function pValue(value: number): string {
  if (value <= 1e-300) return "<1 × 10⁻³⁰⁰";
  if (value < 0.001) return value.toExponential(3).replace("e", " × 10^");
  return value.toPrecision(4);
}

function roleScore(value: number): string {
  if (Math.abs(value) >= 1000 || (Math.abs(value) > 0 && Math.abs(value) < 0.001)) {
    return value.toExponential(2);
  }
  return value.toFixed(3).replace(/\.?0+$/, "");
}

function confidenceCoordinate(value: number): string {
  return value === -1 ? "—" : Math.abs(value).toLocaleString();
}

function breakpointUncertaintyReason(
  name: "beginning" | "ending",
  boundary: BreakpointErasureBoundary,
): string {
  const reasons: string[] = [];
  if (boundary.uncertainDueToErasure) {
    if (boundary.erasureAdjacent) {
      reasons.push(`The ${name} breakpoint touches sequence erased during an earlier cyclic pass.`);
    } else if (boundary.nearestErasureInformativeSites === 0) {
      reasons.push(
        `The ${name} breakpoint has no intervening RDP information-rich positions before earlier erased sequence.`,
      );
    } else {
      reasons.push(
        `The ${name} breakpoint's supplied CheckEnds range reaches earlier erased sequence within ${
          boundary.nearestErasureInformativeSites ?? boundary.rdpWindowInformativeSites
        } RDP information-rich position${
          boundary.nearestErasureInformativeSites === 1 ? "" : "s"
        }.`,
      );
    }
  }
  if (boundary.inputMissingDataInCheckRange) {
    reasons.push("The same native range contains a source-shaped input MissingData run.");
  }
  if (boundary.linearEdgeWithinRdpWindow) {
    reasons.push(
      `A complete RDP window cannot be formed ${name === "beginning" ? "before" : "after"} the breakpoint at the linear alignment edge.`,
    );
  }
  if (!boundary.informationProfileAvailable) {
    reasons.push("No usable information-rich position profile is available for the native check.");
  }
  return reasons.join(" ");
}

function EventSchematic({ event, alignmentLength }: { event: ReconciledEvent; alignmentLength: number }) {
  const percent = (position: number) => `${Math.max(0, Math.min(100, (position / alignmentLength) * 100))}%`;
  const singleWidth = Math.max(0, event.ending - event.beginning);
  return (
    <div className="event-schematic" aria-label="Schematic recombinant region">
      <div className="schematic-labels">
        <span>1</span>
        <strong>{event.recombinantName}</strong>
        <span>{alignmentLength.toLocaleString()}</span>
      </div>
      <div className="sequence-track">
        <span className="major-track" />
        {event.wrapsOrigin ? (
          <>
            <span className="minor-track" style={{ left: 0, width: percent(event.ending) }} />
            <span className="minor-track" style={{ left: percent(event.beginning), right: 0 }} />
          </>
        ) : (
          <span className="minor-track" style={{ left: percent(event.beginning), width: percent(singleWidth) }} />
        )}
        <span className="breakpoint" style={{ left: percent(event.beginning) }} />
        <span className="breakpoint" style={{ left: percent(event.ending) }} />
      </div>
      <div className="schematic-key">
        <span className="key-major">Major-parent-like</span>
        <span className="key-minor">Minor-parent-like</span>
        <span>
          {event.beginning.toLocaleString()} → {event.ending.toLocaleString()}
          {event.wrapsOrigin ? " · wraps origin" : ""}
        </span>
      </div>
    </div>
  );
}

export function ReviewStep({
  results,
  alignmentLength,
  sequences,
  onGetPlot,
  onGetEventAlignment,
  onGetEventTrees,
  onGetEventPhylpro,
  onReviewState,
  onUpdateEvent,
  onUpdateEventGroup,
  onReconcileAfter,
  reconciling,
  onSaveProject,
  checkpointDirty,
  checkpointSaving,
  onBack,
  onExport,
}: ReviewStepProps) {
  const [selectedId, setSelectedId] = useState(results.events[0]?.id ?? -1);
  const [plot, setPlot] = useState<SignalPlotData | null>(null);
  const [plotLoading, setPlotLoading] = useState(false);
  const [alignmentOpen, setAlignmentOpen] = useState(false);
  const [alignmentView, setAlignmentView] = useState<EventAlignmentView | null>(null);
  const [alignmentLoading, setAlignmentLoading] = useState(false);
  const [alignmentError, setAlignmentError] = useState("");
  const [alignmentFlankSites, setAlignmentFlankSites] = useState(30);
  const [treesOpen, setTreesOpen] = useState(false);
  const [treeView, setTreeView] = useState<EventTreeView | null>(null);
  const [treeLoading, setTreeLoading] = useState(false);
  const [treeError, setTreeError] = useState("");
  const [phylproOpen, setPhylproOpen] = useState(false);
  const [phylproView, setPhylproView] = useState<EventPhylproView | null>(null);
  const [phylproLoading, setPhylproLoading] = useState(false);
  const [phylproError, setPhylproError] = useState("");
  const [phylproWindowSites, setPhylproWindowSites] = useState(60);
  const [phylproGapMode, setPhylproGapMode] = useState<PhylproGapMode>(
    "ignore-missing-pairwise",
  );
  const [phylproIncludeSelf, setPhylproIncludeSelf] = useState(false);
  const [editing, setEditing] = useState(false);
  const [editError, setEditError] = useState("");
  const [editingGroup, setEditingGroup] = useState(false);
  const [groupDraft, setGroupDraft] = useState<number[]>([]);
  const [groupSearch, setGroupSearch] = useState("");
  const [groupError, setGroupError] = useState("");
  const maskedSequences = useMemo(
    () => new Set(results.maskedSequenceIndices),
    [results.maskedSequenceIndices],
  );
  const disabledSequences = useMemo(
    () => new Set(results.disabledSequenceIndices),
    [results.disabledSequenceIndices],
  );
  const selectedIndex = results.events.findIndex((event) => event.id === selectedId);
  const selected = results.events[selectedIndex];
  const anchor = selected
    ? results.signals.find((signal) => signal.id === selected.anchorSignalId)
    : undefined;
  const [draft, setDraft] = useState<EventEdit>({
    recombinant: selected?.recombinant ?? 0,
    majorParent: selected?.majorParent ?? 1,
    minorParent: selected?.minorParent ?? 2,
    beginning: selected?.beginning ?? 1,
    ending: selected?.ending ?? 1,
  });

  useEffect(() => {
    if (!selected) return;
    setDraft({
      recombinant: selected.recombinant,
      majorParent: selected.majorParent,
      minorParent: selected.minorParent,
      beginning: selected.beginning,
      ending: selected.ending,
    });
    setEditing(false);
    setEditError("");
    setEditingGroup(false);
    setGroupDraft(selected.coRecombinantSequenceIndices);
    setGroupSearch("");
    setGroupError("");
  }, [
    selected?.id,
    selected?.manualAdjusted,
    selected?.recombinant,
    selected?.majorParent,
    selected?.minorParent,
    selected?.beginning,
    selected?.ending,
    selected?.groupManualAdjusted,
    selected?.coRecombinantSequenceIndices,
  ]);

  useEffect(() => {
    if (!anchor) return;
    let live = true;
    setPlotLoading(true);
    setPlot(null);
    onGetPlot(anchor.id)
      .then((value) => {
        if (live) setPlot(value);
      })
      .catch(() => {
        if (live) setPlot(null);
      })
      .finally(() => {
        if (live) setPlotLoading(false);
      });
    return () => {
      live = false;
    };
  }, [anchor?.id, onGetPlot]);

  useEffect(() => {
    setAlignmentOpen(false);
    setAlignmentView(null);
    setAlignmentError("");
    setTreesOpen(false);
    setTreeView(null);
    setTreeError("");
    setPhylproOpen(false);
    setPhylproView(null);
    setPhylproError("");
  }, [selected?.id]);

  useEffect(() => {
    if (!alignmentOpen || !selected) return;
    let live = true;
    setAlignmentLoading(true);
    setAlignmentView(null);
    setAlignmentError("");
    onGetEventAlignment(selected.id, alignmentFlankSites, 28)
      .then((value) => {
        if (live) setAlignmentView(value);
      })
      .catch((caught: unknown) => {
        if (!live) return;
        setAlignmentError(
          caught instanceof Error ? caught.message : "Breakpoint alignment data was not returned.",
        );
      })
      .finally(() => {
        if (live) setAlignmentLoading(false);
      });
    return () => {
      live = false;
    };
  }, [
    alignmentFlankSites,
    alignmentOpen,
    onGetEventAlignment,
    selected?.automaticCoRecombinantSequenceIndices,
    selected?.beginning,
    selected?.coRecombinantSequenceIndices,
    selected?.ending,
    selected?.id,
    selected?.majorParent,
    selected?.minorParent,
    selected?.recombinant,
  ]);

  useEffect(() => {
    if (!treesOpen || !selected) return;
    let live = true;
    setTreeLoading(true);
    setTreeView(null);
    setTreeError("");
    onGetEventTrees(selected.id)
      .then((value) => {
        if (live) setTreeView(value);
      })
      .catch((caught: unknown) => {
        if (!live) return;
        setTreeError(caught instanceof Error ? caught.message : "Regional tree data was not returned.");
      })
      .finally(() => {
        if (live) setTreeLoading(false);
      });
    return () => {
      live = false;
    };
  }, [
    onGetEventTrees,
    selected?.automaticCoRecombinantSequenceIndices,
    selected?.beginning,
    selected?.coRecombinantSequenceIndices,
    selected?.ending,
    selected?.id,
    selected?.majorParent,
    selected?.minorParent,
    selected?.recombinant,
    treesOpen,
  ]);

  useEffect(() => {
    if (!phylproOpen || !selected) return;
    if (!Number.isInteger(phylproWindowSites) || phylproWindowSites < 10 || phylproWindowSites > 5000) {
      setPhylproLoading(false);
      setPhylproView(null);
      setPhylproError("Choose a PHYLPRO window between 10 and 5,000 sites.");
      return;
    }
    let live = true;
    setPhylproLoading(true);
    setPhylproView(null);
    setPhylproError("");
    // Changing a numeric window can produce several intermediate values. Wait
    // briefly so a large alignment does not queue redundant O(L*N) profiles.
    const timer = window.setTimeout(() => {
      onGetEventPhylpro(
        selected.id,
        phylproWindowSites,
        phylproGapMode,
        phylproIncludeSelf,
      )
        .then((value) => {
          if (live) setPhylproView(value);
        })
        .catch((caught: unknown) => {
          if (!live) return;
          setPhylproError(
            caught instanceof Error ? caught.message : "PHYLPRO review data were not returned.",
          );
        })
        .finally(() => {
          if (live) setPhylproLoading(false);
        });
    }, 180);
    return () => {
      live = false;
      window.clearTimeout(timer);
    };
  }, [
    onGetEventPhylpro,
    phylproGapMode,
    phylproIncludeSelf,
    phylproOpen,
    phylproWindowSites,
    selected?.beginning,
    selected?.ending,
    selected?.id,
    selected?.majorParent,
    selected?.minorParent,
    selected?.recombinant,
  ]);

  const counts = useMemo(
    () => ({
      accepted: results.events.filter((event) => event.reviewState === "accepted").length,
      rejected: results.events.filter((event) => event.reviewState === "rejected").length,
      unreviewed: results.events.filter((event) => event.reviewState === "unreviewed").length,
    }),
    [results.events],
  );

  const roleChoices = useMemo(() => {
    if (!selected) return [];
    const indices = new Set([
      selected.recombinant,
      selected.majorParent,
      selected.minorParent,
      ...selected.detectableSequenceIndices,
      ...selected.coRecombinantSequenceIndices,
      selected.roleConsensus.recommendedRecombinant,
      selected.roleConsensus.recommendedMajorParent,
      selected.roleConsensus.recommendedMinorParent,
    ]);
    return [...indices]
      .sort((left, right) => left - right)
      .map((index) => sequences[index])
      .filter((sequence): sequence is SequenceSummary =>
        sequence !== undefined && !disabledSequences.has(sequence.index));
  }, [disabledSequences, selected, sequences]);

  const groupMatches = useMemo(() => {
    if (!selected) return [];
    const query = groupSearch.trim().toLocaleLowerCase();
    const selectedSet = new Set(groupDraft);
    const automaticSet = new Set(selected.automaticCoRecombinantSequenceIndices);
    return sequences
      .filter((sequence) =>
        !disabledSequences.has(sequence.index) &&
        sequence.index !== selected.majorParent &&
        sequence.index !== selected.minorParent &&
        (!query ||
          sequence.name.toLocaleLowerCase().includes(query) ||
          String(sequence.index + 1).includes(query)),
      )
      .sort((left, right) => {
        const selectedDifference = Number(selectedSet.has(right.index)) - Number(selectedSet.has(left.index));
        if (selectedDifference) return selectedDifference;
        const automaticDifference = Number(automaticSet.has(right.index)) - Number(automaticSet.has(left.index));
        return automaticDifference || left.index - right.index;
      })
      .slice(0, 40);
  }, [disabledSequences, groupDraft, groupSearch, selected, sequences]);

  if (!selected) {
    return (
      <section className="step-page empty-results">
        <div className="empty-icon"><GitCompareArrows size={28} /></div>
        <h1>No significant recombination events</h1>
        <p>
          No {results.discoveryMethods.join(", ")} signal passed the current window, threshold, and correction settings
          {results.analysisMode === "query-reference" ? " in the constrained query-vs-reference triplets" : ""}.
        </p>
        {results.maxChiPeakLimitTriplets > 0 ? (
          <p>
            MaxChi reached its supplied 100-peak retry bound for {results.maxChiPeakLimitTriplets.toLocaleString()}
            triplet{results.maxChiPeakLimitTriplets === 1 ? "" : "s"}; additional raw peaks in those
            profiles were not explored.
          </p>
        ) : null}
        {results.chimaeraPeakLimitTargets > 0 ? (
          <p>
            CHIMAERA reached its supplied 100-peak retry bound for {results.chimaeraPeakLimitTargets.toLocaleString()}
            target profile{results.chimaeraPeakLimitTargets === 1 ? "" : "s"}; additional raw peaks
            in those profiles were not explored.
          </p>
        ) : null}
        <div>
          <button className="button button-quiet" type="button" onClick={onBack}>Change scan settings</button>
          <button className="button button-primary" type="button" onClick={onExport}>Export the null result</button>
        </div>
      </section>
    );
  }

  const move = (amount: number) => {
    const next = (selectedIndex + amount + results.events.length) % results.events.length;
    setSelectedId(results.events[next].id);
  };

  const saveCorrection = async () => {
    setEditError("");
    if (new Set([draft.recombinant, draft.majorParent, draft.minorParent]).size !== 3) {
      setEditError("The three sequence roles must be distinct.");
      return;
    }
    if (
      draft.beginning < 1 ||
      draft.ending < 1 ||
      draft.beginning > alignmentLength ||
      draft.ending > alignmentLength
    ) {
      setEditError(`Breakpoints must be between 1 and ${alignmentLength.toLocaleString()}.`);
      return;
    }
    try {
      await onUpdateEvent(selected.id, draft);
      setEditing(false);
    } catch {
      setEditError("The correction was not saved. Check the event values and try again.");
    }
  };

  const saveGroupCorrection = async () => {
    setGroupError("");
    const sequenceIndices = [...new Set([...groupDraft, selected.recombinant])]
      .filter((index) => index !== selected.majorParent && index !== selected.minorParent)
      .sort((left, right) => left - right);
    try {
      await onUpdateEventGroup(selected.id, sequenceIndices, true);
      setEditingGroup(false);
    } catch {
      setGroupError("The co-recombinant group correction could not be saved.");
    }
  };

  const restoreAutomaticGroup = async () => {
    setGroupError("");
    try {
      await onUpdateEventGroup(
        selected.id,
        selected.automaticCoRecombinantSequenceIndices,
        false,
      );
      setEditingGroup(false);
    } catch {
      setGroupError("The automatic two-of-three group could not be restored.");
    }
  };

  const plottedSignal = anchor
    ? {
        ...anchor,
        beginning: selected.beginning,
        ending: selected.ending,
        wrapsOrigin: selected.wrapsOrigin,
      }
    : null;
  const pendingEvent = results.downstreamReconciliationRequiredAfter;
  const canReconcile = pendingEvent === selected.id && selected.reviewState !== "unreviewed";
  const nextUnreviewedEvent = results.events.find(
    (event) => event.reviewState === "unreviewed",
  )?.id ?? null;
  const decisionBlocked = pendingEvent !== null
    ? pendingEvent !== selected.id
    : nextUnreviewedEvent !== null && selected.id > nextUnreviewedEvent;
  const currentHypothesis = selected.roleHypotheses[0];
  const breakpointConfidence = selected.breakpointConfidence;
  const maxChiRecheck = selected.maxChiTripletRecheck;
  const chimaeraRecheck = selected.chimaeraTripletRecheck;
  const geneconvRecheck = selected.geneconvTripletRecheck;
  const threeSeqRecheck = selected.threeSeqTripletRecheck;
  const siscanRecheck = selected.siscanTripletRecheck;
  const geneconvRecheckRecombinantName = geneconvRecheck.recombinantLocal === null
    ? null
    : [selected.recombinantName, selected.majorParentName, selected.minorParentName][geneconvRecheck.recombinantLocal];
  const chimaeraRecheckTargetName = chimaeraRecheck.bestTarget === null
    ? null
    : [selected.recombinantName, selected.majorParentName, selected.minorParentName][chimaeraRecheck.bestTarget];
  const threeSeqRecheckTargetName = threeSeqRecheck.bestTarget === null
    ? null
    : [selected.recombinantName, selected.majorParentName, selected.minorParentName][threeSeqRecheck.bestTarget];
  const maxChiDiscovery = anchor?.maxChiDiscovery ?? null;
  const chimaeraDiscovery = anchor?.chimaeraDiscovery ?? null;
  const geneconvDiscovery = anchor?.geneconvDiscovery ?? null;
  const threeSeqDiscovery = anchor?.threeSeqDiscovery ?? null;
  const bootscanDiscovery = anchor?.bootscanDiscovery ?? null;
  const siscanDiscovery = anchor?.siscanDiscovery ?? null;
  const chimaeraParentOneLocal = chimaeraDiscovery
    ? ([1, 2, 0] as const)[chimaeraDiscovery.targetLocal]
    : null;
  const selectedRecombinantReferenceGroup = selected.queryReferenceInputRole === "reference"
    ? selected.referenceGroup ?? 0
    : 0;
  const beginningBreakpointUncertain =
    selected.breakpointContext.beginning.nativeCheckEndsWarning;
  const endingBreakpointUncertain =
    selected.breakpointContext.ending.nativeCheckEndsWarning;
  const breakpointUncertaintyMessage =
    beginningBreakpointUncertain && endingBreakpointUncertain
      ? "Both breakpoints trigger the supplied RDP CheckEnds uncertainty check."
      : beginningBreakpointUncertain
        ? breakpointUncertaintyReason("beginning", selected.breakpointContext.beginning)
        : breakpointUncertaintyReason("ending", selected.breakpointContext.ending);
  const recommendationDiffers =
    selected.roleConsensus.recommendedRecombinant !== selected.recombinant ||
    selected.roleConsensus.recommendedMajorParent !== selected.majorParent ||
    selected.roleConsensus.recommendedMinorParent !== selected.minorParent;
  const matrixPairLabel = (pair: number | null) => {
    if (pair === 0) return "5′ breakpoint";
    if (pair === 1) return "3′ breakpoint";
    if (pair === 2) return "tract / outside";
    return "insufficient sites";
  };
  const maxChiPairLabel = (pair: number | null) => {
    if (pair === 0) return `${selected.recombinantName} / ${selected.majorParentName}`;
    if (pair === 1) return `${selected.recombinantName} / ${selected.minorParentName}`;
    if (pair === 2) return `${selected.majorParentName} / ${selected.minorParentName}`;
    return "No screened peak";
  };
  const maxChiDiscoveryPairLabel = (pair: number | null) => {
    if (!anchor || pair === null) return "No screened peak";
    if (pair === 0) return `${anchor.tripletNames[0]} / ${anchor.tripletNames[1]}`;
    if (pair === 1) return `${anchor.tripletNames[0]} / ${anchor.tripletNames[2]}`;
    if (pair === 2) return `${anchor.tripletNames[1]} / ${anchor.tripletNames[2]}`;
    return "No screened peak";
  };
  const inputRoleLabel = (sequence: number) => {
    if (results.analysisMode !== "query-reference") return null;
    const group = results.referenceGroupIndices[sequence] ?? 0;
    return group > 0 ? `Reference group ${group}` : "Query input";
  };
  const applyRecommendation = async () => {
    setEditError("");
    try {
      await onUpdateEvent(selected.id, {
        recombinant: selected.roleConsensus.recommendedRecombinant,
        majorParent: selected.roleConsensus.recommendedMajorParent,
        minorParent: selected.roleConsensus.recommendedMinorParent,
        beginning: selected.beginning,
        ending: selected.ending,
      });
    } catch {
      setEditError("The role recommendation could not be applied.");
    }
  };

  return (
    <section className="step-page review-page" aria-labelledby="review-title">
      <header className="page-heading review-heading">
        <div>
          <span className="eyebrow">04 · Review</span>
          <h1 id="review-title">Refine the event hypothesis</h1>
          <p>
            Work in event order. Accept correct calls, repair the first material error, then
            re-identify only the later events—the review loop described in the RDP5 manual.
          </p>
        </div>
        <div className="review-heading-actions">
          <div className="review-counts">
            <span><i className="status-accepted" />{counts.accepted} accepted</span>
            <span><i className="status-rejected" />{counts.rejected} rejected</span>
            <span><i className="status-unreviewed" />{counts.unreviewed} unreviewed</span>
          </div>
          <button
            className="button button-secondary"
            type="button"
            disabled={checkpointSaving}
            onClick={onSaveProject}
            title={checkpointDirty ? "Download the latest review state" : "Download another copy of the current checkpoint"}
          >
            <Save size={15} />
            {checkpointSaving
              ? "Saving checkpoint…"
              : checkpointDirty
                ? "Save project checkpoint"
                : "Checkpoint current"}
          </button>
        </div>
      </header>

      {results.analysisMode === "query-reference" ? (
        <div className="notice notice-blue">
          <Layers3 size={18} />
          <p>
            This analysis used the manual’s automated query-vs-reference scheme: {results.queryReference.querySequenceCount.toLocaleString()} enabled queries were screened with {results.queryReference.referenceSequenceCount.toLocaleString()} references across {results.queryReference.referenceGroupCount.toLocaleString()} groups. Every primary triplet contained one query and two differently grouped references; role inference was not forced to call the query recombinant. The initial scan plan recorded {results.correctionTests.toLocaleString()} group-pair × query opportunities{results.correction === "bonferroni" ? " for Bonferroni correction" : ", although correction was disabled"}; the supplied MakeMCCorrection factor remains fixed across cyclic rounds.
          </p>
        </div>
      ) : null}

      <div className="notice notice-blue">
        <AlertTriangle size={18} />
        <p>
          {results.discoveryMethods.join(", ")} signals were ranked together in strongest-first cyclic passes, with each inferred co-group tract erased
          and re-entered as a gap-padded fragment before the next {results.analysisMode === "query-reference" ? "constrained" : "exploratory"} screen. Unchanged triplets reuse their XOverList-style summaries while affected rows and new fragments run fresh kernels. Three evidence sets are
          evaluated for every role; native PhPr, leave-one-out, displacement, collapsed-tree, and
          TrpScore decision contributions are auditable below. The project correction factor is the
          initial scan plan’s {results.correctionTests.toLocaleString()} opportunities and stays fixed
          even when fragment re-entry changes a later round’s actual workload.
        </p>
      </div>

      {results.fragmentReentryCapped ? (
        <div className="notice notice-amber">
          <AlertTriangle size={18} />
          <p>
            The {results.fragmentSequenceCap}-fragment browser safety cap was reached. Review later
            events cautiously because additional native fragment copies were not retained for re-screening.
          </p>
        </div>
      ) : null}

      {results.maxChiPeakLimitTriplets > 0 ? (
        <div className="notice notice-amber">
          <AlertTriangle size={18} />
          <p>
            MaxChi reached the supplied 100-peak retry bound for {results.maxChiPeakLimitTriplets.toLocaleString()}
            triplet{results.maxChiPeakLimitTriplets === 1 ? "" : "s"} while positive raw peaks
            remained. Calls are retained, but those profiles may contain additional unexplored peaks.
          </p>
        </div>
      ) : null}

      {results.chimaeraPeakLimitTargets > 0 ? (
        <div className="notice notice-amber">
          <AlertTriangle size={18} />
          <p>
            CHIMAERA reached the supplied 100-peak retry bound for {results.chimaeraPeakLimitTargets.toLocaleString()}
            target profile{results.chimaeraPeakLimitTargets === 1 ? "" : "s"} while positive raw
            peaks remained. Calls are retained, but those profiles may contain additional unexplored peaks.
          </p>
        </div>
      ) : null}

      {results.threeSeqApproximateEvaluations > 0 ? (
        <div className="notice notice-blue">
          <AlertTriangle size={18} />
          <p>
            {results.threeSeqApproximateEvaluations.toLocaleString()} 3SEQ orientation
            {results.threeSeqApproximateEvaluations === 1 ? " used" : "s used"} the supplied
            Siegmund discrete approximation because the exact hypergeometric state space exceeded
            the browser bound. Each anchor records which probability route produced it.
          </p>
        </div>
      ) : null}

      {!results.fragmentReentry ? (
        <div className="notice notice-amber">
          <AlertTriangle size={18} />
          <p>
            This alignment meets the supplied desktop source’s {results.fragmentReentryAlignmentLengthLimit.toLocaleString()}-site
            cutoff for suppressing synthetic fragment copies. Tracts were still erased between cyclic passes.
          </p>
        </div>
      ) : null}

      {pendingEvent !== null ? (
        <div className="notice notice-amber reconciliation-notice">
          <RefreshCw size={18} />
          <p>
            Event {pendingEvent + 1} was corrected or rejected. Record that decision, then
            re-identify the downstream chain before continuing the ordered review. Later rows are
            retained only as stale audit context until that rebuild replaces them.
          </p>
        </div>
      ) : null}

      {pendingEvent === null &&
      nextUnreviewedEvent !== null &&
      selected.id > nextUnreviewedEvent ? (
        <div className="notice notice-blue">
          <Layers3 size={18} />
          <p>
            Event {nextUnreviewedEvent + 1} is the next undecided event. You can inspect this call
            now, but record decisions in analysis order so later calls retain a valid history.
          </p>
        </div>
      ) : null}

      <div className="review-workspace">
        <aside className="event-list" aria-label="Reconciled recombination events">
          <div className="event-list-heading">
            <span className="eyebrow">Analysis order</span>
            <strong>{results.events.length} events</strong>
          </div>
          <div className="event-list-scroll">
            {results.events.map((event) => {
              const stale = pendingEvent !== null && event.id > pendingEvent;
              const recombinantReferenceGroup = event.queryReferenceInputRole === "reference"
                ? event.referenceGroup ?? 0
                : 0;
              return (
                <button
                  type="button"
                  key={event.id}
                  className={`event-list-item${selected.id === event.id ? " is-selected" : ""}${stale ? " is-stale" : ""}${recombinantReferenceGroup > 0 ? " is-reference-recombinant" : ""}`}
                  onClick={() => setSelectedId(event.id)}
                >
                  <span className={`event-status status-${event.reviewState}`} />
                  <span>
                    <strong>Event {event.id + 1}</strong>
                    <small title={event.recombinantName}>
                      {event.recombinantName}{recombinantReferenceGroup > 0 ? ` · ref ${recombinantReferenceGroup}` : ""}
                    </small>
                  </span>
                  <span>
                    <strong>{stale ? "stale" : pValue(event.bestCorrectedPValue)}</strong>
                    <small>
                      {stale
                        ? "awaiting rebuild"
                        : `${event.detectionMethods.join(" + ")} · ${event.supportSignalIds.length} signal${event.supportSignalIds.length === 1 ? "" : "s"}`}
                    </small>
                  </span>
                </button>
              );
            })}
          </div>
        </aside>

        <div className="event-detail">
          <div className="event-toolbar">
            <div>
              <span className="eyebrow">Event {selected.id + 1} of {results.events.length} · {selected.detectionMethods.join(" + ")}</span>
              <h2>{selected.recombinantName}</h2>
              {selectedRecombinantReferenceGroup > 0 ? (
                <span className="reference-recombinant-badge">
                  Reference recombinant · group {selectedRecombinantReferenceGroup}
                </span>
              ) : null}
            </div>
            <div className="event-navigation">
              <button type="button" onClick={() => move(-1)} aria-label="Previous event"><ChevronLeft /></button>
              <button type="button" onClick={() => move(1)} aria-label="Next event"><ChevronRight /></button>
            </div>
          </div>

          <EventSchematic event={selected} alignmentLength={alignmentLength} />

          {selectedRecombinantReferenceGroup > 0 ? (
            <div className="event-reference-warning">
              <AlertTriangle size={16} />
              <p>
                This event currently identifies a reference sequence as recombinant. That is an
                allowed and explicitly documented query-vs-reference outcome; inspect its parents
                and co-recombinant group rather than forcing the query into the recombinant role.
              </p>
            </div>
          ) : null}

          {selected.maxChiChimaeraOnlySupport ? (
            <div className="event-breakpoint-warning">
              <AlertTriangle size={16} />
              <p>
                This event is supported by MaxChi and CHIMAERA only. The RDP5 manual treats them as
                closely related methods, not independent confirmation; inspect the breakpoint,
                role, alignment, and tree evidence before accepting the call.
              </p>
            </div>
          ) : null}

          {beginningBreakpointUncertain || endingBreakpointUncertain ? (
            <div className="event-breakpoint-warning">
              <AlertTriangle size={16} />
              <p>
                {breakpointUncertaintyMessage} The RDP5 workflow treats the affected position as
                uncertain; open the original-alignment context before accepting this call.
              </p>
            </div>
          ) : null}

          <section className={`breakpoint-confidence-card${breakpointConfidence.available ? "" : " is-unavailable"}`}>
            <header>
              <div>
                <span className="eyebrow">Statistical breakpoint confidence</span>
                <h3>BURT / BenHMM polishing</h3>
              </div>
              <span className="fidelity-badge">
                {breakpointConfidence.available
                  ? `${breakpointConfidence.informationRichSites.toLocaleString()} information-rich sites`
                  : "No usable interval"}
              </span>
            </header>
            {breakpointConfidence.available ? (
              <>
                <div className="breakpoint-confidence-grid">
                  {breakpointConfidence.boundaries.map((boundary) => (
                    <article key={boundary.name} className={boundary.intervalAvailable ? "" : "is-unavailable"}>
                      <div>
                        <span>{boundary.name === "beginning" ? "Beginning boundary" : "Ending boundary"}</span>
                        <strong>
                          {boundary.inputCoordinate.toLocaleString()} → {boundary.polishedCoordinate.toLocaleString()}
                        </strong>
                        <small>
                          {boundary.intervalAvailable
                            ? boundary.sourceIntervalContainsInput
                              ? "input inside selected interval"
                              : "nearest interval did not contain input"
                            : "no matched state transition"}
                        </small>
                      </div>
                      <dl>
                        <div>
                          <dt>99% source range</dt>
                          <dd>
                            {boundary.intervalAvailable
                              ? `${confidenceCoordinate(boundary.confidence99.beginning)} → ${confidenceCoordinate(boundary.confidence99.ending)}${boundary.confidence99.wrapsOrigin ? " · wraps" : ""}`
                              : "—"}
                          </dd>
                        </div>
                        <div>
                          <dt>HMM position</dt>
                          <dd>{boundary.intervalAvailable ? confidenceCoordinate(boundary.hmmCoordinate) : "—"}</dd>
                        </div>
                        <div>
                          <dt>95% source range</dt>
                          <dd>
                            {boundary.intervalAvailable
                              ? `${confidenceCoordinate(boundary.confidence95.beginning)} → ${confidenceCoordinate(boundary.confidence95.ending)}${boundary.confidence95.wrapsOrigin ? " · wraps" : ""}`
                              : "—"}
                          </dd>
                        </div>
                      </dl>
                      <p>
                        {boundary.repositioned ? "Breakpoint repositioned" : "Breakpoint retained"}
                        {boundary.missingDataAdjusted ? " · input MissingData adjustment" : ""}
                        {boundary.finalGapAdjusted ? " · final gap relocation" : ""}
                      </p>
                    </article>
                  ))}
                </div>
                <footer>
                  <span>
                    {breakpointConfidence.candidateIntervalCount} candidate HMM interval{breakpointConfidence.candidateIntervalCount === 1 ? "" : "s"}
                    {breakpointConfidence.singleTransitionAssignment ? " · single-transition endpoint assignment" : ""}
                    {breakpointConfidence.insufficientInsideOrOutsideReverted ? " · coordinates reverted by the three-usable-site guard" : ""}
                  </span>
                  <span>
                    21 seeded starts · strict 0.995 / 0.999 posterior bounds · MSVC rand adapter
                  </span>
                </footer>
              </>
            ) : (
              <p className="breakpoint-confidence-empty">
                {breakpointConfidence.unavailableReason === "disabled"
                  ? "BURT polishing was disabled in the scan settings; the detected or manually edited coordinates were preserved."
                  : `BURT was attempted for this representative triplet but did not produce a matched interval (${breakpointConfidence.unavailableReason?.replaceAll("-", " ") ?? "unknown reason"}). The detected coordinates remain available for manual review.`}
              </p>
            )}
          </section>

          {maxChiDiscovery ? (
            <section className="maxchi-recheck-card">
              <header>
                <div>
                  <span className="eyebrow">Anchor discovery trace</span>
                  <h3>MaxChi peak-to-tract construction</h3>
                </div>
                <span className="maxchi-status is-hit">Corrected discovery hit</span>
              </header>
              <div className="maxchi-recheck-grid">
                <div>
                  <span>Raw peak</span>
                  <strong>{roleScore(maxChiDiscovery.maximumChiSquare)}</strong>
                  <small>{maxChiDiscoveryPairLabel(maxChiDiscovery.peakPair)}</small>
                </div>
                <div>
                  <span>Peak attempt</span>
                  <strong>{maxChiDiscovery.peakAttempt}</strong>
                  <small>Raw-χ² strongest-first heap order</small>
                </div>
                <div>
                  <span>Grown window</span>
                  <strong>{maxChiDiscovery.grownHalfWindow}</strong>
                  <small>Initial half-window {maxChiDiscovery.initialHalfWindow}</small>
                </div>
                <div>
                  <span>Selected tract</span>
                  <strong>{maxChiDiscovery.tractSide === "left" ? "Left" : "Right"}</strong>
                  <small>
                    flank χ² {roleScore(maxChiDiscovery.leftFlankChiSquare)} / {roleScore(maxChiDiscovery.rightFlankChiSquare)}
                  </small>
                </div>
                <div>
                  <span>Raw tail</span>
                  <strong>{pValue(maxChiDiscovery.rawPValue)}</strong>
                  <small>Before positional correction</small>
                </div>
                <div>
                  <span>Within triplet</span>
                  <strong>{pValue(maxChiDiscovery.withinTripletPValue)}</strong>
                  <small>{maxChiDiscovery.variableSites.toLocaleString()} positions × three profiles</small>
                </div>
                <div>
                  <span>Project corrected</span>
                  <strong>{pValue(maxChiDiscovery.correctedPValue)}</strong>
                  <small>Discovery threshold passed</small>
                </div>
                <div>
                  <span>Peak position</span>
                  <strong>{maxChiDiscovery.peakAlignmentPosition.toLocaleString()}</strong>
                  <small>Critical difference {maxChiDiscovery.criticalDifference}</small>
                </div>
              </div>
              <footer>
                <span>MCXoverF multi-peak discovery · source 12-term/11-divisor smoothing used for basin destruction only</span>
                <span>
                  {maxChiDiscovery.missingDataWindowFilterApplied ? "MissingData/erasure filter applied" : "No MissingData window bans"}
                  {maxChiDiscovery.linearEdgeWindowFilterApplied ? " · linear-edge filter applied" : ""}
                </span>
              </footer>
              <p className="maxchi-scope-note">
                This trace authored the anchor signal before ordinary event reconciliation. Final
                sequence roles and BURT coordinates may differ after consensus and polishing.
              </p>
            </section>
          ) : null}

          {chimaeraDiscovery && anchor && chimaeraParentOneLocal !== null ? (
            <section className="maxchi-recheck-card">
              <header>
                <div>
                  <span className="eyebrow">Anchor discovery trace</span>
                  <h3>CHIMAERA target profile and tract construction</h3>
                </div>
                <span className="maxchi-status is-hit">Corrected discovery hit</span>
              </header>
              <div className="maxchi-recheck-grid">
                <div>
                  <span>Candidate recombinant</span>
                  <strong>{anchor.tripletNames[chimaeraDiscovery.targetLocal]}</strong>
                  <small>Triplet member {chimaeraDiscovery.targetLocal + 1} of 3 rotations</small>
                </div>
                <div>
                  <span>Raw peak</span>
                  <strong>{roleScore(chimaeraDiscovery.maximumChiSquare)}</strong>
                  <small>
                    target / parent-one: {anchor.tripletNames[chimaeraDiscovery.targetLocal]} / {anchor.tripletNames[chimaeraParentOneLocal]}
                  </small>
                </div>
                <div>
                  <span>Peak attempt</span>
                  <strong>{chimaeraDiscovery.peakAttempt}</strong>
                  <small>Raw-χ² strongest-first order for this target</small>
                </div>
                <div>
                  <span>Grown window</span>
                  <strong>{chimaeraDiscovery.grownHalfWindow}</strong>
                  <small>Initial half-window {chimaeraDiscovery.initialHalfWindow}</small>
                </div>
                <div>
                  <span>Selected tract</span>
                  <strong>{chimaeraDiscovery.tractSide === "left" ? "Left" : "Right"}</strong>
                  <small>
                    flank χ² {roleScore(chimaeraDiscovery.leftFlankChiSquare)} / {roleScore(chimaeraDiscovery.rightFlankChiSquare)}
                  </small>
                </div>
                <div>
                  <span>Raw tail</span>
                  <strong>{pValue(chimaeraDiscovery.rawPValue)}</strong>
                  <small>Before positional correction</small>
                </div>
                <div>
                  <span>Within triplet</span>
                  <strong>{pValue(chimaeraDiscovery.withinTripletPValue)}</strong>
                  <small>{chimaeraDiscovery.informationRichSites.toLocaleString()} target positions × three rotations</small>
                </div>
                <div>
                  <span>Project corrected</span>
                  <strong>{pValue(chimaeraDiscovery.correctedPValue)}</strong>
                  <small>Discovery threshold passed</small>
                </div>
                <div>
                  <span>Peak position</span>
                  <strong>{chimaeraDiscovery.peakAlignmentPosition.toLocaleString()}</strong>
                  <small>Critical difference {chimaeraDiscovery.criticalDifference}</small>
                </div>
                <div>
                  <span>Parent-one contrast</span>
                  <strong>
                    {roleScore(chimaeraDiscovery.insideParentOneMatchRate)} / {roleScore(chimaeraDiscovery.outsideParentOneMatchRate)}
                  </strong>
                  <small>Inside / outside tract match rate</small>
                </div>
              </div>
              <footer>
                <span>AlistChi → FastRecCheckChim → CXoverA · source smoothing is used only for peak-basin destruction</span>
                <span>
                  {chimaeraDiscovery.missingDataWindowFilterApplied ? "MissingData/erasure filter applied" : "No MissingData window bans"}
                  {chimaeraDiscovery.linearEdgeWindowFilterApplied ? " · linear-edge filter applied" : ""}
                </span>
              </footer>
              <p className="maxchi-scope-note">
                This target-specific binary profile authored the anchor before shared role consensus.
                The retained target is provisional; final parent order and BURT coordinates may differ.
              </p>
            </section>
          ) : null}

          {geneconvDiscovery && anchor ? (
            <section className="maxchi-recheck-card">
              <header>
                <div>
                  <span className="eyebrow">Anchor discovery trace</span>
                  <h3>GENECONV six-track fragment score</h3>
                </div>
                <span className="maxchi-status is-hit">Corrected KA hit</span>
              </header>
              <div className="maxchi-recheck-grid">
                <div>
                  <span>Fragment track</span>
                  <strong>{geneconvDiscovery.track}</strong>
                  <small>
                    {geneconvDiscovery.track < 3 ? "Inner pair-match run" : "Outer discordant-sequence run"}
                  </small>
                </div>
                <div>
                  <span>Provisional roles</span>
                  <strong>{anchor.recombinantName}</strong>
                  <small>{anchor.minorParentName} minor · {anchor.majorParentName} major</small>
                </div>
                <div>
                  <span>Fragment score</span>
                  <strong>{geneconvDiscovery.fragmentScore.toLocaleString()}</strong>
                  <small>Strict critical score {geneconvDiscovery.criticalScore.toLocaleString()}</small>
                </div>
                <div>
                  <span>Signed evidence</span>
                  <strong>{geneconvDiscovery.positiveSites.toLocaleString()} / {geneconvDiscovery.discordantSites.toLocaleString()}</strong>
                  <small>Positive / discordant sites · penalty {geneconvDiscovery.mismatchPenalty}</small>
                </div>
                <div>
                  <span>Raw KA P</span>
                  <strong>{pValue(geneconvDiscovery.rawPValue)}</strong>
                  <small>GCCalcPValP2, before the initial scan-plan factor</small>
                </div>
                <div>
                  <span>Project corrected</span>
                  <strong>{pValue(geneconvDiscovery.correctedPValue)}</strong>
                  <small>
                    {results.correction === "bonferroni"
                      ? `RDP5 XOverList equivalent · × ${anchor.correctionTests.toLocaleString()}`
                      : "RDP5 XOverList equivalent · correction disabled"}
                  </small>
                </div>
                <div>
                  <span>Lambda</span>
                  <strong>{roleScore(geneconvDiscovery.lambda)}</strong>
                  <small>Bounded supplied Newton path</small>
                </div>
                <div>
                  <span>Karlin–Altschul K</span>
                  <strong>{roleScore(geneconvDiscovery.karlinAltschulK)}</strong>
                  <small>{geneconvDiscovery.polymorphicSites.toLocaleString()} non-monomorphic sites</small>
                </div>
              </div>
              <footer>
                <span>FindSubSeqGCAP6 → GetFragsP → GetMaxFragScoreP → CalcKMaxP → GCCalcPValP2 → GCXoverD</span>
                <span>Ignored indels · stable lowest-P ordering · overlap limit {results.geneconvMaxOverlaps}</span>
              </footer>
              <p className="maxchi-scope-note">
                The automated supplied path leaves the minimum fragment-length, polymorphism, and
                pair-score controls inactive. This KA fragment authored the anchor before shared
                role consensus and BURT breakpoint polishing. BURT can move the displayed event
                boundaries later, but it does not recalculate either probability shown above.
              </p>
            </section>
          ) : null}

          {bootscanDiscovery && anchor ? (
            <section className="maxchi-recheck-card">
              <header>
                <div>
                  <span className="eyebrow">Anchor discovery trace</span>
                  <h3>BootScan closest-pair bootstrap support</h3>
                </div>
                <span className="maxchi-status is-hit">Corrected binomial hit</span>
              </header>
              <div className="maxchi-recheck-grid">
                <div>
                  <span>Supported pair</span>
                  <strong>{maxChiDiscoveryPairLabel(bootscanDiscovery.supportedPair)}</strong>
                  <small>Strict unique closest-pair votes only</small>
                </div>
                <div>
                  <span>Provisional roles</span>
                  <strong>{anchor.recombinantName}</strong>
                  <small>{anchor.minorParentName} minor · {anchor.majorParentName} major</small>
                </div>
                <div>
                  <span>Peak / mean support</span>
                  <strong>{(bootscanDiscovery.maximumPairSupport * 100).toFixed(1)}% / {(bootscanDiscovery.meanPairSupport * 100).toFixed(1)}%</strong>
                  <small>{bootscanDiscovery.usableWindows} usable of {bootscanDiscovery.windowsScored} tract windows</small>
                </div>
                <div>
                  <span>Bootstrap trace P</span>
                  <strong>{pValue(bootscanDiscovery.bootstrapPValue)}</strong>
                  <small>1 − mean supported-pair fraction</small>
                </div>
                <div>
                  <span>Raw binomial P</span>
                  <strong>{pValue(bootscanDiscovery.rawPValue)}</strong>
                  <small>{bootscanDiscovery.tractPairMatches} matches in {bootscanDiscovery.tractInformativeSites} tract sites</small>
                </div>
                <div>
                  <span>Project corrected</span>
                  <strong>{pValue(bootscanDiscovery.correctedPValue)}</strong>
                  <small>
                    {results.correction === "bonferroni"
                      ? `MakeScoresBS equivalent · × ${anchor.correctionTests.toLocaleString()}`
                      : "MakeScoresBS equivalent · correction disabled"}
                  </small>
                </div>
              </div>
              <footer>
                <span>BSXoverR → SEQBOOT2 → FastBootDist → GetPltVal → ScanBSPlots → MakeBSEvent</span>
                <span>
                  JC distance mode · {results.bootscanWindowSites}/{results.bootscanStepSites} sites · {results.bootscanBootstrapReplicates} replicates
                </span>
              </footer>
              <p className="maxchi-scope-note">
                This is the supplied automated distance/binomial path. Tree, alternative
                substitution-model, and permutation modes are not represented by this trace.
                The binomial P is fixed from the detected pre-BURT bounds and is not recalculated
                after polishing. Its log-domain tail avoids factorial underflow, so an extremely
                small positive value can survive where desktop arithmetic becomes zero. Final roles
                and coordinates may differ after shared consensus and BURT polishing.
              </p>
            </section>
          ) : null}

          {siscanDiscovery && anchor ? (
            <section className="maxchi-recheck-card">
              <header>
                <div>
                  <span className="eyebrow">Anchor discovery trace</span>
                  <h3>SISCAN sister-pair permutation switch</h3>
                </div>
                <span className="maxchi-status is-hit">
                  {results.correction === "bonferroni" ? "Corrected SISCAN hit" : "Uncorrected SISCAN hit"}
                </span>
              </header>
              <div className="maxchi-recheck-grid">
                <div>
                  <span>Background sister pair</span>
                  <strong>{maxChiDiscoveryPairLabel(siscanDiscovery.globalPair)}</strong>
                  <small>Closest pair across the full triplet</small>
                </div>
                <div>
                  <span>Tract sister pair</span>
                  <strong>{maxChiDiscoveryPairLabel(siscanDiscovery.candidatePair)}</strong>
                  <small>{siscanDiscovery.windowsInRegion} consecutive pair-switch window{siscanDiscovery.windowsInRegion === 1 ? "" : "s"}</small>
                </div>
                <div>
                  <span>Nearest outlier</span>
                  <strong>{sequences[siscanDiscovery.outlierSequence]?.name ?? `Sequence ${siscanDiscovery.outlierSequence + 1}`}</strong>
                  <small>GetSSOL on the round-cached source WPGMA context</small>
                </div>
                <div>
                  <span>Maximum Z</span>
                  <strong>{roleScore(siscanDiscovery.maximumZ)}</strong>
                  <small>{siscanDiscovery.selectedScoreFamily} score {siscanDiscovery.selectedScore}</small>
                </div>
                <div>
                  <span>Normal tail P</span>
                  <strong>{pValue(siscanDiscovery.normalTailPValue)}</strong>
                  <small>Supplied NormalZ approximation</small>
                </div>
                <div>
                  <span>Region / window adjusted</span>
                  <strong>{pValue(siscanDiscovery.regionLengthAdjustedPValue)} / {pValue(siscanDiscovery.windowAdjustedPValue)}</strong>
                  <small>Shrunken length, then alignment/window opportunity factor</small>
                </div>
                <div>
                  <span>Project corrected</span>
                  <strong>{pValue(siscanDiscovery.correctedPValue)}</strong>
                  <small>{results.correction === "bonferroni" ? `${anchor.correctionTests.toLocaleString()} initial scan-plan opportunities` : "Project correction disabled"}</small>
                </div>
                <div>
                  <span>Permutation work</span>
                  <strong>{siscanDiscovery.permutationDraws.toLocaleString()}</strong>
                  <small>{siscanDiscovery.informativeSites.toLocaleString()} retained variable-pattern sites</small>
                </div>
              </div>
              <footer>
                <span>SSXoverC → GetSSOL → Get3Score/GetPScores2 → DoPerms3 → MakeZValue2/DoSums → FindMaxZ → ShrinkRegionC</span>
                <span>{results.siscanWindowSites}/{results.siscanStepSites} sites · {results.siscanScanPermutations}/{results.siscanPValuePermutations} scan/final permutations · seed {results.siscanRandomSeed}</span>
              </footer>
              <p className="maxchi-scope-note">
                The browser retains the supplied QuickCheckB control-flow quirk and one seeded
                MakeVRand-style flat random prefix. The WPGMA context and random prefix are shared
                across triplets; prior-event erasure invalidates only the state-dependent context.
                The displayed P-value stages remain separate so native differences can be localized.
              </p>
            </section>
          ) : null}

          {threeSeqDiscovery && anchor ? (
            <section className="maxchi-recheck-card">
              <header>
                <div>
                  <span className="eyebrow">Anchor discovery trace</span>
                  <h3>3SEQ hypergeometric random-walk excursion</h3>
                </div>
                <span className="maxchi-status is-hit">
                  {results.correction === "bonferroni" ? "Corrected 3SEQ hit" : "Uncorrected 3SEQ hit"}
                </span>
              </header>
              <div className="maxchi-recheck-grid">
                <div>
                  <span>Candidate recombinant</span>
                  <strong>{anchor.tripletNames[threeSeqDiscovery.targetLocal]}</strong>
                  <small>Triplet member {threeSeqDiscovery.targetLocal + 1} of 3 rotations</small>
                </div>
                <div>
                  <span>Selected walk</span>
                  <strong>{threeSeqDiscovery.walkDirection === "ascent" ? "Ascent" : "Descent"}</strong>
                  <small>Strict lower-P orientation; equal tails retain descent</small>
                </div>
                <div>
                  <span>Boundary excursion</span>
                  <strong>{threeSeqDiscovery.maximumExcursion.toLocaleString()}</strong>
                  <small>After source CheckwrapC tract extension</small>
                </div>
                <div>
                  <span>Probability excursion</span>
                  <strong>{threeSeqDiscovery.probabilityExcursion.toLocaleString()}</strong>
                  <small>
                    {threeSeqDiscovery.missingDataSplitApplied
                      ? "CheckSplit3Seq sub-tract nK used by GetTSPVal"
                      : "Pre-wrap nK used by Seq3PVals/GetTSPVal"}
                  </small>
                </div>
                <div>
                  <span>Parent-match balance</span>
                  <strong>
                    {threeSeqDiscovery.parentOneMatches.toLocaleString()} / {threeSeqDiscovery.parentTwoMatches.toLocaleString()}
                  </strong>
                  <small>Selected parent one / parent two</small>
                </div>
                <div>
                  <span>Raw hypergeometric P</span>
                  <strong>{pValue(threeSeqDiscovery.rawPValue)}</strong>
                  <small>Before the project opportunity factor</small>
                </div>
                <div>
                  <span>{results.correction === "bonferroni" ? "Dunn–Šidák corrected" : "Threshold P"}</span>
                  <strong>{pValue(threeSeqDiscovery.correctedPValue)}</strong>
                  <small>
                    {results.correction === "bonferroni"
                      ? `${anchor.correctionTests.toLocaleString()} project opportunities`
                      : "Project multiple-comparison correction disabled"}
                  </small>
                </div>
                <div>
                  <span>Probability route</span>
                  <strong>{threeSeqDiscovery.exactProbability ? "Exact" : "Siegmund"}</strong>
                  <small>
                    {threeSeqDiscovery.siegmundFallback
                      ? "Supplied large-profile fallback"
                      : "Finite hypergeometric walk distribution"}
                  </small>
                </div>
                <div>
                  <span>Information-rich sites</span>
                  <strong>{threeSeqDiscovery.informationRichSites.toLocaleString()}</strong>
                  <small>Parents differ; target matches exactly one</small>
                </div>
              </div>
              <footer>
                <span>FindSubSeqTS → Seq3PVals/GetTSPVal → CheckwrapC → TSXOver</span>
                <span>
                  {threeSeqDiscovery.missingDataSplitApplied
                    ? "Post-erasure CheckSplit3Seq trim and re-probability applied"
                    : "No prior-erasure missing run changed this tract"}
                </span>
              </footer>
              <p className="maxchi-scope-note">
                This target-specific walk authored the anchor before shared role consensus and
                BURT polishing. The bounded exact evaluator replaces the supplied desktop lookup
                table without rescanning alignment bytes; native saved-output validation is still
                required.
              </p>
            </section>
          ) : null}

          <section className={`maxchi-recheck-card${geneconvRecheck.profileAvailable ? "" : " is-unavailable"}`}>
            <header>
              <div>
                <span className="eyebrow">Six signed tracks · secondary corroboration</span>
                <h3>GENECONV recheck</h3>
              </div>
              <span className={`maxchi-status${geneconvRecheck.sourceRecheckHit ? " is-hit" : ""}`}>
                {geneconvRecheck.profileAvailable
                  ? geneconvRecheck.bestTrack === null
                    ? "No qualifying fragment"
                    : geneconvRecheck.sourceRecheckHit ? "Corrected hit" : "No corrected hit"
                  : "Profile unavailable"}
              </span>
            </header>
            {geneconvRecheck.profileAvailable ? (
              <>
                <div className="maxchi-recheck-grid">
                  <div>
                    <span>Best track</span>
                    <strong>{geneconvRecheck.bestTrack ?? "—"}</strong>
                    <small>{geneconvRecheck.bestTrack === null ? "No surviving fragment" : geneconvRecheck.bestTrack < 3 ? "Inner pair-match" : "Outer discordant-sequence"}</small>
                  </div>
                  <div>
                    <span>Provisional recombinant</span>
                    <strong>{geneconvRecheckRecombinantName ?? "—"}</strong>
                    <small>Representative-triplet role only</small>
                  </div>
                  <div>
                    <span>Raw KA P</span>
                    <strong>{geneconvRecheck.rawPValue === null ? "—" : pValue(geneconvRecheck.rawPValue)}</strong>
                    <small>Before the initial scan-plan factor</small>
                  </div>
                  <div>
                    <span>Project corrected</span>
                    <strong>{geneconvRecheck.correctedPValue === null ? "—" : pValue(geneconvRecheck.correctedPValue)}</strong>
                    <small>{geneconvRecheck.bonferroniApplied ? `${geneconvRecheck.correctionTests.toLocaleString()} initial scan-plan opportunities` : "Correction disabled"}</small>
                  </div>
                  <div>
                    <span>Fragment score</span>
                    <strong>{geneconvRecheck.fragmentScore ?? "—"}</strong>
                    <small>Strict critical score {geneconvRecheck.criticalScore ?? "—"}</small>
                  </div>
                  <div>
                    <span>Candidate tract</span>
                    <strong>{geneconvRecheck.beginning === null || geneconvRecheck.ending === null ? "—" : `${geneconvRecheck.beginning}–${geneconvRecheck.ending}`}</strong>
                    <small>{geneconvRecheck.wrapsOrigin ? "Origin-wrapping" : "Linear coordinates"}</small>
                  </div>
                </div>
                <footer>
                  <span>
                    {geneconvRecheck.polymorphicSites.toLocaleString()} polymorphic sites · {geneconvRecheck.tracksScreened} tracks · {geneconvRecheck.fragmentsScored.toLocaleString()} positive starts
                  </span>
                  <span>
                    {geneconvRecheck.qualifiedFragments.toLocaleString()} above critical · {geneconvRecheck.overlapRejectedFragments.toLocaleString()} overlap rejection{geneconvRecheck.overlapRejectedFragments === 1 ? "" : "s"}
                    {geneconvRecheck.numericalFallbackTracks ? ` · ${geneconvRecheck.numericalFallbackTracks} bounded numerical fallback(s)` : ""}
                  </span>
                </footer>
              </>
            ) : (
              <p>
                The representative triplet did not retain a usable ordinary ignored-indel
                GENECONV profile after gaps, prior erasures, all-one-pair, and source skew gates.
              </p>
            )}
            <p className="maxchi-scope-note">
              This reuses the supplied GCXoverD six-track KA kernel as a non-coordinate-changing
              representative check. It is separate from event discovery; permutation, manual-pair,
              alternative-indel, and full native late event-reconstruction modes remain pending.
            </p>
          </section>

          <section className={`maxchi-recheck-card${threeSeqRecheck.profileAvailable ? "" : " is-unavailable"}`}>
            <header>
              <div>
                <span className="eyebrow">Two orientations · inverse interval copies · secondary corroboration</span>
                <h3>3SEQ Findall recheck</h3>
              </div>
              <span className={`maxchi-status${threeSeqRecheck.sourceRecheckHit ? " is-hit" : ""}`}>
                {threeSeqRecheck.profileAvailable
                  ? threeSeqRecheck.sourceRecheckHit ? "Corrected hit" : "No corrected hit"
                  : "Profile unavailable"}
              </span>
            </header>
            {threeSeqRecheck.profileAvailable ? (
              <>
                <div className="maxchi-recheck-grid">
                  <div>
                    <span>Best target</span>
                    <strong>{threeSeqRecheckTargetName ?? "—"}</strong>
                    <small>{threeSeqRecheck.bestDirection ?? "No qualifying orientation"}</small>
                  </div>
                  <div>
                    <span>Raw hypergeometric P</span>
                    <strong>{threeSeqRecheck.rawPValue === null ? "—" : pValue(threeSeqRecheck.rawPValue)}</strong>
                    <small>{threeSeqRecheck.exactProbability ? "Exact Single-state DP" : threeSeqRecheck.siegmundFallback ? "Siegmund fallback" : "No qualifying call"}</small>
                  </div>
                  <div>
                    <span>Project corrected</span>
                    <strong>{threeSeqRecheck.correctedPValue === null ? "—" : pValue(threeSeqRecheck.correctedPValue)}</strong>
                    <small>{threeSeqRecheck.correctionApplied ? `${threeSeqRecheck.correctionTests.toLocaleString()} initial scan-plan opportunities` : "Correction disabled"}</small>
                  </div>
                  <div>
                    <span>Candidate tract</span>
                    <strong>{threeSeqRecheck.beginning === null || threeSeqRecheck.ending === null ? "—" : `${threeSeqRecheck.beginning}–${threeSeqRecheck.ending}`}</strong>
                    <small>{threeSeqRecheck.wrapsOrigin ? "Origin-wrapping" : "Linear coordinates"}</small>
                  </div>
                  <div>
                    <span>Walk evidence</span>
                    <strong>{threeSeqRecheck.parentOneMatches} / {threeSeqRecheck.parentTwoMatches}</strong>
                    <small>{threeSeqRecheck.informationRichSites.toLocaleString()} information-rich sites</small>
                  </div>
                  <div>
                    <span>Findall output</span>
                    <strong>{threeSeqRecheck.qualifyingOrientations} / {threeSeqRecheck.sourceListEntries}</strong>
                    <small>Qualifying orientations / source list entries</small>
                  </div>
                </div>
                <footer>
                  <span>{threeSeqRecheck.targetProfilesScanned} target profiles · {threeSeqRecheck.exactProbabilityEvaluations} exact evaluations · {threeSeqRecheck.approximateProbabilityEvaluations} fallback evaluations</span>
                  <span>{threeSeqRecheck.missingDataSplitApplied ? "CheckSplit3Seq/SubPVal changed the selected evidence" : "No selected split evidence"}</span>
                </footer>
              </>
            ) : (
              <p>No target rotation retained the supplied four information-rich sites after gaps and prior erasures.</p>
            )}
            <p className="maxchi-scope-note">
              This is the supplied TSXOver(1) two-orientation Findall shape, including the
              inverse-interval list copy for each qualifying call. It corroborates the representative
              triplet without moving reconciled coordinates and remains native-golden unvalidated.
            </p>
          </section>

          <section className={`maxchi-recheck-card${siscanRecheck.profileAvailable ? "" : " is-unavailable"}`}>
            <header>
              <div>
                <span className="eyebrow">Nearest fourth sequence · fixed event bounds · secondary corroboration</span>
                <h3>SISCAN fixed-region recheck</h3>
              </div>
              <span className={`maxchi-status${siscanRecheck.sourceRecheckHit ? " is-hit" : ""}`}>
                {siscanRecheck.representativeSkipped
                  ? "Representative skipped"
                  : siscanRecheck.profileAvailable
                    ? siscanRecheck.sourceRecheckHit ? "Corrected hit" : "No corrected hit"
                    : siscanRecheck.status === "not-requested" ? "Not requested" : "Profile unavailable"}
              </span>
            </header>
            {siscanRecheck.profileAvailable ? (
              <>
                <div className="maxchi-recheck-grid">
                  <div>
                    <span>Nearest outlier</span>
                    <strong>{siscanRecheck.outlierSequence === null ? "—" : sequences[siscanRecheck.outlierSequence]?.name ?? `Sequence ${siscanRecheck.outlierSequence + 1}`}</strong>
                    <small>Source WPGMA cophenetic selection</small>
                  </div>
                  <div>
                    <span>Scored sister pair</span>
                    <strong>{maxChiDiscoveryPairLabel(siscanRecheck.scoredPair)}</strong>
                    <small>Background pair {maxChiDiscoveryPairLabel(siscanRecheck.globalPair)}</small>
                  </div>
                  <div>
                    <span>Maximum Z</span>
                    <strong>{roleScore(siscanRecheck.maximumZ)}</strong>
                    <small>{siscanRecheck.selectedScoreFamily} score {siscanRecheck.selectedScore}</small>
                  </div>
                  <div>
                    <span>Normal / region P</span>
                    <strong>{pValue(siscanRecheck.normalTailPValue)} / {pValue(siscanRecheck.regionLengthAdjustedPValue)}</strong>
                    <small>Before alignment/window and project factors</small>
                  </div>
                  <div>
                    <span>Window-adjusted P</span>
                    <strong>{pValue(siscanRecheck.windowAdjustedPValue)}</strong>
                    <small>Fixed reconciled event region</small>
                  </div>
                  <div>
                    <span>Project corrected</span>
                    <strong>{pValue(siscanRecheck.correctedPValue)}</strong>
                    <small>{siscanRecheck.bonferroniApplied ? `${siscanRecheck.correctionTests.toLocaleString()} initial scan-plan opportunities` : "Correction disabled"}</small>
                  </div>
                </div>
                <footer>
                  <span>{siscanRecheck.informativeSites.toLocaleString()} retained sites · {siscanRecheck.permutationDraws.toLocaleString()} vertical-permutation draws</span>
                  <span>GetSSOL → Get3Score/GetPScores2 → DoPerms3P → MakeZValue2 → DoSums</span>
                </footer>
              </>
            ) : (
              <p>
                {siscanRecheck.representativeSkipped
                  ? "This finalized-list row is the event representative; its triplet result is shown once in the representative evidence."
                  : siscanRecheck.status === "not-requested"
                    ? "SISCAN confirmation was disabled for this analysis."
                    : "No eligible distinct-origin fourth sequence or usable fixed-region variable-pattern profile remained."}
              </p>
            )}
            <p className="maxchi-scope-note">
              This source-shaped confirmation scores the already reconciled bounds and cannot move
              the event. It is enabled by default to match the ordinary RDP5 confirmation workflow;
              native saved-output comparison is still required.
            </p>
          </section>

          <section className={`maxchi-recheck-card${maxChiRecheck.profileAvailable ? "" : " is-unavailable"}`}>
            <header>
              <div>
                <span className="eyebrow">Secondary triplet corroboration</span>
                <h3>MaxChi recheck</h3>
              </div>
              <span className={`maxchi-status${maxChiRecheck.sourceRecheckHit ? " is-hit" : ""}`}>
                {maxChiRecheck.profileAvailable
                  ? maxChiRecheck.bestPair === null
                    ? "No qualifying peak"
                    : maxChiRecheck.sourceRecheckHit ? "Corrected hit" : "No corrected hit"
                  : "Profile unavailable"}
              </span>
            </header>
            {maxChiRecheck.profileAvailable ? (
              <>
                <div className="maxchi-recheck-grid">
                  <div>
                    <span>Maximum χ²</span>
                    <strong>{maxChiRecheck.maximumChiSquare === null ? "—" : roleScore(maxChiRecheck.maximumChiSquare)}</strong>
                    <small>{maxChiPairLabel(maxChiRecheck.bestPair)}</small>
                  </div>
                  <div>
                    <span>Raw χ² tail</span>
                    <strong>{maxChiRecheck.localPValue === null ? "—" : pValue(maxChiRecheck.localPValue)}</strong>
                    <small>Before position or triplet correction</small>
                  </div>
                  <div>
                    <span>Within triplet</span>
                    <strong>{maxChiRecheck.withinTripletPValue === null ? "—" : pValue(maxChiRecheck.withinTripletPValue)}</strong>
                    <small>Variable positions × three pair profiles</small>
                  </div>
                  <div>
                    <span>Project corrected</span>
                    <strong>{maxChiRecheck.correctedPValue === null ? "—" : pValue(maxChiRecheck.correctedPValue)}</strong>
                    <small>{maxChiRecheck.bonferroniApplied ? `${maxChiRecheck.correctionTests.toLocaleString()} initial scan-plan opportunities` : "Correction disabled"}</small>
                  </div>
                </div>
                <footer>
                  <span>
                    {maxChiRecheck.variableSites.toLocaleString()} MaxChi-variable sites · initial half-window {maxChiRecheck.halfWindow} · grown {maxChiRecheck.bestPair === null ? "—" : maxChiRecheck.grownHalfWindow}
                  </span>
                  <span>
                    Peak at alignment position {maxChiRecheck.peakAlignmentPosition ?? "—"} · critical match difference {maxChiRecheck.criticalDifference}
                    {maxChiRecheck.missingDataWindowFilterApplied ? " · MissingData/erasure filter applied" : ""}
                    {maxChiRecheck.linearEdgeWindowFilterApplied ? " · linear-edge filter applied" : ""}
                  </span>
                </footer>
              </>
            ) : (
              <p>
                The representative triplet did not retain the seven variable sites and six-site
                half-window required by the supplied MaxChi path after gaps, missing data, and
                earlier erased tracts were applied.
              </p>
            )}
            <p className="maxchi-scope-note">
              This is the source-shaped FastRecCheckMC2 strongest-peak statistic. It corroborates
              the finalized event hypothesis without moving its coordinates; it is distinct from
              the active MCXoverF MaxChi discovery path shown above when MaxChi supplied the anchor.
            </p>
          </section>

          <section className={`maxchi-recheck-card${chimaeraRecheck.profileAvailable ? "" : " is-unavailable"}`}>
            <header>
              <div>
                <span className="eyebrow">Three target rotations · secondary corroboration</span>
                <h3>CHIMAERA recheck</h3>
              </div>
              <span className={`maxchi-status${chimaeraRecheck.sourceRecheckHit ? " is-hit" : ""}`}>
                {chimaeraRecheck.profileAvailable
                  ? chimaeraRecheck.bestTarget === null
                    ? "No qualifying peak"
                    : chimaeraRecheck.sourceRecheckHit ? "Corrected hit" : "No corrected hit"
                  : "Profile unavailable"}
              </span>
            </header>
            {chimaeraRecheck.profileAvailable ? (
              <>
                <div className="maxchi-recheck-grid">
                  <div>
                    <span>Maximum χ²</span>
                    <strong>{chimaeraRecheck.maximumChiSquare === null ? "—" : roleScore(chimaeraRecheck.maximumChiSquare)}</strong>
                    <small>{chimaeraRecheckTargetName ? `Target ${chimaeraRecheckTargetName}` : "No screened target peak"}</small>
                  </div>
                  <div>
                    <span>Raw χ² tail</span>
                    <strong>{chimaeraRecheck.localPValue === null ? "—" : pValue(chimaeraRecheck.localPValue)}</strong>
                    <small>Before position or triplet correction</small>
                  </div>
                  <div>
                    <span>Within triplet</span>
                    <strong>{chimaeraRecheck.withinTripletPValue === null ? "—" : pValue(chimaeraRecheck.withinTripletPValue)}</strong>
                    <small>Target positions × three rotations</small>
                  </div>
                  <div>
                    <span>Project corrected</span>
                    <strong>{chimaeraRecheck.correctedPValue === null ? "—" : pValue(chimaeraRecheck.correctedPValue)}</strong>
                    <small>{chimaeraRecheck.bonferroniApplied ? `${chimaeraRecheck.correctionTests.toLocaleString()} initial scan-plan opportunities` : "Correction disabled"}</small>
                  </div>
                </div>
                <footer>
                  <span>
                    {chimaeraRecheck.targetProfilesScanned} target profile{chimaeraRecheck.targetProfilesScanned === 1 ? "" : "s"} · {chimaeraRecheck.informationRichSites.toLocaleString()} selected-target information-rich sites · initial half-window {chimaeraRecheck.halfWindow} · grown {chimaeraRecheck.bestTarget === null ? "—" : chimaeraRecheck.grownHalfWindow}
                  </span>
                  <span>
                    Peak at alignment position {chimaeraRecheck.peakAlignmentPosition ?? "—"} · critical match difference {chimaeraRecheck.criticalDifference}
                    {chimaeraRecheck.missingDataWindowFilterApplied ? " · MissingData/erasure filter applied" : ""}
                    {chimaeraRecheck.linearEdgeWindowFilterApplied ? " · linear-edge filter applied" : ""}
                  </span>
                </footer>
              </>
            ) : (
              <p>
                None of the three target rotations retained enough information-rich positions for
                the supplied CHIMAERA window and critical-difference screen.
              </p>
            )}
            <p className="maxchi-scope-note">
              This is the source-shaped FastRecCheckChim three-target strongest-peak statistic. It
              corroborates the finalized event hypothesis without moving its coordinates and is
              separate from exploratory CXoverA event discovery. The manual treats CHIMAERA and
              MaxChi as closely related methods, not independent confirmation.
            </p>
          </section>

          <div className="role-grid">
            <div className="role-recombinant">
              <span>Current recombinant</span>
              <strong>{selected.recombinantName}</strong>
              {inputRoleLabel(selected.recombinant) ? <small>{inputRoleLabel(selected.recombinant)}</small> : null}
            </div>
            <div className="role-major">
              <span>Major-parent-like</span>
              <strong>{selected.majorParentName}</strong>
              {inputRoleLabel(selected.majorParent) ? <small>{inputRoleLabel(selected.majorParent)}</small> : null}
            </div>
            <div className="role-minor">
              <span>Minor-parent-like</span>
              <strong>{selected.minorParentName}</strong>
              {inputRoleLabel(selected.minorParent) ? <small>{inputRoleLabel(selected.minorParent)}</small> : null}
            </div>
            <div>
              <span>Best corrected p-value</span>
              <strong>{pValue(selected.bestCorrectedPValue)}</strong>
            </div>
          </div>

          <div className="evidence-strip">
            <span><strong>{selected.detectionRound}</strong> detection round</span>
            <span><strong>{selected.erasedNucleotideSites.toLocaleString()}</strong> sites erased</span>
            <span><strong>{selected.fragmentSequencesAdded}</strong> fragments re-entered</span>
            {selected.fragmentAssistedDetection ? <span className="manual-badge">Fragment-assisted</span> : null}
            <span><Layers3 size={14} /><strong>{selected.supportSignalIds.length}</strong> discovery signals</span>
            <span><strong>{currentHypothesis.detectableSignalSetIndices.length}</strong> detectable</span>
            <span><strong>{currentHypothesis.distanceCorrelationSetIndices.length}</strong> distance-correlated</span>
            <span><strong>{currentHypothesis.phylogeneticCorrelationSetIndices.length}</strong> tree-correlated</span>
            <span><strong>{selected.coRecombinantSequenceIndices.length}</strong> co-recombinant</span>
            <span><strong>{selected.traceEvidence.length}</strong> masked traces</span>
            <span>
              <GitCompareArrows size={14} />
              <strong>{breakpointConfidence.candidateIntervalCount}</strong> BURT intervals
            </span>
            <span className={maxChiRecheck.sourceRecheckHit ? "maxchi-hit-badge" : undefined}>
              <strong>MaxChi</strong> {maxChiRecheck.profileAvailable
                ? maxChiRecheck.bestPair === null
                  ? "no qualifying peak"
                  : maxChiRecheck.sourceRecheckHit ? "corrected hit" : "no corrected hit"
                : "unavailable"}
            </span>
            <span className={chimaeraRecheck.sourceRecheckHit ? "maxchi-hit-badge" : undefined}>
              <strong>CHIMAERA</strong> {chimaeraRecheck.profileAvailable
                ? chimaeraRecheck.bestTarget === null
                  ? "no qualifying peak"
                  : chimaeraRecheck.sourceRecheckHit ? "corrected hit" : "no corrected hit"
                : "unavailable"}
            </span>
            <span className={threeSeqRecheck.sourceRecheckHit ? "maxchi-hit-badge" : undefined}>
              <strong>3SEQ</strong> {threeSeqRecheck.profileAvailable
                ? threeSeqRecheck.sourceRecheckHit ? "corrected hit" : "no corrected hit"
                : "unavailable"}
            </span>
            {beginningBreakpointUncertain ? (
              <span className="uncertainty-badge">Beginning boundary uncertain</span>
            ) : null}
            {endingBreakpointUncertain ? (
              <span className="uncertainty-badge">Ending boundary uncertain</span>
            ) : null}
            {selected.groupManualAdjusted ? <span className="manual-badge">Manual group</span> : null}
            {selected.manualAdjusted ? <span className="manual-badge">Manually corrected</span> : null}
          </div>

          <section className="role-consensus-card">
            <div className="role-consensus-heading">
              <span className="role-consensus-icon"><GitBranch size={18} /></span>
              <div>
                <span className="eyebrow">Recombinant identification</span>
                <h3>
                  {selected.roleConsensus.informative
                    ? `${selected.roleConsensus.recommendedRecombinantName} is the weighted recommendation`
                    : "The available role metrics are not informative"}
                </h3>
                <p>
                  {selected.roleConsensus.informative
                    ? `${(selected.roleConsensus.confidence * 100).toFixed(0)}% decision-score margin · major-parent-like ${selected.roleConsensus.recommendedMajorParentName} · minor-parent-like ${selected.roleConsensus.recommendedMinorParentName}`
                    : "Keep the current roles and inspect the profile manually."}
                </p>
              </div>
              {selected.roleConsensus.informative && recommendationDiffers ? (
                <button
                  className="button button-secondary"
                  type="button"
                  onClick={applyRecommendation}
                  disabled={decisionBlocked}
                >
                  Apply recommendation
                </button>
              ) : selected.roleConsensus.informative ? (
                <span className="consensus-aligned"><Check size={14} /> Current roles agree</span>
              ) : null}
            </div>
            <div className="role-metric-row">
              {selected.roleConsensus.metrics.map((metric) => (
                <span className={metric.informative ? "is-informative" : ""} key={metric.method}>
                  <b>
                    {metric.method}
                    <i>{metric.weight > 0 ? `${metric.weight}-point` : "context"}</i>
                  </b>
                  <strong>
                    {metric.winningRole === null
                      ? "not decisive"
                      : selected.roleHypotheses[metric.winningRole]?.presumedRecombinantName ?? "—"}
                  </strong>
                  <small>{metric.scores.map(roleScore).join(" · ")}</small>
                </span>
              ))}
            </div>
          </section>

          {editing ? (
            <div className="event-editor">
              <div className="card-heading split-heading">
                <div>
                  <span className="eyebrow">Manual correction</span>
                  <h3>Repair roles or breakpoints</h3>
                </div>
                <button className="button button-quiet" type="button" onClick={() => setEditing(false)}>Cancel</button>
              </div>
              <div className="editor-grid">
                {([
                  ["Recombinant", "recombinant"],
                  ["Major-parent-like", "majorParent"],
                  ["Minor-parent-like", "minorParent"],
                ] as const).map(([label, key]) => (
                  <label key={key}>
                    <span>{label}</span>
                    <select
                      value={draft[key]}
                      onChange={(event) => setDraft((current) => ({ ...current, [key]: Number(event.target.value) }))}
                    >
                      {roleChoices.map((sequence) => (
                        <option value={sequence.index} key={sequence.index}>
                          {sequence.name}{inputRoleLabel(sequence.index) ? ` — ${inputRoleLabel(sequence.index)}` : ""}
                        </option>
                      ))}
                    </select>
                  </label>
                ))}
                <label>
                  <span>Beginning</span>
                  <input
                    type="number"
                    min={1}
                    max={alignmentLength}
                    value={draft.beginning}
                    onChange={(event) => setDraft((current) => ({ ...current, beginning: Number(event.target.value) }))}
                  />
                </label>
                <label>
                  <span>Ending</span>
                  <input
                    type="number"
                    min={1}
                    max={alignmentLength}
                    value={draft.ending}
                    onChange={(event) => setDraft((current) => ({ ...current, ending: Number(event.target.value) }))}
                  />
                </label>
              </div>
              {editError ? <p className="editor-error" role="alert">{editError}</p> : null}
              <button className="button button-primary" type="button" onClick={saveCorrection}>
                <Save size={16} /> Save correction
              </button>
            </div>
          ) : (
            <button
              className="edit-event-button"
              type="button"
              onClick={() => setEditing(true)}
              disabled={decisionBlocked}
              title={decisionBlocked ? "Finish the earlier workflow decision first" : undefined}
            >
              <Pencil size={15} /> Correct roles or breakpoints
            </button>
          )}

          <div className="plot-card">
            <div className="card-heading split-heading">
              <div>
                <span className="eyebrow">
                  Strongest supporting {plottedSignal?.method ?? selected.anchorMethod} profile
                </span>
                <h3>
                  {plottedSignal?.method === "MAXCHI"
                    ? "Raw chi-square pair profiles"
                    : plottedSignal?.method === "CHIMAERA"
                      ? "Target-specific CHIMAERA chi-square profile"
                      : plottedSignal?.method === "GENECONV"
                        ? "GENECONV Karlin–Altschul fragment envelope"
                      : plottedSignal?.method === "3SEQ"
                          ? "Target-specific hypergeometric random walks"
                      : plottedSignal?.method === "BOOTSCAN"
                            ? "Strict closest-pair bootstrap support"
                          : plottedSignal?.method === "SISCAN"
                            ? "SISCAN sister-pair permutation Z scores"
                          : "Information-rich sliding window"}
                </h3>
              </div>
              <span className="fidelity-badge">
                {plottedSignal?.method === "MAXCHI"
                  ? `${results.maxChiWindowSites} target variable sites`
                  : plottedSignal?.method === "CHIMAERA"
                    ? `${results.chimaeraWindowSites} target information-rich sites`
                    : plottedSignal?.method === "GENECONV"
                      ? `G${results.geneconvMismatchScale} · overlap ${results.geneconvMaxOverlaps}`
                      : plottedSignal?.method === "3SEQ"
                        ? `${plottedSignal.threeSeqDiscovery?.informationRichSites ?? plottedSignal.informativeSites} information-rich sites`
                        : plottedSignal?.method === "BOOTSCAN"
                          ? `${results.bootscanWindowSites}/${results.bootscanStepSites} sites · ${results.bootscanBootstrapReplicates} replicates`
                        : plottedSignal?.method === "SISCAN"
                          ? `${results.siscanWindowSites}/${results.siscanStepSites} sites · ${results.siscanScanPermutations}/${results.siscanPValuePermutations} permutations`
                        : `${results.windowSites} information-rich sites`}
              </span>
            </div>
            {plottedSignal ? (
              <SignalPlot plot={plot} signal={plottedSignal} loading={plotLoading} />
            ) : (
              <div className="plot-placeholder">The anchor signal is unavailable.</div>
            )}
          </div>

          {alignmentOpen ? (
            <EventAlignmentInspector
              event={selected}
              view={alignmentView}
              loading={alignmentLoading}
              error={alignmentError}
              flankSites={alignmentFlankSites}
              onFlankSites={setAlignmentFlankSites}
              onClose={() => setAlignmentOpen(false)}
            />
          ) : (
            <section className="alignment-launch-card">
              <div>
                <span className="eyebrow">Graphical breakpoint check</span>
                <h3>Inspect the original alignment at both boundaries</h3>
                <p>
                  Load a bounded role-, group-, and evidence-focused view only when you need it;
                  the full alignment remains inside the analysis worker. This checkpoint applies
                  the supplied CheckEnds warning path and keeps the independent BURT/BenHMM
                  99%/95% ranges visible alongside the parent-pattern review brackets.
                </p>
              </div>
              <button
                className="button button-secondary"
                type="button"
                onClick={() => setAlignmentOpen(true)}
              >
                <GitCompareArrows size={16} /> Open alignment context
              </button>
            </section>
          )}

          {phylproOpen ? (
            <EventPhylproInspector
              event={selected}
              alignmentLength={alignmentLength}
              view={phylproView}
              loading={phylproLoading}
              error={phylproError}
              windowSites={phylproWindowSites}
              gapMode={phylproGapMode}
              includeSelf={phylproIncludeSelf}
              onWindowSites={setPhylproWindowSites}
              onGapMode={setPhylproGapMode}
              onIncludeSelf={setPhylproIncludeSelf}
              onClose={() => setPhylproOpen(false)}
            />
          ) : (
            <section className="alignment-launch-card">
              <div>
                <span className="eyebrow">PHYLPRO breakpoint review</span>
                <h3>Compare left/right phylogenetic-profile correlations</h3>
                <p>
                  Lazily run the supplied PHYLPRO distance-regression calculation for the current
                  recombinant and parents. It is a diagnostic profile: RDP5 does not implement a
                  PHYLPRO significance test, so no synthetic p-value or discovery event is added.
                </p>
              </div>
              <button
                className="button button-secondary"
                type="button"
                onClick={() => setPhylproOpen(true)}
              >
                <GitCompareArrows size={16} /> Open PHYLPRO profile
              </button>
            </section>
          )}

          <div className="event-evidence-grid">
            <section className="co-group-card">
              <div className="co-group-heading">
                <div>
                  <span className="eyebrow">Current co-recombinant group</span>
                  <h3>
                    {selected.groupManualAdjusted
                      ? "Manually corrected descendants"
                      : "Automatic two-of-three group"}
                  </h3>
                </div>
                {!editingGroup ? (
                  <button
                    className="button button-quiet button-compact"
                    type="button"
                    onClick={() => {
                      setGroupDraft(selected.coRecombinantSequenceIndices);
                      setEditingGroup(true);
                    }}
                    disabled={decisionBlocked}
                  >
                    <Pencil size={14} /> Correct group
                  </button>
                ) : null}
              </div>

              {editingGroup ? (
                <div className="co-group-editor">
                  <p>
                    Select sequences descended from the same ancestral recombinant. The current
                    recombinant stays included; parent representatives are unavailable.
                  </p>
                  <label className="co-group-search">
                    <Search size={15} />
                    <input
                      type="search"
                      value={groupSearch}
                      placeholder="Find a sequence name or row"
                      onChange={(event) => setGroupSearch(event.target.value)}
                    />
                  </label>
                  <div className="co-group-options">
                    {groupMatches.map((sequence) => {
                      const checked = groupDraft.includes(sequence.index) ||
                        sequence.index === selected.recombinant;
                      const automatic = selected.automaticCoRecombinantSequenceIndices.includes(
                        sequence.index,
                      );
                      return (
                        <label className={checked ? "is-selected" : ""} key={sequence.index}>
                          <input
                            type="checkbox"
                            checked={checked}
                            disabled={sequence.index === selected.recombinant}
                            onChange={(event) => setGroupDraft((current) =>
                              event.target.checked
                                ? [...new Set([...current, sequence.index])]
                                : current.filter((index) => index !== sequence.index),
                            )}
                          />
                          <span title={sequence.name}>{sequence.name}</span>
                          <small>
                            {sequence.index === selected.recombinant
                              ? "recombinant"
                              : automatic
                                ? "2-of-3"
                                : maskedSequences.has(sequence.index)
                                  ? "masked"
                                  : `row ${sequence.index + 1}`}
                          </small>
                        </label>
                      );
                    })}
                    {!groupMatches.length ? (
                      <p className="empty-evidence">No eligible sequence matches this search.</p>
                    ) : null}
                  </div>
                  {groupError ? <p className="editor-error" role="alert">{groupError}</p> : null}
                  <div className="co-group-editor-actions">
                    <button className="button button-primary" type="button" onClick={saveGroupCorrection}>
                      <Save size={15} /> Save group
                    </button>
                    {selected.groupManualAdjusted ? (
                      <button className="button button-secondary" type="button" onClick={restoreAutomaticGroup}>
                        <RotateCcw size={15} /> Restore automatic
                      </button>
                    ) : null}
                    <button className="button button-quiet" type="button" onClick={() => setEditingGroup(false)}>
                      Cancel
                    </button>
                  </div>
                </div>
              ) : (
                <>
                  <div className="sequence-chips">
                    {selected.coRecombinantSequenceNames.map((name, index) => (
                      <span className="consensus-chip" key={`${index}-${name}`}>{name}</span>
                    ))}
                  </div>
                  {selected.groupManualAdjusted ? (
                    <p className="co-group-baseline">
                      Automatic evidence suggested {selected.automaticCoRecombinantSequenceNames.length}
                      {" "}sequence{selected.automaticCoRecombinantSequenceNames.length === 1 ? "" : "s"};
                      the saved manual group drives re-screening and alignment exports.
                    </p>
                  ) : null}
                </>
              )}
            </section>
            <section>
              <span className="eyebrow">Masked sequence follow-up</span>
              <h3>Similar RDP profiles</h3>
              {selected.traceEvidence.length ? (
                <div className="trace-list">
                  {selected.traceEvidence.map((trace) => (
                    <div key={trace.sequenceIndex}>
                      <span title={trace.sequenceName}>{trace.sequenceName}</span>
                      <strong>{trace.significant ? "significant" : "trace"}</strong>
                      <small>{pValue(trace.correctedPValue)}</small>
                    </div>
                  ))}
                </div>
              ) : <p className="empty-evidence">No structurally matching masked-sequence trace was found.</p>}
            </section>
          </div>

          <section className="hypothesis-card">
            <div className="card-heading split-heading">
              <div>
                <span className="eyebrow">Three role hypotheses</span>
                <h3>Reconciliation evidence by presumed recombinant</h3>
              </div>
              <span className="fidelity-badge">Three evidence sets · 2-of-3 rule</span>
            </div>
            <div className="hypothesis-grid">
              {selected.roleHypotheses.map((hypothesis, index) => (
                <article className={index === 0 ? "is-current" : ""} key={hypothesis.presumedRecombinant}>
                  <span>{index === 0 ? "Current role" : "Alternative role"}</span>
                  <strong title={hypothesis.presumedRecombinantName}>{hypothesis.presumedRecombinantName}</strong>
                  <small>parents: {hypothesis.parentOneName} · {hypothesis.parentTwoName}</small>
                  <div>
                    <span><b>{hypothesis.detectableSignalSetIndices.length}</b> detectable</span>
                    <span><b>{hypothesis.distanceCorrelationSetIndices.length}</b> correlated</span>
                    <span><b>{hypothesis.phylogeneticCorrelationSetIndices.length}</b> tree</span>
                    <span><b>{hypothesis.completeTwoOfThreeSetIndices.length}</b> 2 of 3</span>
                  </div>
                </article>
              ))}
            </div>
          </section>

          <section className="distance-correlation-card">
            <div className="card-heading split-heading">
              <div>
                <span className="eyebrow">Distance-correlation set</span>
                <h3>Six-value paired-matrix evidence</h3>
              </div>
              <span className="fidelity-badge">
                {currentHypothesis.validSequences} / {currentHypothesis.testedSequences} testable
              </span>
            </div>
            <p className="evidence-method-note">
              Each row compares the candidate with the presumed recombinant across the two sides
              of each breakpoint and the tract/outside pair. The strongest positive direct-polarity
              or category-relabelled correlation is shown. Native dominant-category warnings suppress
              ambiguous breakpoint pairs; MakeACOR topology affinity (or the exact dual-correlation
              override) gates the positive P-score aggregate before MakeRList group membership is assigned.
              FinalTrim then removes duplicated correlation pairs, performs its fixed-point and expansion
              passes, and hands active membership to ConsensusOK before the two-of-three group is formed.
            </p>
            <div className="distance-evidence-list">
              {currentHypothesis.distanceCorrelationEvidence.slice(0, 16).map((evidence) => {
                const pair = evidence.bestMatrixPair;
                const correlation = pair === null ? null : evidence.correlations[pair];
                const probability = pair === null ? null : evidence.pValues[pair];
                return (
                  <div key={evidence.sequenceIndex}>
                    <span title={evidence.sequenceName}>{evidence.sequenceName}</span>
                    <small>
                      {pair === null && evidence.warningFiltered.some(Boolean)
                        ? "warning-filtered"
                        : matrixPairLabel(pair)}
                    </small>
                    <strong>{correlation === null ? "—" : `r ${correlation.toFixed(3)}`}</strong>
                    <small>
                      {probability === null
                        ? "not testable"
                        : evidence.inversionCodes[pair!]
                          ? `class ${evidence.inversionCodes[pair!]} · P ${pValue(probability)}`
                          : `Σ ${evidence.aggregateScore.toFixed(2)} / ${evidence.aggregateTarget.toFixed(2)}`}
                    </small>
                    <small>
                      {evidence.duplicateFiltered.some(Boolean)
                        ? evidence.duplicateCleanedSupport
                          ? "duplicate pair suppressed"
                          : "duplicate-only direct evidence"
                        : evidence.inverseSupport
                        ? evidence.strippedInverseOnly
                          ? "inverse stripped"
                          : "inverse + direct"
                        : evidence.strongCorrelationOverride
                          ? "dual-r override"
                          : evidence.acceptableAffinity
                            ? "affinity pass"
                            : "affinity blocked"}
                    </small>
                    <i className={evidence.significant ? "is-significant" : "is-detectable-only"}>
                      {evidence.significant
                        ? "correlated"
                        : evidence.strippedInverseOnly
                          ? "inverse only"
                        : evidence.overlapEligible
                          ? "signal only"
                          : "low overlap"}
                    </i>
                  </div>
                );
              })}
            </div>
            {currentHypothesis.distanceCorrelationEvidence.length > 16 ? (
              <p className="evidence-overflow">
                {currentHypothesis.distanceCorrelationEvidence.length - 16} additional evidence rows are retained in the project export.
              </p>
            ) : null}
          </section>

          <section className="phylogenetic-card">
            <div className="card-heading split-heading">
              <div>
                <span className="eyebrow">Phylogenetic-correlation set</span>
                <h3>Six-region bootstrap neighbour joining</h3>
              </div>
              <span className="fidelity-badge">
                {selected.treePanel.sequenceCount} sequences · 10 replicates
              </span>
            </div>
            <p className="evidence-method-note">
              Each candidate must group more closely with the presumed recombinant than either
              parent in both trees of a paired region. Internal branches below 50% bootstrap support
              are collapsed. The analytical comparisons use the desktop RDP calculation path and its
              midpoint-rooted ultrametric topology ranks after Clearcut float NJ and SEQBOOT2
              resampling (seed {selected.treePanel.randomSeed}); breakpoint flanks target{" "}
              {selected.treePanel.flankVariableSiteTarget} informative variable sites on each side.{" "}
              {selected.treePanel.subsampled
                ? `The closest ${selected.treePanel.sequenceCap} sequences form the tree panel; remaining active sequences use the marked Jukes–Cantor fallback.`
                : "All active sequences are in the tree panel."}
            </p>
            <div className="tree-region-strip">
              {selected.treePanel.regions.map((region) => (
                <span className={region.usable ? "is-usable" : ""} key={region.name}>
                  <b>{region.name.replaceAll("-", " ")}</b>
                  {region.usable
                    ? `${region.supportedInternalBranches}/${region.internalBranches} supported · ${region.collapsedDistanceRankLevels} ranks`
                    : `${region.sites} sites · fallback`}
                </span>
              ))}
            </div>
            {treesOpen ? (
              <EventTreeInspector
                event={selected}
                view={treeView}
                loading={treeLoading}
                error={treeError}
                onClose={() => setTreesOpen(false)}
              />
            ) : (
              <button
                className="tree-inspector-launch"
                type="button"
                onClick={() => setTreesOpen(true)}
              >
                <GitBranch size={16} />
                <span>
                  <strong>Compare the six regional topologies</strong>
                  <small>Lazy edge-list view · inside/outside and both breakpoint pairs</small>
                </span>
                <ArrowRight size={15} />
              </button>
            )}
            <div className="distance-evidence-list phylogenetic-evidence-list">
              {currentHypothesis.phylogeneticCorrelationEvidence.slice(0, 16).map((evidence) => (
                <div key={evidence.sequenceIndex}>
                  <span title={evidence.sequenceName}>{evidence.sequenceName}</span>
                  <small>{matrixPairLabel(evidence.bestTreePair)}</small>
                  <strong>{evidence.supportingTreePairs} / 3 pairs</strong>
                  <small>
                    {evidence.disabledExcluded
                      ? "disabled from secondary evidence"
                      : evidence.distanceFallback
                        ? "JC fallback"
                        : "bootstrap NJ"}
                  </small>
                  <i className={evidence.included ? "is-significant" : "is-detectable-only"}>
                    {evidence.included ? "tree-correlated" : "other evidence"}
                  </i>
                </div>
              ))}
            </div>
            {currentHypothesis.phylogeneticCorrelationEvidence.length > 16 ? (
              <p className="evidence-overflow">
                {currentHypothesis.phylogeneticCorrelationEvidence.length - 16} additional tree-evidence rows are retained in the project export.
              </p>
            ) : null}
          </section>

          <section className="late-matrix-card">
            <div className="card-heading split-heading">
              <div>
                <span className="eyebrow">Late native consensus · active primary-RDP recheck</span>
                <h3>FinalTrim and ConsensusOK membership</h3>
              </div>
              <span className="fidelity-badge">OKSeq 0–18 · six late rechecks</span>
            </div>
            <p className="evidence-method-note">
              The supplied list-build and selected-role paths now drive grouping: duplicate cleanup and the
              nearest-nonrecombinant fixed point feed both ordered FinalTrim expansion passes,
              selected-role pruning, OKSeq 15, completed CScore, and ConsensusOK’s primary,
              distance-equivalence, and straggler rebuilds, followed by the shared selected-role
              raw/direct-tree outlier cleanup and strict inlier admission. Empty-role fallback,
              swap-last list removals, the inherited third-list loop index, exact matrix comparisons,
              and the native Long accumulator rounding are preserved. OKSeq 10 and 11 remain
              explicit source-zero slots. An asterisk on a breakpoint class marks rejection by the
              supplied opening topology-consistency check. Each finalized list candidate is then
              rerun through the primary RDP profile with the supplied LowP × 100,000 threshold lift;
              emitted, event-overlapping, and ordinarily significant results remain distinct. The
              same finalized list also receives a missing-data-aware FastRecCheckMC2 strongest-peak
              MaxChi recheck, retaining raw, within-triplet, and project-corrected probabilities as
              separate values, a fused FastRecCheckChim pass across all three target rotations, and
              the ordinary six-track GENECONV KA fragment pass with its active overlap rule. The
              supplied TSXOver(1) 3SEQ Findall pass evaluates both walk orientations and retains
              its inverse-interval list copies. A fixed-region SISCAN pass retains the nearest
              WPGMA outlier, vertical-permutation Z score, and every probability-adjustment stage.
              These late rechecks do not change reconciled event coordinates; they are separate from
              the active discovery scheduler. The MaxChi and CHIMAERA results remain
              visibly related-method evidence rather than two independent confirmations. The
              source-shaped 3SEQ recheck is also auditable here but remains golden-unvalidated.
            </p>
            <div className="late-matrix-scroll" tabIndex={0}>
              <div className="late-matrix-heading" aria-hidden="true">
                <span>Sequence</span>
                <span>Collapsed</span>
                <span>Raw tree</span>
                <span>Whole JC</span>
                <span>5′ JC</span>
                <span>3′ JC</span>
                <span>Event JC</span>
                <span>Match ×</span>
                <span>BP class</span>
                <span>Raw matrix Σ</span>
                <span>OK15</span>
                <span>CScore</span>
                <span>Final list</span>
                <span>RDP recheck</span>
                <span>MaxChi</span>
                <span>CHIMAERA</span>
                <span>GENECONV</span>
                <span>3SEQ</span>
                <span>SISCAN</span>
              </div>
              <div className="late-matrix-list">
                {currentHypothesis.distanceCorrelationEvidence
                  .filter((evidence) => evidence.finalTrimMatrix.appliesToNonrepresentative)
                  .slice(0, 16)
                  .map((evidence) => {
                    const matrix = evidence.finalTrimMatrix;
                    const match = evidence.calcMatch;
                    const consensusScore = evidence.consensusScore;
                    const rdpRecheck = evidence.postGroupRdpRecheck;
                    const maxChi = evidence.postGroupMaxChiRecheck;
                    const chimaera = evidence.postGroupChimaeraRecheck;
                    const geneconv = evidence.postGroupGeneconvRecheck;
                    const threeSeq = evidence.postGroupThreeSeqRecheck;
                    const siscan = evidence.postGroupSiscanRecheck;
                    const checkpointSummary = match.status === "complete-active-rff0"
                      ? `Beginning L/C/R ${roleScore(match.checkpointMatches[0])} / ${roleScore(match.checkpointMatches[4])} / ${roleScore(match.checkpointMatches[1])}; ending L/C/R ${roleScore(match.checkpointMatches[2])} / ${roleScore(match.checkpointMatches[5])} / ${roleScore(match.checkpointMatches[3])}; raw class ${match.rawBreakpointMatchClass}; ${match.topologyFiltered ? "rejected by the ConsensusOK topology check" : "topology-consistent"}${match.topologyDistanceFallback ? " with bounded JC fallback" : ""}; smoothing half-window ${match.smoothingHalfWindow}`
                      : "CalcMatchY unavailable: insufficient variable sites or the supplied three-alignment-length fragment bound was reached";
                    const detectedRegionSummary = matrix.detectedEventMatch
                      ? `Best direct event signal ${matrix.detectedEventSignalId === null ? "unknown" : matrix.detectedEventSignalId + 1}, tract ${matrix.detectedEventBeginning}–${matrix.detectedEventEnding}, ${(matrix.detectedEventOverlap * 100).toFixed(1)}% symmetric overlap. MatchMat order: anchor/candidate, anchor/parent 0, active bare-CompMat parent-1 lookup, anchor/parent 1, candidate/parent 0, candidate/parent 1, parent pair = ${matrix.detectedRegionMatchDistances.map(roleScore).join(" / ")}.${matrix.detectedRegionSaturated ? " Anchor/candidate comparison saturated at source sentinel 3." : ""}`
                      : `No qualifying direct event-catalogue tract for this pre-StripDupInv RList row. MatchMat fallback values = ${matrix.detectedRegionMatchDistances.map(roleScore).join(" / ")}.`;
                    const consensusStage = consensusScore.consensusPrimaryMember
                      ? "primary"
                      : consensusScore.consensusEquivalentMember
                        ? "equivalent"
                        : consensusScore.consensusStragglerMember
                          ? "straggler"
                          : consensusScore.consensusRebuiltMember
                            ? "restored"
                            : "out";
                    const finalListStage = consensusScore.finalDistanceMember
                      ? consensusScore.selectedTreeCleanupAdded
                        ? "tree add"
                        : consensusStage
                      : consensusScore.selectedTreeCleanupPrunedOut
                        ? "tree out"
                        : "out";
                    const consensusScoreSummary = `Completed native CScore: corrected P (OK0) ${consensusScore.correctedCorrelationPValue.toExponential(3)}; event overlap (OK1) ${consensusScore.detectedEventOverlap.toFixed(3)}; detectable set (OK2) ${consensusScore.detectableSetMember ? "yes" : "no"}; pattern (OK3) ${roleScore(consensusScore.patternScore)}; initial RList (OK4) ${consensusScore.initialRListMember ? "yes" : "no"}; duplicate-cleaned list (OK5) ${consensusScore.duplicateCleanedMember ? "yes" : "no"}; nearest-nonrecombinant membership (OK6) ${consensusScore.nearestNonrecombinantMember ? "yes" : "no"}; first expansion ${consensusScore.finalTrimFirstExpansionAdded ? "added" : "no"}; second expansion ${consensusScore.finalTrimSecondExpansionAdded ? "added" : "no"}; selected-role prune ${consensusScore.selectedRolePrunedOut ? "removed" : "no"}; final FinalTrim membership (OK15) ${consensusScore.finalTrimMember ? "yes" : "no"}; base before OK15 ${roleScore(consensusScore.baseScoreBeforeFinalMembership)}; after OK15 ${roleScore(consensusScore.scoreAfterFinalMembership)}; RCorrX ${roleScore(consensusScore.maximumDirectCorrelation)}; after RCorrX ${roleScore(consensusScore.scoreAfterRcorrx)}; source Long NS multiplier ${roleScore(consensusScore.sourceLongMatrixMultiplier)}; final CScore ${roleScore(consensusScore.finalScore)}; ConsensusOK stage ${consensusStage}${consensusScore.consensusFallbackRestored ? "; one empty role restored all three pre-rebuild lists" : ""}; selected-tree cleanup ${consensusScore.selectedTreeCleanupPrunedOut ? "removed" : consensusScore.selectedTreeCleanupAdded ? "added" : "unchanged"}; final distance membership ${consensusScore.finalDistanceMember ? "yes" : "no"}.`;
                    const rdpRecheckLabel = rdpRecheck.status === "representative-skipped"
                      ? "rep"
                      : rdpRecheck.status === "not-in-final-distance-list"
                        ? "—"
                        : rdpRecheck.status === "profile-unavailable"
                          ? "n/a"
                          : rdpRecheck.eventRedetected
                            ? rdpRecheck.significant ? "sig" : "trace"
                            : "none";
                    const rdpRecheckSummary = rdpRecheck.status === "complete"
                      ? `Primary RDP post-group recheck completed with local-P cutoff ${rdpRecheck.localPValueCutoff.toExponential(3)}: ${rdpRecheck.emittedSignalCount} emitted signal(s), ${rdpRecheck.candidateSignalCount} with this candidate as recombinant, and ${rdpRecheck.overlappingSignalCount} overlapping the event by more than 30%. ${rdpRecheck.eventRedetected ? `Best tract ${rdpRecheck.bestBeginning}–${rdpRecheck.bestEnding}${rdpRecheck.bestWrapsOrigin ? " (origin-wrapping)" : ""}, overlap ${(rdpRecheck.bestOverlap * 100).toFixed(1)}%, local P ${rdpRecheck.bestLocalPValue?.toExponential(3)}, corrected P ${rdpRecheck.bestCorrectedPValue?.toExponential(3)} (${rdpRecheck.significant ? "ordinary project cutoff passed" : "retained as a loosened-threshold trace"}).` : "The RDP profile was available but did not redetect an overlapping candidate-recombinant tract."}`
                      : rdpRecheck.status === "profile-unavailable"
                        ? "The native-shaped post-group RDP recheck was requested, but this triplet did not have enough informative sites for an RDP profile."
                        : rdpRecheck.status === "representative-skipped"
                          ? "The supplied loop skips the role representative itself."
                          : "This sequence was not in the finalized distance list, so the supplied post-group loop does not recheck it.";
                    const maxChiLabel = maxChi.status === "representative-skipped"
                      ? "rep"
                      : maxChi.status === "not-in-final-distance-list" ? "—"
                      : maxChi.status === "profile-unavailable"
                        ? "n/a"
                        : maxChi.sourceRecheckHit ? "hit" : "none";
                    const maxChiSummary = maxChi.status === "complete-active-unvalidated"
                      ? maxChi.bestPair === null
                        ? `The MaxChi profile retained ${maxChi.variableSites} variable sites and a ${maxChi.halfWindow}-site half-window, but no peak passed the source critical match-difference screen.`
                        : `FastRecCheckMC2 strongest-peak recheck: ${maxChi.variableSites} variable sites; half-window ${maxChi.halfWindow}, grown to ${maxChi.grownHalfWindow}; pair ${maxChi.bestPair}; peak alignment position ${maxChi.peakAlignmentPosition}; χ² ${roleScore(maxChi.maximumChiSquare ?? 0)}; raw tail ${maxChi.localPValue?.toExponential(3)}; within-triplet P ${maxChi.withinTripletPValue?.toExponential(3)}; corrected P ${maxChi.correctedPValue?.toExponential(3)}${maxChi.bonferroniApplied ? ` across ${maxChi.correctionTests.toLocaleString()} initial scan-plan opportunities` : " with project correction disabled"} (${maxChi.sourceRecheckHit ? "source recheck cutoff passed" : "cutoff not passed"})${maxChi.missingDataWindowFilterApplied ? "; MissingData/prior-erasure windows filtered" : ""}${maxChi.linearEdgeWindowFilterApplied ? "; linear-edge windows filtered" : ""}. This late-list corroboration is separate from MaxChi event discovery.`
                      : maxChi.status === "profile-unavailable"
                        ? `The MaxChi recheck was requested, but only ${maxChi.variableSites} usable variable sites remained after gaps, missing data, and earlier erasures.`
                        : maxChi.status === "representative-skipped"
                          ? "The supplied loop skips the role representative itself."
                        : "This sequence was not in the finalized distance list, so no late MaxChi recheck was requested.";
                    const chimaeraLabel = chimaera.status === "representative-skipped"
                      ? "rep"
                      : chimaera.status === "not-in-final-distance-list" ? "—"
                      : chimaera.status === "profile-unavailable"
                        ? "n/a"
                        : chimaera.sourceRecheckHit ? "hit" : "none";
                    const chimaeraSummary = chimaera.status === "complete-active-unvalidated"
                      ? chimaera.bestTarget === null
                        ? `FastRecCheckChim screened ${chimaera.targetProfilesScanned} usable target profile(s), but no peak passed the source critical match-difference screen.`
                        : `FastRecCheckChim three-target strongest-peak recheck: target ${chimaera.bestTarget}; ${chimaera.informationRichSites} information-rich sites; half-window ${chimaera.halfWindow}, grown to ${chimaera.grownHalfWindow}; peak alignment position ${chimaera.peakAlignmentPosition}; χ² ${roleScore(chimaera.maximumChiSquare ?? 0)}; raw tail ${chimaera.localPValue?.toExponential(3)}; within-triplet P ${chimaera.withinTripletPValue?.toExponential(3)}; corrected P ${chimaera.correctedPValue?.toExponential(3)}${chimaera.bonferroniApplied ? ` across ${chimaera.correctionTests.toLocaleString()} initial scan-plan opportunities` : " with project correction disabled"} (${chimaera.sourceRecheckHit ? "source recheck cutoff passed" : "cutoff not passed"})${chimaera.missingDataWindowFilterApplied ? "; MissingData/prior-erasure windows filtered" : ""}${chimaera.linearEdgeWindowFilterApplied ? "; linear-edge windows filtered" : ""}. This is related-method evidence beside MaxChi, not independent confirmation.`
                      : chimaera.status === "profile-unavailable"
                        ? "The CHIMAERA recheck was requested, but none of the target rotations retained a usable information-rich profile."
                        : chimaera.status === "representative-skipped"
                          ? "The supplied loop skips the role representative itself."
                          : "This sequence was not in the finalized distance list, so no late CHIMAERA recheck was requested.";
                    const geneconvLabel = geneconv.status === "representative-skipped"
                      ? "rep"
                      : geneconv.status === "not-in-final-distance-list" ? "—"
                      : geneconv.status === "profile-unavailable"
                        ? "n/a"
                        : geneconv.sourceRecheckHit ? "hit" : "none";
                    const geneconvSummary = geneconv.status === "complete-active-unvalidated"
                      ? geneconv.bestTrack === null
                        ? `The ordinary GENECONV profile screened ${geneconv.tracksScreened} track(s) and ${geneconv.fragmentsScored} positive start(s), but no fragment survived the critical, probability, and overlap gates.`
                        : `GCXoverD six-track ordinary-kernel recheck: track ${geneconv.bestTrack}; provisional recombinant local role ${geneconv.recombinantLocal}; tract ${geneconv.beginning}–${geneconv.ending}${geneconv.wrapsOrigin ? " (origin-wrapping)" : ""}; fragment score ${geneconv.fragmentScore} above strict critical ${geneconv.criticalScore}; raw KA P ${geneconv.rawPValue?.toExponential(3)}; corrected P ${geneconv.correctedPValue?.toExponential(3)}${geneconv.bonferroniApplied ? ` across ${geneconv.correctionTests.toLocaleString()} initial scan-plan opportunities` : " with project correction disabled"} (${geneconv.sourceRecheckHit ? "source recheck cutoff passed" : "cutoff not passed"}); ${geneconv.polymorphicSites} polymorphic sites, ${geneconv.qualifiedFragments} above-critical fragments, ${geneconv.overlapRejectedFragments} overlap rejection(s)${geneconv.numericalFallbackTracks ? `, ${geneconv.numericalFallbackTracks} bounded numerical fallback(s)` : ""}. This corroboration does not move event coordinates.`
                      : geneconv.status === "profile-unavailable"
                        ? `The GENECONV recheck was requested, but no usable six-track profile remained${geneconv.sourceSkewFilterRejected ? " after the supplied skew gate" : " after gaps and prior erasures"}.`
                        : geneconv.status === "representative-skipped"
                          ? "The supplied loop skips the role representative itself."
                          : "This sequence was not in the finalized distance list, so no late GENECONV recheck was requested.";
                    const threeSeqLabel = threeSeq.status === "representative-skipped"
                      ? "rep"
                      : threeSeq.status === "not-in-final-distance-list" ? "—"
                      : threeSeq.status === "profile-unavailable"
                        ? "n/a"
                        : threeSeq.sourceRecheckHit ? "hit" : "none";
                    const threeSeqSummary = threeSeq.status === "complete-active-unvalidated"
                      ? threeSeq.bestTarget === null
                        ? `TSXOver(1) screened ${threeSeq.targetProfilesScanned} usable target profile(s), but neither Findall orientation cleared the corrected source gate.`
                        : `TSXOver(1) Findall recheck: target ${threeSeq.bestTarget}, ${threeSeq.bestDirection} walk; tract ${threeSeq.beginning}–${threeSeq.ending}${threeSeq.wrapsOrigin ? " (origin-wrapping)" : ""}; ${threeSeq.informationRichSites} information-rich sites; parent matches ${threeSeq.parentOneMatches}/${threeSeq.parentTwoMatches}; probability/boundary excursions ${threeSeq.probabilityExcursion}/${threeSeq.maximumExcursion}; raw P ${threeSeq.rawPValue?.toExponential(3)}; corrected P ${threeSeq.correctedPValue?.toExponential(3)}${threeSeq.correctionApplied ? ` across ${threeSeq.correctionTests.toLocaleString()} initial scan-plan opportunities` : " with project correction disabled"}; ${threeSeq.qualifyingOrientations} qualifying orientation(s), ${threeSeq.sourceListEntries} source list entries${threeSeq.missingDataSplitApplied ? "; CheckSplit3Seq/SubPVal applied" : ""}; ${threeSeq.exactProbability ? "exact Single-state DP" : "Siegmund fallback"}. This corroboration does not move event coordinates.`
                      : threeSeq.status === "profile-unavailable"
                        ? "The 3SEQ Findall recheck was requested, but no target retained four information-rich sites."
                        : threeSeq.status === "representative-skipped"
                          ? "The supplied loop skips the role representative itself."
                        : "This sequence was not in the finalized distance list, so no late 3SEQ recheck was requested.";
                    const siscanLabel = siscan.status === "representative-skipped"
                      ? "rep"
                      : siscan.status === "not-requested" ? "—"
                      : siscan.profileAvailable
                        ? siscan.sourceRecheckHit ? "hit" : "none"
                        : "n/a";
                    const siscanSummary = siscan.status === "complete-active-unvalidated"
                      ? `Fixed-region SISCAN recheck: nearest outlier ${siscan.outlierSequence === null ? "unavailable" : siscan.outlierSequence + 1}; background/scored pair ${siscan.globalPair}/${siscan.scoredPair ?? "none"}; ${siscan.informativeSites} retained sites; ${siscan.permutationDraws} vertical-permutation draws; maximum Z ${roleScore(siscan.maximumZ)}; normal-tail P ${siscan.normalTailPValue.toExponential(3)}; region-adjusted P ${siscan.regionLengthAdjustedPValue.toExponential(3)}; window-adjusted P ${siscan.windowAdjustedPValue.toExponential(3)}; corrected P ${siscan.correctedPValue.toExponential(3)}${siscan.bonferroniApplied ? ` across ${siscan.correctionTests.toLocaleString()} initial scan-plan opportunities` : " with project correction disabled"} (${siscan.sourceRecheckHit ? "source recheck cutoff passed" : "cutoff not passed"}). This corroboration does not move event coordinates.`
                      : siscan.status === "representative-skipped"
                        ? "The supplied loop skips the role representative itself."
                        : siscan.status === "not-requested"
                          ? "This sequence was not in the finalized distance list, or SISCAN confirmation was disabled."
                          : siscan.status === "outlier-unavailable"
                            ? "No eligible distinct-origin fourth sequence remained for the SISCAN recheck."
                            : "The fixed event region did not retain a usable SISCAN variable-pattern profile.";
                    return (
                      <div key={evidence.sequenceIndex}>
                        <span title={evidence.sequenceName}>
                          {evidence.sequenceName}
                          <small>
                            {matrix.treeDistanceFallback ? "bounded JC fallback" : "saved tree panel"}
                          </small>
                        </span>
                        <strong>{roleScore(matrix.collapsedTreePositionScore)}</strong>
                        <strong>{roleScore(matrix.rawTreePositionScore)}</strong>
                        <strong>{roleScore(matrix.relativeDistanceScore)}</strong>
                        <strong className={!matrix.breakpointScoreAvailable[0] ? "is-unavailable" : ""}>
                          {matrix.breakpointScoreAvailable[0]
                            ? roleScore(matrix.breakpointDistanceScores[0])
                            : "warn"}
                        </strong>
                        <strong className={!matrix.breakpointScoreAvailable[1] ? "is-unavailable" : ""}>
                          {matrix.breakpointScoreAvailable[1]
                            ? roleScore(matrix.breakpointDistanceScores[1])
                            : "warn"}
                        </strong>
                        <strong title={detectedRegionSummary}>
                          {roleScore(matrix.detectedRegionDistanceScore)}
                        </strong>
                        <strong
                          className={match.status === "unavailable" ? "is-unavailable" : ""}
                          title={checkpointSummary}
                        >
                          {match.status === "complete-active-rff0"
                            ? roleScore(match.regionalMatchScore)
                            : "n/a"}
                        </strong>
                        <strong
                          className={match.status === "unavailable" ? "is-unavailable" : "calc-match-class"}
                          title={checkpointSummary}
                        >
                          {match.status === "complete-active-rff0"
                            ? `${match.breakpointMatchClass > 0 ? "+" : ""}${match.breakpointMatchClass}${match.topologyFiltered ? "*" : ""}`
                            : "n/a"}
                        </strong>
                        <b>{roleScore(matrix.activeConsensusMatrixScore)}</b>
                        <strong
                          className={consensusScore.finalTrimMember ? "membership-yes" : "membership-no"}
                          title={consensusScoreSummary}
                        >
                          {consensusScore.finalTrimMember ? "in" : "out"}
                        </strong>
                        <b title={consensusScoreSummary}>
                          {roleScore(consensusScore.finalScore)}
                        </b>
                        <strong
                          className={consensusScore.finalDistanceMember ? "membership-yes" : "membership-no"}
                          title={consensusScoreSummary}
                        >
                          {finalListStage}
                        </strong>
                        <strong
                          className={rdpRecheck.significant ? "membership-yes" : rdpRecheck.eventRedetected ? "is-detectable-only" : "membership-no"}
                          title={rdpRecheckSummary}
                        >
                          {rdpRecheckLabel}
                        </strong>
                        <strong
                          className={maxChi.sourceRecheckHit ? "membership-yes" : maxChi.profileAvailable ? "membership-no" : "is-unavailable"}
                          title={maxChiSummary}
                        >
                          {maxChiLabel}
                        </strong>
                        <strong
                          className={chimaera.sourceRecheckHit ? "membership-yes" : chimaera.profileAvailable ? "membership-no" : "is-unavailable"}
                          title={chimaeraSummary}
                        >
                          {chimaeraLabel}
                        </strong>
                        <strong
                          className={geneconv.sourceRecheckHit ? "membership-yes" : geneconv.profileAvailable ? "membership-no" : "is-unavailable"}
                          title={geneconvSummary}
                        >
                          {geneconvLabel}
                        </strong>
                        <strong
                          className={threeSeq.sourceRecheckHit ? "membership-yes" : threeSeq.profileAvailable ? "membership-no" : "is-unavailable"}
                          title={threeSeqSummary}
                        >
                          {threeSeqLabel}
                        </strong>
                        <strong
                          className={siscan.sourceRecheckHit ? "membership-yes" : siscan.profileAvailable ? "membership-no" : "is-unavailable"}
                          title={siscanSummary}
                        >
                          {siscanLabel}
                        </strong>
                      </div>
                    );
                  })}
              </div>
            </div>
            {currentHypothesis.distanceCorrelationEvidence.filter(
              (evidence) => evidence.finalTrimMatrix.appliesToNonrepresentative,
            ).length > 16 ? (
              <p className="evidence-overflow">
                Additional mapped matrix-score rows for this and the two alternative roles remain
                in project JSON.
              </p>
            ) : null}
          </section>

          <div className="review-actions">
            <div>
              <button
                type="button"
                className={`button review-accept${selected.reviewState === "accepted" ? " is-selected" : ""}`}
                onClick={() => onReviewState(selected.id, "accepted")}
                disabled={decisionBlocked}
              >
                <Check size={17} /> Accept event
              </button>
              <button
                type="button"
                className={`button review-reject${selected.reviewState === "rejected" ? " is-selected" : ""}`}
                onClick={() => onReviewState(selected.id, "rejected")}
                disabled={decisionBlocked}
              >
                <X size={17} /> Reject event
              </button>
              {selected.reviewState !== "unreviewed" ? (
                <button
                  className="button button-quiet"
                  type="button"
                  onClick={() => onReviewState(selected.id, "unreviewed")}
                  disabled={decisionBlocked}
                >
                  <RotateCcw size={16} /> Clear decision
                </button>
              ) : null}
            </div>
            <button
              className={canReconcile ? "button button-secondary" : "button button-disabled"}
              type="button"
              disabled={!canReconcile || reconciling}
              onClick={() => onReconcileAfter(selected.id)}
              title={
                pendingEvent === null
                  ? "No changed event is waiting for downstream reconciliation"
                  : selected.id !== pendingEvent
                    ? `Return to event ${pendingEvent + 1}`
                    : selected.reviewState === "unreviewed"
                      ? "Accept or reject the changed event before re-identifying later events"
                      : undefined
              }
            >
              <RefreshCw className={reconciling ? "spin" : ""} size={16} />
              {reconciling ? "Re-identifying…" : "Re-identify later events"}
            </button>
          </div>
        </div>
      </div>

      <footer className="step-actions">
        <button className="button button-quiet" type="button" onClick={onBack}>
          <ArrowLeft size={16} /> Scan summary
        </button>
        <button className="button button-primary" type="button" onClick={onExport}>
          Export analysis <ArrowRight size={16} />
        </button>
      </footer>
    </section>
  );
}
