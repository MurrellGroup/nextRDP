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
  EventEdit,
  ReconciledEvent,
  ReviewState,
  ScanResults,
  SequenceSummary,
  SignalPlot as SignalPlotData,
} from "../lib/types";
import { SignalPlot } from "./SignalPlot";

interface ReviewStepProps {
  results: ScanResults;
  alignmentLength: number;
  sequences: SequenceSummary[];
  onGetPlot: (signalId: number) => Promise<SignalPlotData>;
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
  onReviewState,
  onUpdateEvent,
  onUpdateEventGroup,
  onReconcileAfter,
  reconciling,
  onSaveProject,
  onBack,
  onExport,
}: ReviewStepProps) {
  const [selectedId, setSelectedId] = useState(results.events[0]?.id ?? -1);
  const [plot, setPlot] = useState<SignalPlotData | null>(null);
  const [plotLoading, setPlotLoading] = useState(false);
  const [editing, setEditing] = useState(false);
  const [editError, setEditError] = useState("");
  const [editingGroup, setEditingGroup] = useState(false);
  const [groupDraft, setGroupDraft] = useState<number[]>([]);
  const [groupSearch, setGroupSearch] = useState("");
  const [groupError, setGroupError] = useState("");
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
      .filter((sequence): sequence is SequenceSummary => sequence !== undefined);
  }, [selected, sequences]);

  const groupMatches = useMemo(() => {
    if (!selected) return [];
    const query = groupSearch.trim().toLocaleLowerCase();
    const selectedSet = new Set(groupDraft);
    const automaticSet = new Set(selected.automaticCoRecombinantSequenceIndices);
    return sequences
      .filter((sequence) =>
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
  }, [groupDraft, groupSearch, selected, sequences]);

  if (!selected) {
    return (
      <section className="step-page empty-results">
        <div className="empty-icon"><GitCompareArrows size={28} /></div>
        <h1>No significant RDP events</h1>
        <p>No primary signal passed the current window, threshold, and correction settings.</p>
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
          <button className="button button-secondary" type="button" onClick={onSaveProject}>
            <Save size={15} /> Save project checkpoint
          </button>
        </div>
      </header>

      <div className="notice notice-blue">
        <AlertTriangle size={18} />
        <p>
          Events were found in strongest-first cyclic passes, with each inferred co-group tract erased
          and re-entered as a gap-padded fragment before the next full screen. Three evidence sets are
          evaluated for every role; native PhPr, leave-one-out, displacement, collapsed-tree, and
          TrpScore decision contributions are auditable below.
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
        <aside className="event-list" aria-label="Reconciled RDP events">
          <div className="event-list-heading">
            <span className="eyebrow">Analysis order</span>
            <strong>{results.events.length} events</strong>
          </div>
          <div className="event-list-scroll">
            {results.events.map((event) => {
              const stale = pendingEvent !== null && event.id > pendingEvent;
              return (
                <button
                  type="button"
                  key={event.id}
                  className={`event-list-item${selected.id === event.id ? " is-selected" : ""}${stale ? " is-stale" : ""}`}
                  onClick={() => setSelectedId(event.id)}
                >
                  <span className={`event-status status-${event.reviewState}`} />
                  <span>
                    <strong>Event {event.id + 1}</strong>
                    <small title={event.recombinantName}>{event.recombinantName}</small>
                  </span>
                  <span>
                    <strong>{stale ? "stale" : pValue(event.bestCorrectedPValue)}</strong>
                    <small>
                      {stale
                        ? "awaiting rebuild"
                        : `${event.supportSignalIds.length} signal${event.supportSignalIds.length === 1 ? "" : "s"}`}
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
              <span className="eyebrow">Event {selected.id + 1} of {results.events.length}</span>
              <h2>{selected.recombinantName}</h2>
            </div>
            <div className="event-navigation">
              <button type="button" onClick={() => move(-1)} aria-label="Previous event"><ChevronLeft /></button>
              <button type="button" onClick={() => move(1)} aria-label="Next event"><ChevronRight /></button>
            </div>
          </div>

          <EventSchematic event={selected} alignmentLength={alignmentLength} />

          <div className="role-grid">
            <div className="role-recombinant">
              <span>Current recombinant</span>
              <strong>{selected.recombinantName}</strong>
            </div>
            <div className="role-major">
              <span>Major-parent-like</span>
              <strong>{selected.majorParentName}</strong>
            </div>
            <div className="role-minor">
              <span>Minor-parent-like</span>
              <strong>{selected.minorParentName}</strong>
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
            <span><Layers3 size={14} /><strong>{selected.supportSignalIds.length}</strong> primary signals</span>
            <span><strong>{currentHypothesis.detectableSignalSetIndices.length}</strong> detectable</span>
            <span><strong>{currentHypothesis.distanceCorrelationSetIndices.length}</strong> distance-correlated</span>
            <span><strong>{currentHypothesis.phylogeneticCorrelationSetIndices.length}</strong> tree-correlated</span>
            <span><strong>{selected.coRecombinantSequenceIndices.length}</strong> co-recombinant</span>
            <span><strong>{selected.traceEvidence.length}</strong> masked traces</span>
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
                        <option value={sequence.index} key={sequence.index}>{sequence.name}</option>
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
                <span className="eyebrow">Strongest supporting RDP profile</span>
                <h3>Information-rich sliding window</h3>
              </div>
              <span className="fidelity-badge">{results.windowSites} variable sites</span>
            </div>
            {plottedSignal ? (
              <SignalPlot plot={plot} signal={plottedSignal} loading={plotLoading} />
            ) : (
              <div className="plot-placeholder">The anchor signal is unavailable.</div>
            )}
          </div>

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
                                : sequence.masked
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
              The active first FinalTrim pass also marks direct-correlation pairs duplicated across
              competing role lists; those diagnostics do not yet prune the two-of-three group.
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
              are collapsed. {selected.treePanel.subsampled
                ? `The closest ${selected.treePanel.sequenceCap} sequences form the tree panel; remaining active sequences use the marked Jukes–Cantor fallback.`
                : "All active sequences are in the tree panel."}
            </p>
            <div className="tree-region-strip">
              {selected.treePanel.regions.map((region) => (
                <span className={region.usable ? "is-usable" : ""} key={region.name}>
                  <b>{region.name.replaceAll("-", " ")}</b>
                  {region.usable
                    ? `${region.supportedInternalBranches}/${region.internalBranches} supported`
                    : `${region.sites} sites · fallback`}
                </span>
              ))}
            </div>
            <div className="distance-evidence-list phylogenetic-evidence-list">
              {currentHypothesis.phylogeneticCorrelationEvidence.slice(0, 16).map((evidence) => (
                <div key={evidence.sequenceIndex}>
                  <span title={evidence.sequenceName}>{evidence.sequenceName}</span>
                  <small>{matrixPairLabel(evidence.bestTreePair)}</small>
                  <strong>{evidence.supportingTreePairs} / 3 pairs</strong>
                  <small>
                    {evidence.maskedExcluded
                      ? "masked from trees"
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
