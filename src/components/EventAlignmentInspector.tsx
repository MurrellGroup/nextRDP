import type { CSSProperties } from "react";

import type {
  EventAlignmentPanel,
  EventAlignmentRole,
  EventAlignmentView,
  ReconciledEvent,
} from "../lib/types";

interface EventAlignmentInspectorProps {
  event: ReconciledEvent;
  view: EventAlignmentView | null;
  loading: boolean;
  error: string;
  flankSites: number;
  onFlankSites: (value: number) => void;
  onClose: () => void;
}

const roleLabels: Record<EventAlignmentRole, string> = {
  recombinant: "recombinant",
  "major-parent": "major parent",
  "minor-parent": "minor parent",
  "co-recombinant": "co-recombinant",
  evidence: "supporting evidence",
};

function normaliseBase(base: string | undefined): string {
  return (base ?? "").toLocaleUpperCase();
}

function baseClass(
  base: string,
  recombinantBase: string,
  majorBase: string,
  minorBase: string,
): string {
  const value = normaliseBase(base);
  const recombinant = normaliseBase(recombinantBase);
  const major = normaliseBase(majorBase);
  const minor = normaliseBase(minorBase);
  if (!/^[ACGTU]$/.test(value)) return "base-missing";
  if (value === major && value === minor) return "base-shared";
  if (value === major && value !== minor) return "base-major";
  if (value === minor && value !== major) return "base-minor";
  if (value === recombinant) return "base-recombinant";
  return "base-other";
}

function nativeBreakpointWarning(panel: EventAlignmentPanel): string {
  const reasons: string[] = [];
  if (panel.uncertainDueToErasure) {
    const eventLabel =
      "prior erased event" +
      (panel.uncertainErasureEventIds.length === 1 ? " " : "s ") +
      panel.uncertainErasureEventIds.map((id) => id + 1).join(", ");
    if (panel.erasureAdjacent) {
      reasons.push("touches " + eventLabel);
    } else if (panel.nearestErasureInformativeSites === 0) {
      reasons.push("has no intervening RDP information-rich positions before " + eventLabel);
    } else {
      reasons.push(
        "native range reaches " +
          eventLabel +
          " within " +
          (panel.nearestErasureInformativeSites ?? panel.rdpWindowInformativeSites) +
          " RDP information-rich positions",
      );
    }
  }
  if (panel.inputMissingDataInCheckRange) {
    reasons.push("native range contains an input MissingData run");
  }
  if (panel.linearEdgeWithinRdpWindow) {
    reasons.push("full RDP window is unavailable at the linear edge");
  }
  if (!panel.informationProfileAvailable) {
    reasons.push("information-rich profile unavailable");
  }
  return reasons.join(" · ");
}

function confidenceRange(beginning: number, ending: number, wraps: boolean): string {
  return `${Math.abs(beginning).toLocaleString()} → ${Math.abs(ending).toLocaleString()}${wraps ? " · wraps origin" : ""}`;
}

function CoordinateRuler({ panel }: { panel: EventAlignmentPanel }) {
  const style = {
    gridTemplateColumns: `210px repeat(${panel.coordinates.length}, 16px)`,
  } satisfies CSSProperties;
  return (
    <div className="alignment-matrix-row alignment-ruler" style={style}>
      <span className="alignment-sticky-cell">Alignment coordinate</span>
      {panel.coordinates.map((coordinate, index) => {
        const labelled = index === panel.centerIndex || coordinate === 1 || coordinate % 10 === 0;
        const classNames = [
          index === panel.centerIndex ? "is-breakpoint" : "",
          coordinate === panel.parentTransition.leftInformativeCoordinate
            ? "is-transition-left"
            : "",
          coordinate === panel.parentTransition.rightInformativeCoordinate
            ? "is-transition-right"
            : "",
        ].filter(Boolean).join(" ");
        return (
          <span
            className={classNames}
            key={`${coordinate}-${index}`}
            title={`Alignment site ${coordinate.toLocaleString()}`}
          >
            {labelled ? coordinate.toLocaleString() : ""}
          </span>
        );
      })}
    </div>
  );
}

function AlignmentPanelView({
  view,
  panel,
  panelIndex,
}: {
  view: EventAlignmentView;
  panel: EventAlignmentPanel;
  panelIndex: 0 | 1;
}) {
  const recombinant = view.rows.find((row) => row.role === "recombinant")?.panels[panelIndex] ?? "";
  const major = view.rows.find((row) => row.role === "major-parent")?.panels[panelIndex] ?? "";
  const minor = view.rows.find((row) => row.role === "minor-parent")?.panels[panelIndex] ?? "";
  const style = {
    gridTemplateColumns: `210px repeat(${panel.coordinates.length}, 16px)`,
  } satisfies CSSProperties;

  return (
    <section className="alignment-panel" aria-labelledby={`alignment-${panel.name}-title`}>
      <div className="alignment-panel-heading">
        <div>
          <span>{panel.name === "beginning" ? "Left boundary" : "Right boundary"}</span>
          <h4 id={`alignment-${panel.name}-title`}>
            Breakpoint {panel.center.toLocaleString()}
          </h4>
        </div>
        <div className="alignment-panel-context">
          {panel.statisticalConfidence.intervalAvailable ? (
            <small className="is-statistical-confidence">
              BURT 99% {confidenceRange(
                panel.statisticalConfidence.confidence99.beginning,
                panel.statisticalConfidence.confidence99.ending,
                panel.statisticalConfidence.confidence99.wrapsOrigin,
              )} · 95% {confidenceRange(
                panel.statisticalConfidence.confidence95.beginning,
                panel.statisticalConfidence.confidence95.ending,
                panel.statisticalConfidence.confidence95.wrapsOrigin,
              )}
            </small>
          ) : null}
          {panel.nativeCheckEndsWarning ? (
            <small className="is-warning">{nativeBreakpointWarning(panel)}</small>
          ) : panel.parentTransition.supported ? (
            <small className="is-supported">
              expected parent switch bracketed at{" "}
              {panel.parentTransition.leftInformativeCoordinate?.toLocaleString()} →{" "}
              {panel.parentTransition.rightInformativeCoordinate?.toLocaleString()}
            </small>
          ) : (
            <small>expected parent switch not bracketed in this window</small>
          )}
          {panel.nativeCheckEndsApplied && panel.nativeCheckRange.coordinateCount > 0 ? (
            <small>
              CheckEnds range {panel.nativeCheckRange.beginning.toLocaleString()} →{" "}
              {panel.nativeCheckRange.ending.toLocaleString()}
              {panel.nativeCheckRange.wrapsOrigin ? " · wraps origin" : ""} ·{" "}
              {panel.rdpWindowInformativeSites} information-rich-site window
            </small>
          ) : null}
          {view.circular && panel.coordinates[0] > panel.coordinates.at(-1)! ? (
            <small>window crosses coordinate 1</small>
          ) : null}
        </div>
      </div>
      <div className="alignment-matrix-scroll" tabIndex={0}>
        <div className="alignment-matrix">
          <CoordinateRuler panel={panel} />
          {view.rows.map((row) => {
            const sequence = row.panels[panelIndex];
            return (
              <div className="alignment-matrix-row" style={style} key={row.sequenceIndex}>
                <span className={`alignment-sticky-cell role-${row.role}`}>
                  <span title={row.sequenceName}>{row.sequenceName}</span>
                  <small>{roleLabels[row.role]}</small>
                  <i>
                    {row.currentGroupMember ? "current group" : ""}
                    {row.trace ? `${row.currentGroupMember ? " · " : ""}masked trace` : ""}
                    {row.masked && !row.trace ? `${row.currentGroupMember ? " · " : ""}masked` : ""}
                    {row.disabled ? `${row.currentGroupMember || row.masked ? " · " : ""}disabled` : ""}
                  </i>
                </span>
                {[...sequence].map((base, index) => (
                  <code
                    className={`${baseClass(base, recombinant[index], major[index], minor[index])}${
                      index === panel.centerIndex ? " is-breakpoint" : ""
                    }`}
                    key={`${panel.coordinates[index]}-${index}`}
                    title={`${row.sequenceName} · site ${panel.coordinates[index]?.toLocaleString()} · ${normaliseBase(base) || "blank"}`}
                  >
                    {normaliseBase(base) || "·"}
                  </code>
                ))}
              </div>
            );
          })}
        </div>
      </div>
    </section>
  );
}

export function EventAlignmentInspector({
  event,
  view,
  loading,
  error,
  flankSites,
  onFlankSites,
  onClose,
}: EventAlignmentInspectorProps) {
  return (
    <section className="breakpoint-alignment-card">
      <div className="card-heading split-heading">
        <div>
          <span className="eyebrow">Graphical breakpoint check</span>
          <h3>Original-alignment context</h3>
          <p>
            Compare the recombinant and parent patterns at both inferred boundaries before
            accepting or correcting event {event.id + 1}.
          </p>
        </div>
        <div className="alignment-card-actions">
          <div className="alignment-window-control" aria-label="Alignment context on each side">
            {[15, 30, 60].map((value) => (
              <button
                className={flankSites === value ? "is-selected" : ""}
                type="button"
                onClick={() => onFlankSites(value)}
                key={value}
              >
                ±{value}
              </button>
            ))}
          </div>
          <button className="button button-quiet button-compact" type="button" onClick={onClose}>
            Close
          </button>
        </div>
      </div>

      {loading ? (
        <div className="alignment-loading" role="status">
          <span />
          Preparing bounded breakpoint windows…
        </div>
      ) : error ? (
        <div className="alignment-error" role="alert">{error}</div>
      ) : view ? (
        <>
          {view.panels.some((panel) => panel.nativeCheckEndsWarning) ? (
            <div className="alignment-context-warning">
              At least one boundary triggers the supplied RDP CheckEnds warning because its native
              information-rich-site range reaches erased sequence or input MissingData, cannot form
              a full window at a linear edge, or lacks a usable profile. RDP5 treats that position as
              uncertain; inspect the parent-pattern transition before accepting it.
            </div>
          ) : null}
          <div className="alignment-legend" aria-label="Alignment base colour legend">
            <span className="base-major">Major-parent match</span>
            <span className="base-minor">Minor-parent match</span>
            <span className="base-shared">Parents agree</span>
            <span className="base-recombinant">Recombinant-only state</span>
            <span className="base-missing">Gap or ambiguity</span>
            <span className="transition-left-key">Expected state before boundary</span>
            <span className="transition-right-key">Expected state after boundary</span>
          </div>
          <div className="alignment-panels">
            <AlignmentPanelView view={view} panel={view.panels[0]} panelIndex={0} />
            <AlignmentPanelView view={view} panel={view.panels[1]} panelIndex={1} />
          </div>
          <p className="alignment-footnote">
            Showing {view.rows.length} of {view.candidateRowCount} role, group, trace, and evidence
            rows from the original alignment{view.omittedRowCount
              ? ` · ${view.omittedRowCount} lower-priority row${view.omittedRowCount === 1 ? "" : "s"} omitted`
            : ""}. Coordinates are 1-based{view.circular ? " and wrap at the origin" : ""}.
            Parent-switch brackets are review intervals from the closest expected informative
            states, not statistical confidence intervals. The separate BURT/BenHMM card retains
            the source-labelled 99% and 95% statistical ranges, HMM positions, and breakpoint
            repositioning state. Erasure proximity is counted in RDP information-rich triplet
            positions rather than raw alignment columns{view.fragmentAssisted
              ? " · the selected call used a retained working fragment"
              : ""}.
          </p>
        </>
      ) : null}
    </section>
  );
}
