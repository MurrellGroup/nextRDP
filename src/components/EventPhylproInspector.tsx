import { Info, X } from "lucide-react";

import type {
  EventPhylproView,
  PhylproGapMode,
  ReconciledEvent,
} from "../lib/types";

interface EventPhylproInspectorProps {
  event: ReconciledEvent;
  alignmentLength: number;
  view: EventPhylproView | null;
  loading: boolean;
  error: string;
  windowSites: number;
  gapMode: PhylproGapMode;
  includeSelf: boolean;
  onWindowSites: (value: number) => void;
  onGapMode: (value: PhylproGapMode) => void;
  onIncludeSelf: (value: boolean) => void;
  onClose: () => void;
}

const width = 920;
const height = 310;
const margin = { top: 24, right: 22, bottom: 44, left: 54 };

function correlation(value: number): string {
  return Number.isFinite(value) ? value.toFixed(4) : "—";
}

export function EventPhylproInspector({
  event,
  alignmentLength,
  view,
  loading,
  error,
  windowSites,
  gapMode,
  includeSelf,
  onWindowSites,
  onGapMode,
  onIncludeSelf,
  onClose,
}: EventPhylproInspectorProps) {
  const points = view?.points ?? [];
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;
  const rawMinimum = view ? Math.min(1, view.minimumValue) : 0;
  const rawMaximum = view ? Math.max(rawMinimum + 0.01, view.maximumValue) : 1;
  const padding = Math.max(0.02, (rawMaximum - rawMinimum) * 0.08);
  // PPRegression narrows its ratio without clamping. Preserve any tiny
  // source-rounding excursion rather than clipping it at the SVG boundary.
  const yMinimum = rawMinimum - padding;
  const yMaximum = rawMaximum + padding;
  const ySpan = Math.max(0.01, yMaximum - yMinimum);
  const x = (value: number) =>
    margin.left + ((value - 1) / Math.max(1, alignmentLength - 1)) * innerWidth;
  const y = (value: number) =>
    margin.top + ((yMaximum - value) / ySpan) * innerHeight;
  const path = (key: "recombinant" | "majorParent" | "minorParent") =>
    points
      .map((point, index) =>
        `${index ? "L" : "M"}${x(point.alignmentPosition).toFixed(2)},${y(point[key]).toFixed(2)}`
      )
      .join(" ");
  const highlight = (start: number, end: number, key: string) => (
    <rect
      key={key}
      x={x(start)}
      y={margin.top}
      width={Math.max(2, x(end) - x(start))}
      height={innerHeight}
      className="phylpro-event-region"
    />
  );

  return (
    <section className="phylpro-inspector" aria-labelledby="phylpro-title">
      <header>
        <div>
          <span className="eyebrow">Supplied PHYLPRO breakpoint review</span>
          <h3 id="phylpro-title">Left/right phylogenetic-profile correlation</h3>
          <p>
            Each line correlates one role sequence&apos;s Hamming-distance vector on the two sides
            of a moving central partition. A local drop indicates changed phylogenetic affinity.
          </p>
        </div>
        <button className="button button-quiet button-compact" type="button" onClick={onClose}>
          <X size={14} /> Close
        </button>
      </header>

      <div className="phylpro-controls">
        <label>
          <span>Total window</span>
          <input
            type="number"
            min={10}
            max={5000}
            step={2}
            value={windowSites}
            onChange={(input) => onWindowSites(Number(input.target.value))}
          />
        </label>
        <label>
          <span>Missing-data policy</span>
          <select
            value={gapMode}
            onChange={(input) => onGapMode(input.target.value as PhylproGapMode)}
          >
            <option value="ignore-missing-pairwise">Ignore pairwise</option>
            <option value="strip-any-missing-column">Strip whole column</option>
          </select>
        </label>
        <label className="phylpro-checkbox">
          <input
            type="checkbox"
            checked={includeSelf}
            onChange={(input) => onIncludeSelf(input.target.checked)}
          />
          <span>Include zero self-distance</span>
        </label>
      </div>

      {loading ? (
        <div className="plot-placeholder">Calculating the three source PHYLPRO profiles…</div>
      ) : error ? (
        <div className="plot-placeholder is-error">{error}</div>
      ) : view && points.length ? (
        <>
          <div className="phylpro-plot-wrap">
            <svg className="phylpro-plot" viewBox={`0 0 ${width} ${height}`} role="img">
              <title>PHYLPRO correlation profiles for event {event.id + 1}</title>
              <desc>
                Pearson correlations between pairwise Hamming-distance vectors in the left and
                right half-windows for the recombinant, major parent, and minor parent.
              </desc>
              {[0, 0.25, 0.5, 0.75, 1].map((tick) => {
                const value = yMinimum + tick * ySpan;
                return (
                  <g key={tick}>
                    <line
                      x1={margin.left}
                      x2={width - margin.right}
                      y1={y(value)}
                      y2={y(value)}
                      className="plot-grid"
                    />
                    <text x={margin.left - 9} y={y(value) + 4} textAnchor="end" className="plot-label">
                      {value.toFixed(2)}
                    </text>
                  </g>
                );
              })}
              {event.wrapsOrigin ? (
                <>
                  {highlight(1, event.ending, "left")}
                  {highlight(event.beginning, alignmentLength, "right")}
                </>
              ) : highlight(event.beginning, event.ending, "single")}
              <line
                x1={x(event.beginning)}
                x2={x(event.beginning)}
                y1={margin.top}
                y2={height - margin.bottom}
                className="phylpro-breakpoint"
              />
              <line
                x1={x(event.ending)}
                x2={x(event.ending)}
                y1={margin.top}
                y2={height - margin.bottom}
                className="phylpro-breakpoint"
              />
              <path d={path("recombinant")} className="phylpro-line is-recombinant" />
              <path d={path("majorParent")} className="phylpro-line is-major" />
              <path d={path("minorParent")} className="phylpro-line is-minor" />
              <text x={margin.left} y={height - 13} className="plot-label">1</text>
              <text x={width - margin.right} y={height - 13} textAnchor="end" className="plot-label">
                {alignmentLength.toLocaleString()} · alignment position
              </text>
            </svg>
          </div>
          <div className="phylpro-legend">
            <span className="is-recombinant">Recombinant · {view.sequenceNames[0]}</span>
            <span className="is-major">Major parent · {view.sequenceNames[1]}</span>
            <span className="is-minor">Minor parent · {view.sequenceNames[2]}</span>
          </div>
          <div className="phylpro-summary">
            {view.minimumBySequence.map((minimum, index) => (
              <div key={minimum.sequenceIndex}>
                <span>{index === 0 ? "Recombinant" : index === 1 ? "Major parent" : "Minor parent"} minimum</span>
                <strong>{correlation(minimum.correlation)}</strong>
                <small>at {minimum.position.toLocaleString()}</small>
              </div>
            ))}
            <div>
              <span>Context</span>
              <strong>{view.contextSequences.toLocaleString()} sequences</strong>
              <small>{view.eligibleColumns.toLocaleString()} mapped columns</small>
            </div>
            <div>
              <span>Half-window</span>
              <strong>{view.halfWindowSites.toLocaleString()} sites</strong>
              <small>{view.windowCapped ? "Capped to half the profile" : `${view.windowSites} total requested`}</small>
            </div>
            <div>
              <span>Kernel work</span>
              <strong>{view.evaluatedPoints.toLocaleString()} partitions</strong>
              <small>{view.returnedPoints.toLocaleString()} display points retained</small>
            </div>
          </div>
          <div className="inline-note phylpro-note">
            <Info size={16} />
            <p>
              RDP5 does not implement a PHYLPRO permutation/significance test, so this profile is
              diagnostic and never contributes an invented p-value or event. The supplied active
              calculation maps polymorphic columns after the gap policy despite the manual&apos;s
              broader conceptual description. Masked originals remain context; disabled rows and
              cyclic fragments do not. Browser-normalized ambiguity is treated as missing. Only
              the three displayed distance rows are updated, giving
              the same coefficients as the native all-pairs matrix with linear rather than
              quadratic sequence-count work.
            </p>
          </div>
        </>
      ) : (
        <div className="plot-placeholder">No PHYLPRO profile is available for this event.</div>
      )}
    </section>
  );
}
