import type { RdpSignal, SignalPlot as SignalPlotData } from "../lib/types";

interface SignalPlotProps {
  plot: SignalPlotData | null;
  signal: RdpSignal;
  loading: boolean;
}

const width = 920;
const height = 280;
const margin = { top: 24, right: 22, bottom: 40, left: 48 };

export function SignalPlot({ plot, signal, loading }: SignalPlotProps) {
  if (loading) return <div className="plot-placeholder">Recomputing the selected triplet profile…</div>;
  if (!plot || plot.points.length < 2) return <div className="plot-placeholder">No plot data are available.</div>;

  const xMin = plot.points[0].alignmentPosition;
  const xMax = plot.points[plot.points.length - 1].alignmentPosition;
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;
  const randomWalk = plot.metric === "random-walk-height";
  const unitInterval = plot.metric === "pair-identity" || plot.metric === "bootstrap-support";
  const rawMinimum = randomWalk ? Math.min(0, plot.minimumValue) : 0;
  const rawMaximum = !unitInterval ? Math.max(0, plot.maximumValue) : 1;
  const rawSpan = Math.max(1, rawMaximum - rawMinimum);
  const yMinimum = randomWalk ? rawMinimum - rawSpan * 0.05 : 0;
  const yMaximum = randomWalk
    ? rawMaximum + rawSpan * 0.05
    : !unitInterval ? Math.max(1, rawMaximum * 1.05) : 1;
  const ySpan = Math.max(1, yMaximum - yMinimum);
  const x = (value: number) =>
    margin.left + ((value - xMin) / Math.max(1, xMax - xMin)) * innerWidth;
  const y = (value: number) =>
    margin.top + ((yMaximum - value) / ySpan) * innerHeight;
  const path = (key: "pair12" | "pair13" | "pair23") =>
    plot.points
      .map((point, index) => `${index ? "L" : "M"}${x(point.alignmentPosition).toFixed(2)},${y(point[key]).toFixed(2)}`)
      .join(" ");
  const chimaeraTarget = signal.chimaeraDiscovery?.targetLocal ?? null;
  const chimaeraParentOne = chimaeraTarget === null
    ? null
    : ([1, 2, 0] as const)[chimaeraTarget];
  const chimaeraTrace = chimaeraTarget === 0
    ? "pair12"
    : chimaeraTarget === 1 ? "pair23" : chimaeraTarget === 2 ? "pair13" : null;
  const visibleTraces: readonly ("pair12" | "pair13" | "pair23")[] = chimaeraTrace
    ? [chimaeraTrace]
    : ["pair12", "pair13", "pair23"];

  const highlight = (start: number, end: number, key: string) => (
    <rect
      key={key}
      x={x(start)}
      y={margin.top}
      width={Math.max(2, x(end) - x(start))}
      height={innerHeight}
      rx={4}
      className="plot-event-region"
    />
  );

  return (
    <div className="signal-plot-wrap">
      <svg className="signal-plot" viewBox={`0 0 ${width} ${height}`} role="img">
        <title>
          {signal.method === "CHIMAERA"
            ? "CHIMAERA target χ² profile"
            : signal.method === "GENECONV"
              ? "GENECONV negative log10 KA P fragment envelope"
              : signal.method === "3SEQ"
                ? "3SEQ target-specific hypergeometric random walks"
                : signal.method === "BOOTSCAN"
                  ? "BootScan strict closest-pair bootstrap support"
                : plot.metric === "chi-square" ? "MaxChi χ² profile" : "Sliding-window pairwise identity"} for signal {signal.id + 1}
        </title>
        <desc>
          {signal.method === "CHIMAERA"
            ? "One target-to-parent-one chi-square trace across the candidate recombinant's information-rich binary profile. The highlighted region is the matched tract."
            : signal.method === "GENECONV"
              ? "Three colour-matched inner and outer GENECONV fragment envelopes, measured as negative log10 raw Karlin-Altschul probability. The highlighted region is the selected fragment."
              : signal.method === "3SEQ"
                ? "Three target-specific plus-one/minus-one walks across information-rich sites. Each trace treats one triplet member as the candidate recombinant; the highlighted region is the selected maximum excursion."
                : signal.method === "BOOTSCAN"
                  ? "Three seeded sliding-window bootstrap-support curves. Each replicate votes only for its unique closest pair; the highlighted region is the significant supported-pair tract."
                : plot.metric === "chi-square"
                  ? "Three pairwise maximum chi-square traces across variable sites. The highlighted region is the matched recombinant tract."
                  : "Pairwise identity across information-rich sites for the three sequences used to detect this signal. The highlighted region is bounded by the inferred breakpoints."}
        </desc>
        {[0, 0.25, 0.5, 0.75, 1].map((tick) => (
          <g key={tick}>
            <line
              x1={margin.left}
              x2={width - margin.right}
              y1={y(yMinimum + tick * ySpan)}
              y2={y(yMinimum + tick * ySpan)}
              className="plot-grid"
            />
            <text x={margin.left - 10} y={y(yMinimum + tick * ySpan) + 4} textAnchor="end" className="plot-label">
              {!unitInterval
                ? (yMinimum + tick * ySpan).toFixed(ySpan >= 100 ? 0 : 1)
                : tick.toFixed(2)}
            </text>
          </g>
        ))}
        {signal.wrapsOrigin ? (
          <>
            {highlight(xMin, signal.ending, "left")}
            {highlight(signal.beginning, xMax, "right")}
          </>
        ) : (
          highlight(signal.beginning, signal.ending, "single")
        )}
        {visibleTraces.includes("pair12") ? <path d={path("pair12")} className="plot-line plot-pair-12" /> : null}
        {visibleTraces.includes("pair13") ? <path d={path("pair13")} className="plot-line plot-pair-13" /> : null}
        {visibleTraces.includes("pair23") ? <path d={path("pair23")} className="plot-line plot-pair-23" /> : null}
        <text x={margin.left} y={height - 12} className="plot-label">
          {xMin.toLocaleString()}
        </text>
        <text x={width - margin.right} y={height - 12} textAnchor="end" className="plot-label">
          {xMax.toLocaleString()} · alignment position
        </text>
      </svg>
      <div className="plot-legend" aria-hidden="true">
        {signal.method === "3SEQ" ? (
          <>
            <span className="legend-12">Candidate target {signal.tripletNames[0]}</span>
            <span className="legend-13">Candidate target {signal.tripletNames[1]}</span>
            <span className="legend-23">Candidate target {signal.tripletNames[2]}</span>
          </>
        ) : signal.method === "GENECONV" ? (
          <>
            <span className="legend-12">Inner {signal.tripletNames[0]}:{signal.tripletNames[1]} / outer {signal.tripletNames[2]}</span>
            <span className="legend-13">Inner {signal.tripletNames[0]}:{signal.tripletNames[2]} / outer {signal.tripletNames[1]}</span>
            <span className="legend-23">Inner {signal.tripletNames[1]}:{signal.tripletNames[2]} / outer {signal.tripletNames[0]}</span>
          </>
        ) : chimaeraTarget !== null && chimaeraParentOne !== null ? (
          <span className={chimaeraTrace === "pair12" ? "legend-12" : chimaeraTrace === "pair13" ? "legend-13" : "legend-23"}>
            Target {signal.tripletNames[chimaeraTarget]} : parent-one {signal.tripletNames[chimaeraParentOne]}
          </span>
        ) : (
          <>
            <span className="legend-12">{signal.tripletNames[0]} : {signal.tripletNames[1]}</span>
            <span className="legend-13">{signal.tripletNames[0]} : {signal.tripletNames[2]}</span>
            <span className="legend-23">{signal.tripletNames[1]} : {signal.tripletNames[2]}</span>
          </>
        )}
      </div>
      {!plot.detectionProfileExact ? (
        <p className="plot-context-note">
          Original-alignment reconstruction. The saved statistics and breakpoints retain the
          fragment/erasure-adjusted detection result; this compact checkpoint does not retain every
          historical working-profile point.
        </p>
      ) : null}
    </div>
  );
}
