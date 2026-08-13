import { GitBranch } from "lucide-react";
import { useMemo, useState } from "react";

import type {
  EventAlignmentRole,
  EventTreeEdge,
  EventTreeLeaf,
  EventTreeRegion,
  EventTreeView,
  ReconciledEvent,
} from "../lib/types";

interface EventTreeInspectorProps {
  event: ReconciledEvent;
  view: EventTreeView | null;
  loading: boolean;
  error: string;
  onClose: () => void;
}

type Comparison = "whole-event" | "beginning" | "ending";

const comparisons: Record<
  Comparison,
  { label: string; detail: string; regions: [EventTreeRegion["name"], EventTreeRegion["name"]] }
> = {
  "whole-event": {
    label: "Whole event",
    detail: "outside tract ↔ inside tract",
    regions: ["outside-tract", "inside-tract"],
  },
  beginning: {
    label: "5′ boundary",
    detail: "before ↔ after beginning",
    regions: ["5-prime-outside", "5-prime-inside"],
  },
  ending: {
    label: "3′ boundary",
    detail: "before ↔ after ending",
    regions: ["3-prime-inside", "3-prime-outside"],
  },
};

const roleLabels: Record<EventAlignmentRole, string> = {
  recombinant: "recombinant",
  "major-parent": "major parent",
  "minor-parent": "minor parent",
  "co-recombinant": "co-recombinant",
  evidence: "tree-panel evidence",
};

function compactName(name: string, maximum = 34): string {
  return name.length <= maximum ? name : `${name.slice(0, maximum - 1)}…`;
}

interface DirectedTreeEdge {
  parent: number;
  child: number;
  edge: EventTreeEdge;
}

interface TreeLayout {
  height: number;
  positions: Map<number, { x: number; y: number }>;
  directedEdges: DirectedTreeEdge[];
}

function buildTreeLayout(
  region: EventTreeRegion,
  leaves: EventTreeLeaf[],
  collapseWeak: boolean,
): TreeLayout | null {
  if (!region.usable || region.nodeCount === 0 || region.edges.length === 0) return null;
  const adjacency = Array.from({ length: region.nodeCount }, () => [] as Array<{
    node: number;
    edge: EventTreeEdge;
  }>);
  region.edges.forEach((edge) => {
    if (
      edge.from < 0 ||
      edge.to < 0 ||
      edge.from >= region.nodeCount ||
      edge.to >= region.nodeCount
    ) return;
    adjacency[edge.from].push({ node: edge.to, edge });
    adjacency[edge.to].push({ node: edge.from, edge });
  });

  const root = region.root < region.nodeCount
    ? region.root
    : region.edges.find((edge) => edge.internal)?.from ?? 0;
  const parent = new Array<number>(region.nodeCount).fill(-2);
  const children = Array.from({ length: region.nodeCount }, () => [] as number[]);
  const edgeToParent = new Array<EventTreeEdge | null>(region.nodeCount).fill(null);
  parent[root] = -1;
  const discover = (node: number) => {
    adjacency[node]
      .slice()
      .sort((left, right) => left.node - right.node)
      .forEach(({ node: next, edge }) => {
        if (next === parent[node] || parent[next] !== -2) return;
        parent[next] = node;
        edgeToParent[next] = edge;
        children[node].push(next);
        discover(next);
      });
  };
  discover(root);

  const leafNodes = new Set(leaves.map((leaf) => leaf.node));
  const minimumLeaf = new Array<number | null>(region.nodeCount).fill(null);
  const descendantMinimum = (node: number): number => {
    if (minimumLeaf[node] !== null) return minimumLeaf[node]!;
    const value = leafNodes.has(node)
      ? node
      : children[node].reduce(
          (minimum, child) => Math.min(minimum, descendantMinimum(child)),
          Number.MAX_SAFE_INTEGER,
        );
    minimumLeaf[node] = value;
    return value;
  };
  descendantMinimum(root);
  children.forEach((nodes) => {
    nodes.sort((left, right) => descendantMinimum(left) - descendantMinimum(right));
  });

  const orderedLeaves: number[] = [];
  const collectLeaves = (node: number) => {
    if (leafNodes.has(node)) {
      orderedLeaves.push(node);
      return;
    }
    children[node].forEach(collectLeaves);
  };
  collectLeaves(root);
  const leafSpacing = 24;
  const top = 25;
  const height = Math.max(180, orderedLeaves.length * leafSpacing + 34);
  const y = new Array<number>(region.nodeCount).fill(top);
  orderedLeaves.forEach((node, index) => {
    y[node] = top + index * leafSpacing;
  });
  const placeInternal = (node: number): number => {
    if (leafNodes.has(node) || children[node].length === 0) return y[node];
    const childPositions = children[node].map(placeInternal);
    y[node] = childPositions.reduce((sum, value) => sum + value, 0) / childPositions.length;
    return y[node];
  };
  placeInternal(root);

  const distance = new Array<number>(region.nodeCount).fill(0);
  const directedEdges: DirectedTreeEdge[] = [];
  const measure = (node: number) => {
    children[node].forEach((child) => {
      const edge = edgeToParent[child]!;
      const length = collapseWeak && edge.collapsed ? 0 : Math.max(0, edge.length);
      distance[child] = distance[node] + length;
      directedEdges.push({ parent: node, child, edge });
      measure(child);
    });
  };
  measure(root);
  const maximumDistance = distance.reduce((maximum, value) => Math.max(maximum, value), 0);
  const scale = maximumDistance > 0 ? 485 / maximumDistance : 1;
  const positions = new Map<number, { x: number; y: number }>();
  distance.forEach((value, node) => {
    if (parent[node] !== -2) positions.set(node, { x: 30 + value * scale, y: y[node] });
  });
  return { height, positions, directedEdges };
}

function TreeFigure({
  region,
  leaves,
  collapseWeak,
}: {
  region: EventTreeRegion | undefined;
  leaves: EventTreeLeaf[];
  collapseWeak: boolean;
}) {
  const layout = useMemo(
    () => (region ? buildTreeLayout(region, leaves, collapseWeak) : null),
    [collapseWeak, leaves, region],
  );
  const leavesByNode = useMemo(
    () => new Map(leaves.map((leaf) => [leaf.node, leaf])),
    [leaves],
  );

  if (!region || !layout) {
    return (
      <section className="event-tree-figure is-unavailable">
        <div>
          <span>{region?.name.replaceAll("-", " ") ?? "regional tree"}</span>
          <strong>Tree unavailable</strong>
          <small>
            {region ? `${region.sites} usable region sites` : "No saved region metadata"} ·
            Jukes–Cantor fallback remains in the numeric evidence
          </small>
        </div>
      </section>
    );
  }

  return (
    <section className="event-tree-figure">
      <header>
        <div>
          <span>{region.name.replaceAll("-", " ")}</span>
          <strong>{region.sites.toLocaleString()} sites</strong>
        </div>
        <small>
          {region.supportedInternalBranches}/{region.internalBranches} internal branches ≥50% ·{" "}
          {region.collapsedDistanceRankLevels} distance ranks
        </small>
      </header>
      <div className="event-tree-scroll" tabIndex={0}>
        <svg
          className="event-tree-svg"
          viewBox={`0 0 930 ${layout.height}`}
          style={{ minWidth: 760, height: layout.height }}
          role="img"
          aria-label={`${region.name.replaceAll("-", " ")} neighbour-joining tree`}
        >
          {layout.directedEdges.map(({ parent, child, edge }) => {
            const start = layout.positions.get(parent);
            const end = layout.positions.get(child);
            if (!start || !end) return null;
            const weak = edge.collapsed;
            return (
              <g key={`${parent}-${child}`}>
                <path
                  className={`event-tree-edge${weak ? " is-weak" : ""}${
                    weak && collapseWeak ? " is-collapsed" : ""
                  }`}
                  d={`M ${start.x} ${start.y} V ${end.y} H ${end.x}`}
                />
                {edge.bootstrapSupport !== null ? (
                  <text
                    className={`event-tree-support${weak ? " is-weak" : ""}`}
                    x={Math.max(start.x, end.x - 23)}
                    y={end.y - 4}
                  >
                    {Math.round(edge.bootstrapSupport * 100)}
                  </text>
                ) : null}
              </g>
            );
          })}
          {[...layout.positions.entries()].map(([node, position]) => {
            const leaf = leavesByNode.get(node);
            if (!leaf) return null;
            const suffix = leaf.fragmentEventId === null
              ? ""
              : ` · retained fragment from event ${leaf.fragmentEventId + 1}`;
            const inputRole = leaf.queryReferenceInputRole === "not-applied"
              ? ""
              : leaf.queryReferenceInputRole === "reference"
                ? ` · reference group ${leaf.referenceGroup ?? "?"}`
                : " · query input";
            return (
              <g className={`event-tree-leaf role-${leaf.role}`} key={node}>
                <circle cx={position.x} cy={position.y} r={3.4}>
                  <title>{`${leaf.sequenceName} · ${roleLabels[leaf.role]}${inputRole}${suffix}`}</title>
                </circle>
                <line x1={position.x + 4} y1={position.y} x2={538} y2={position.y} />
                <text x={545} y={position.y + 3.5}>
                  {compactName(leaf.sequenceName)}
                  {leaf.fragmentEventId === null ? "" : ` · fragment e${leaf.fragmentEventId + 1}`}
                </text>
                <text className="event-tree-leaf-meta" x={790} y={position.y + 3.5}>
                  {roleLabels[leaf.role]}
                  {inputRole}
                  {leaf.currentGroupMember ? " · current group" : ""}
                  {leaf.trace ? " · masked trace" : leaf.masked ? " · masked" : ""}
                  {leaf.disabled ? " · disabled tree context" : ""}
                </text>
              </g>
            );
          })}
        </svg>
      </div>
    </section>
  );
}

export function EventTreeInspector({
  event,
  view,
  loading,
  error,
  onClose,
}: EventTreeInspectorProps) {
  const [comparison, setComparison] = useState<Comparison>("whole-event");
  const [collapseWeak, setCollapseWeak] = useState(true);
  const selectedComparison = comparisons[comparison];
  const selectedRegions = selectedComparison.regions.map((name) =>
    view?.regions.find((region) => region.name === name),
  );

  return (
    <section className="event-tree-card">
      <div className="card-heading split-heading">
        <div>
          <span className="eyebrow">Graphical phylogenetic check</span>
          <h3>Regional neighbour-joining topologies</h3>
          <p>
            Compare the saved topology immediately outside and inside the tract or either
            boundary for event {event.id + 1}. Trees are arbitrarily rooted for display only.
          </p>
        </div>
        <button className="button button-quiet button-compact" type="button" onClick={onClose}>
          Close
        </button>
      </div>

      {loading ? (
        <div className="tree-loading" role="status">
          <span /> Preparing the saved regional topologies…
        </div>
      ) : error ? (
        <div className="alignment-error" role="alert">{error}</div>
      ) : view ? (
        <>
          <div className="tree-inspector-toolbar">
            <div className="tree-comparison-tabs" aria-label="Tree region comparison">
              {(Object.keys(comparisons) as Comparison[]).map((key) => (
                <button
                  className={comparison === key ? "is-selected" : ""}
                  type="button"
                  onClick={() => setComparison(key)}
                  key={key}
                >
                  <strong>{comparisons[key].label}</strong>
                  <small>{comparisons[key].detail}</small>
                </button>
              ))}
            </div>
            <label className="tree-collapse-toggle">
              <input
                type="checkbox"
                checked={collapseWeak}
                onChange={(change) => setCollapseWeak(change.target.checked)}
              />
              <span>Collapse branches below 50% bootstrap support</span>
            </label>
          </div>
          <div className="tree-inspector-legend">
            <span className="role-recombinant">Recombinant</span>
            <span className="role-major-parent">Major parent</span>
            <span className="role-minor-parent">Minor parent</span>
            <span className="role-co-recombinant">Co-recombinant</span>
            <span className="tree-weak-key">Bootstrap &lt;50%</span>
          </div>
          <div className="event-tree-pair">
            <TreeFigure region={selectedRegions[0]} leaves={view.leaves} collapseWeak={collapseWeak} />
            <TreeFigure region={selectedRegions[1]} leaves={view.leaves} collapseWeak={collapseWeak} />
          </div>
          <p className="tree-inspector-footnote">
            <GitBranch size={13} /> Supplied Clearcut float NJ · {view.bootstrapReplicates}
            {" "}Microsoft-CRT SEQBOOT2 replicates (seed {view.randomSeed}) · branch labels use the
            desktop base-tree pseudocount percentages. Analysis uses the midpoint-rooted ultrametric
            rank matrices; this drawing uses the saved five-decimal branch lengths. Opening this
            panel transfers only the bounded edge list already produced for reconciliation, not
            distance matrices or alignment rows{view.subsampled
              ? ` · the closest ${view.sequenceCap} working sequences form this panel`
              : ""}{view.fragmentAssisted
              ? " · retained-fragment leaves are labelled with their source event"
              : ""} · breakpoint flanks target {view.flankVariableSiteTarget} informative variable sites.
          </p>
        </>
      ) : null}
    </section>
  );
}
