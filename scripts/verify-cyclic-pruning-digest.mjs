import { createHash } from "node:crypto";

let input = "";
for await (const chunk of process.stdin) input += chunk;

const prefix = "RDP_RESULTS_JSON=";
const resultLine = input.split(/\r?\n/).find((line) => line.startsWith(prefix));
if (!resultLine) throw new Error("The cyclic core regression did not emit its result JSON.");

const results = JSON.parse(resultLine.slice(prefix.length));
const selectedResult = {
  events: results.events.map((event) => [
    event.recombinant,
    event.majorParent,
    event.minorParent,
    event.beginning,
    event.ending,
    event.wrapsOrigin,
    event.detectionRound,
    event.supportSignalIds,
  ]),
  signals: results.signals.map((signal) => [
    signal.method,
    signal.triplet,
    signal.recombinant,
    signal.majorParent,
    signal.minorParent,
    signal.beginning,
    signal.ending,
    signal.correctedPValue,
    signal.eventId,
  ]),
};
const digest = createHash("sha256")
  .update(JSON.stringify(selectedResult))
  .digest("hex");
// The native-parity checkpoint adds the exact Clearcut/Tree2ArrayP2 event-tree
// path and the currently mapped MakeConsensusC statistics before this
// deterministic snapshot is taken.
const expected = "5a68fbe96d81e69e1ae1a837dd9cc38f059f761385a2ed706f29f385305a3232";
if (digest !== expected) {
  throw new Error(`Cyclic pruning changed selected analytical results (${digest}).`);
}

const summary = input.split(/\r?\n/).find((line) =>
  line.startsWith("Cyclic pruning verified:"));
if (summary) console.log(summary);
console.log(`Cyclic selected-result digest verified: ${digest}.`);
