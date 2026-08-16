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
// This digest includes DefineEventP2's decrement-before-check candidate-run
// anchoring. Peaks immediately following a supporting run must attach to the
// preceding run rather than jump forward to the next one.
const expected = "5fd2fffa2413e4aac7ccce896fa9b72d88330e8968a52c985be5903d7766950c";
if (digest !== expected) {
  throw new Error(`Cyclic pruning changed selected analytical results (${digest}).`);
}

const summary = input.split(/\r?\n/).find((line) =>
  line.startsWith("Cyclic pruning verified:"));
if (summary) console.log(summary);
console.log(`Cyclic selected-result digest verified: ${digest}.`);
