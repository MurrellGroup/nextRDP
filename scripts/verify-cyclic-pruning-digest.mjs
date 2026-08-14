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
const expected = "5ad90dbeeecd3ea531d52455dd3ded89498c8d0aeefc5d73c2885e451648e6fa";
if (digest !== expected) {
  throw new Error(`Cyclic pruning changed selected analytical results (${digest}).`);
}

const summary = input.split(/\r?\n/).find((line) =>
  line.startsWith("Cyclic pruning verified:"));
if (summary) console.log(summary);
console.log(`Cyclic selected-result digest verified: ${digest}.`);
