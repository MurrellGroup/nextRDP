import { lstat, readFile, readdir } from "node:fs/promises";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";

const output = resolve(process.cwd(), "dist");
const requiredFiles = [
  "index.html",
  ".nojekyll",
  "wasm/rdp-core.mjs",
  "wasm/rdp-core.wasm",
];

function fail(message) {
  throw new Error(`GitHub Pages artifact check failed: ${message}`);
}

async function requireFile(relativePath) {
  const path = resolve(output, relativePath);
  const metadata = await lstat(path).catch(() => null);
  if (!metadata?.isFile() || (relativePath !== ".nojekyll" && metadata.size === 0)) {
    fail(`${relativePath} is missing or empty`);
  }
  return path;
}

let totalBytes = 0;
async function inspectTree(directory) {
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    const path = resolve(directory, entry.name);
    const metadata = await lstat(path);
    if (metadata.isSymbolicLink()) fail(`${path} is a symbolic link`);
    if (metadata.isDirectory()) {
      await inspectTree(path);
      continue;
    }
    if (!metadata.isFile()) fail(`${path} is not a regular file`);
    if (metadata.nlink > 1) fail(`${path} is hard-linked`);
    totalBytes += metadata.size;
  }
}

const [indexPath, , modulePath, wasmPath] = await Promise.all(
  requiredFiles.map(requireFile),
);

const html = await readFile(indexPath, "utf8");
if (html.includes("/src/main.tsx")) fail("index.html still references development source");
for (const match of html.matchAll(/\b(?:href|src)=["']([^"']+)["']/g)) {
  const url = match[1];
  if (url.startsWith("/") && !url.startsWith("//")) {
    fail(`index.html contains root-relative asset URL ${url}`);
  }
}

const moduleSource = await readFile(modulePath, "utf8");
if (!moduleSource.includes("rdp-core.wasm")) {
  fail("the Emscripten module does not reference rdp-core.wasm");
}

const wasm = await readFile(wasmPath);
if (wasm.length < 8 || !wasm.subarray(0, 4).equals(Buffer.from([0x00, 0x61, 0x73, 0x6d]))) {
  fail("rdp-core.wasm does not have a WebAssembly header");
}

const { default: createRdpModule } = await import(
  `${pathToFileURL(modulePath).href}?verify=${Date.now()}`,
);
const engine = await createRdpModule({
  noInitialRun: true,
  wasmBinary: wasm,
});
if (!(engine.HEAPU8 instanceof Uint8Array)) {
  fail("the Emscripten module does not expose its HEAPU8 upload view");
}

const context = engine._rdp_create();
if (!context) fail("the WASM engine could not create a FASTA smoke-test context");
const fasta = new TextEncoder().encode(
  ">alpha\nACGTACGT\n>beta\nACGTACGA\n>gamma\nACGTTCGA\n",
);
const fastaPointer = engine._malloc(fasta.byteLength);
if (!fastaPointer) {
  engine._rdp_destroy(context);
  fail("the WASM engine could not allocate the FASTA smoke-test input");
}
try {
  engine.HEAPU8.set(fasta, fastaPointer);
  if (engine._rdp_load_alignment(context, fastaPointer, fasta.byteLength) !== 1) {
    const detail = engine.UTF8ToString(engine._rdp_get_error(context));
    fail(`FASTA upload smoke test failed${detail ? `: ${detail}` : ""}`);
  }
  const summaryPointer = engine._rdp_get_summary_json(context);
  const summary = JSON.parse(engine.UTF8ToString(summaryPointer));
  if (summary.sequenceCount !== 3 || summary.alignmentLength !== 8) {
    fail("FASTA upload smoke test returned an unexpected alignment summary");
  }
} finally {
  engine._free(fastaPointer);
  engine._rdp_destroy(context);
}

await inspectTree(output);
console.log(
  `GitHub Pages artifact verified: ${requiredFiles.length} required files, FASTA upload smoke test passed, ${totalBytes.toLocaleString()} total bytes.`,
);
