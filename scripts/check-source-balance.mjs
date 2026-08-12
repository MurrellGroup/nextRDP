import { readFile, readdir } from "node:fs/promises";
import { extname, join, relative, resolve } from "node:path";

const root = process.cwd();
const roots = ["src", "wasm/include", "wasm/src"];
const extensions = new Set([".cpp", ".h", ".hpp", ".ts", ".tsx"]);

async function sourceFiles(directory) {
  const absolute = resolve(root, directory);
  const files = [];
  for (const entry of await readdir(absolute, { withFileTypes: true })) {
    const path = join(absolute, entry.name);
    if (entry.isDirectory()) files.push(...await sourceFiles(relative(root, path)));
    else if (entry.isFile() && extensions.has(extname(entry.name))) files.push(path);
  }
  return files;
}

function fail(path, line, column, message) {
  throw new Error(
    `Source balance check failed: ${relative(root, path)}:${line}:${column}: ${message}`,
  );
}

function checkBalanced(path, source) {
  const opening = new Set(["(", "[", "{"]);
  const expectedOpening = new Map([[")", "("], ["]", "["], ["}", "{"]]);
  const stack = [];
  let mode = "code";
  let quote = "";
  let escaped = false;
  let line = 1;
  let column = 0;

  for (let index = 0; index < source.length; ++index) {
    const character = source[index];
    const next = source[index + 1] ?? "";
    ++column;

    if (character === "\n") {
      ++line;
      column = 0;
      if (mode === "line-comment") mode = "code";
      if (mode === "string" && quote !== "`") {
        fail(path, line - 1, 1, "unterminated quoted string");
      }
      escaped = false;
      continue;
    }
    if (mode === "line-comment") continue;
    if (mode === "block-comment") {
      if (character === "*" && next === "/") {
        mode = "code";
        ++index;
        ++column;
      }
      continue;
    }
    if (mode === "string") {
      if (escaped) {
        escaped = false;
        continue;
      }
      if (character === "\\") {
        escaped = true;
        continue;
      }
      if (character === quote) {
        mode = "code";
        quote = "";
      }
      continue;
    }

    if (character === "/" && next === "/") {
      mode = "line-comment";
      ++index;
      ++column;
      continue;
    }
    if (character === "/" && next === "*") {
      mode = "block-comment";
      ++index;
      ++column;
      continue;
    }
    if (character === "\"" || character === "'" || character === "`") {
      mode = "string";
      quote = character;
      continue;
    }
    if (opening.has(character)) {
      stack.push({ character, line, column });
      continue;
    }
    if (expectedOpening.has(character)) {
      const found = stack.pop();
      const expected = expectedOpening.get(character);
      if (!found || found.character !== expected) {
        fail(path, line, column, `unexpected ${character}`);
      }
    }
  }

  if (mode === "block-comment") fail(path, line, column, "unterminated block comment");
  if (mode === "string") fail(path, line, column, "unterminated string or template literal");
  if (stack.length) {
    const found = stack.at(-1);
    fail(path, found.line, found.column, `unclosed ${found.character}`);
  }
}

const files = (await Promise.all(roots.map(sourceFiles))).flat().sort();
for (const path of files) checkBalanced(path, await readFile(path, "utf8"));
console.log(`Source delimiters verified across ${files.length} C++/TypeScript files.`);
