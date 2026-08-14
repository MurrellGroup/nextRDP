#!/usr/bin/env bash
set -euo pipefail

binary="$(mktemp "${PWD}/wasm/.rdp-cyclic-pruning-core-check.XXXXXX")"
trap 'unlink "$binary" 2>/dev/null || true' EXIT

g++ -std=c++20 -O1 -Wall -Wextra -Wpedantic \
  -Iwasm/src \
  -Iwasm/include \
  scripts/verify-cyclic-pruning-core.cpp \
  wasm/src/alignment.cpp \
  wasm/src/bootscan.cpp \
  wasm/src/burt_confidence.cpp \
  wasm/src/geneconv.cpp \
  wasm/src/maxchi.cpp \
  wasm/src/phylpro.cpp \
  wasm/src/phylogeny.cpp \
  wasm/src/rdp_api.cpp \
  wasm/src/rdp_method.cpp \
  wasm/src/siscan.cpp \
  wasm/src/threeseq.cpp \
  -o "$binary"
RDP_DUMP_RESULTS=1 "$binary" | node scripts/verify-cyclic-pruning-digest.mjs
