#!/usr/bin/env bash
set -euo pipefail

binary="$(mktemp "${PWD}/wasm/.rdp-bootscan-core-check.XXXXXX")"
trap 'rm -f "$binary"' EXIT

g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic \
  -Iwasm/src \
  -Iwasm/include \
  scripts/verify-bootscan-core.cpp \
  wasm/src/alignment.cpp \
  wasm/src/bootscan.cpp \
  wasm/src/burt_confidence.cpp \
  wasm/src/geneconv.cpp \
  wasm/src/maxchi.cpp \
  wasm/src/phylogeny.cpp \
  wasm/src/rdp_api.cpp \
  wasm/src/rdp_method.cpp \
  wasm/src/threeseq.cpp \
  -o "$binary"
"$binary"
