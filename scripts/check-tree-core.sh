#!/usr/bin/env bash
set -euo pipefail

binary="$(mktemp "${PWD}/wasm/.rdp-tree-core-check.XXXXXX")"
trap 'rm -f "$binary"' EXIT

g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic \
  -Iwasm/src \
  scripts/verify-tree-core.cpp \
  wasm/src/alignment.cpp \
  wasm/src/phylogeny.cpp \
  -o "$binary"
"$binary"
