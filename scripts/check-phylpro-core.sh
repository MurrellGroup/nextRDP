#!/usr/bin/env bash
set -euo pipefail

binary="$(mktemp "${PWD}/wasm/.rdp-phylpro-core-check.XXXXXX")"
trap 'unlink "$binary" 2>/dev/null || true' EXIT

g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic \
  -Iwasm/src \
  scripts/verify-phylpro-core.cpp \
  wasm/src/alignment.cpp \
  wasm/src/phylpro.cpp \
  -o "$binary"
"$binary"
