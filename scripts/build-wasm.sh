#!/usr/bin/env bash
set -euo pipefail

mode="${1:-single}"
case "$mode" in
  single)
    thread_option="OFF"
    build_directory="wasm/build-single"
    ;;
  threads)
    thread_option="ON"
    build_directory="wasm/build-threads"
    ;;
  *)
    echo "Usage: $0 [single|threads]" >&2
    exit 2
    ;;
esac

if ! command -v emcmake >/dev/null 2>&1; then
  echo "Emscripten is required (emcmake was not found)." >&2
  exit 1
fi

emcmake cmake -S wasm -B "$build_directory" \
  -DCMAKE_BUILD_TYPE=Release \
  -DRDP_THREADS="$thread_option"
cmake --build "$build_directory" --config Release --parallel
