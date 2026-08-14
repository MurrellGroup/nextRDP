# Session 24 handoff — timed deterministic multicore analysis

Version: `0.24.0-session-24`
Project schema: `org.rdp-web.project/v1alpha19`

Session 24 adds runtime visibility and CPU parallelism without changing the analytical scheduler.
The supplied RDP sources and manual remained read-only and were not compiled.
No alternate RDP implementation was consulted.

## Timing and progress

- A worker-side monotonic clock based on `performance.now()` measures the complete request.
- The result separates setup, opening primary scan, later cyclic rescans, final reconciliation,
  current round, and every completed round. The final display retains millisecond precision.
- Timing is attached to live progress, final `ScanResults`, and the downloaded project analysis.
  Older project schemas remain importable and simply lack historical timing.
- Routine WASM JSON construction, structured-clone posting, and React rendering are limited to one
  update per 100 ms. Forced phase and terminal updates are never suppressed.
- At the zero-triplet boundary of a new cyclic pass, an indeterminate Windows-style marquee keeps
  the progress control visible until determinate progress resumes.

## Multicore boundary

- A persistent `std::thread` executor is compiled only for the pthread WASM target.
- For one triplet, enabled GENECONV, BootScan, MaxChi, CHIMAERA, SISCAN, and 3SEQ method kernels may
  execute concurrently. Each owns its existing method workspace, counters, candidate vector, and
  output vector.
- The RDP profile is prepared once before dispatch. Immutable coordinates, equality tracks,
  similarities, missing-data flags, options, and working alignment are shared read-only.
- Triplet traversal, invalid-schedule filtering, XOverList/BestXOList-style reuse, clean-triplet
  pruning, cyclic event selection, erasure, and cancellation stay on the analysis worker in their
  original serial order.
- Per-method outputs are merged only after all tasks finish and always in the existing method-major
  order. This preserves signal IDs, exact-P tie behavior, cache state, and result serialization.
- The linked `check:multicore-core` fixture enables five discovery families and requires exact
  1-CPU/4-CPU progress and complete results. Its wall-clock figures are informational, not a flaky
  performance threshold.

The browser reads `navigator.hardwareConcurrency`, leaves roughly one quarter of logical CPUs as
headroom, caps the pool at six total CPUs (the number of independent heavy method lanes), and
exposes a 1–N selector. A user can always select one CPU. The active pool is further capped by the
number of heavy methods enabled for discovery. This is method-level parallelism: data sets with
longer alignments and multiple expensive methods benefit most, while a lightweight or single-method
scan may remain synchronization-bound.

## Static hosting and GitHub Pages

Emscripten pthreads require `SharedArrayBuffer` and a cross-origin-isolated page. The production
build now emits both `rdp-core-threads.*` and `rdp-core.*`. Vite development/preview sends COOP and
COEP directly. GitHub Pages cannot configure those headers, so the static bundle registers a
same-origin service worker that adds them to fetched responses and reloads the page once after it
takes control. It does not cache sequence or application data. If service workers or isolation are
unavailable, the worker loads the single-worker module and Settings exposes one CPU.

The Pages verifier now requires both modules, the generated pthread worker helper, bootstrap and
service-worker assets, valid WebAssembly headers, relative project-site URLs, the compatibility ABI
smoke test, and the existing FASTA/cyclic/BootScan/SISCAN/PHYLPRO/graceful-stop checks. GitHub
Actions also runs the deterministic multicore host gate before building the artifact.

## Validation performed in this checkpoint

- Strict TypeScript contracts and production Vite UI build.
- All port translation units linked with `-pthread -DRDP_ENABLE_THREADS=1` in the multicore gate.
- Exact 1-CPU/4-CPU progress and result JSON equality.
- A full ThreadSanitizer run of the heavier one/four-CPU fixture completed without a reported race.
- Existing cyclic selected-result SHA-256 digest unchanged.
- Existing BootScan, SISCAN, event-tree, PHYLPRO, source-contract, and source-balance gates.

The local environment does not provide `emcmake`; therefore the actual dual Emscripten build and
instantiated Pages artifact remain the explicit GitHub Actions gate. This is not a native golden
parity claim.
