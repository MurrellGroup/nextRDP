# Session 7 handoff

This checkpoint is source-only. No C++/Emscripten compilation, TypeScript check, Vite bundle,
preview, server, browser runtime, or project test was invoked. npm was used only with
`--package-lock-only --ignore-scripts`; it installed no dependencies and ran no package scripts.

## Completed in this phase

1. Restored the exact session 6 source checkpoint after workspace maintenance and verified its
   archive SHA-256 before editing.
2. Confirmed that session 6 had no `.github/workflows` path and therefore could not deploy when
   GitHub Pages was configured to use Actions.
3. Added `.github/workflows/deploy-pages.yml`, triggered by manual dispatch or a push to the
   repository's actual default branch. Feature-branch pushes are visible but skip the build.
4. Added the Pages permissions and environment contract: read-only repository contents,
   `pages: write`, `id-token: write`, serialized/cancellable Pages deployments, and the published
   URL from the deployment action.
5. Added a reproducible build environment: Node 20, npm lockfile restore, Emscripten 5.0.1, strict
   TypeScript checking, the existing release C++/WASM build, and the Vite production bundle.
6. Kept Pages on `rdp-core.mjs`/`rdp-core.wasm` only. Pages cannot emit cross-origin-isolation
   headers, so building the pthread variant there would add unusable files without accelerating
   the application. Numerical analysis still runs in its dedicated module worker.
7. Added `scripts/verify-pages-output.mjs` as a deployment gate. It requires nonempty HTML/loader,
   `.nojekyll`, a real WASM magic header, and the loader's external WASM reference; it rejects
   root-relative built assets, development entrypoints, symbolic links, and hard links.
8. Added a pre-build source-contract gate matching the C header, keepalive definitions, CMake
   exports, and worker calls, plus package/lock/native version and project-schema invariants.
9. Preserved hosting at `username.github.io`, `username.github.io/repository/`, and custom domains.
   Vite's relative base and the worker's `document.baseURI` resolution need no repository-name
   environment variable.
10. Added one package-version cache key to both the worker's dynamic `.mjs` import and Emscripten
   `locateFile` requests. Stable-name loader/WASM files therefore cannot be mixed across releases;
   the worker additionally destroys/refuses a context whose native engine version does not match.
11. Added `package-lock.json`, a Node `>=20` declaration, `npm ci` documentation, and a post-build
   verification script in the ordinary `npm run build` chain.
12. Added explicit project-checkpoint lifecycle state. A completed scan and every successful
    accept/reject, role/breakpoint edit, group edit, or downstream re-identification mark the
    project dirty; project import/download marks it current.
13. Added tab-close/reload protection while a scan, reconciliation, or unsaved completed analysis
    exists. Loading another dataset, changing mask/settings, or starting a replacement scan also
    confirms before discarding an uncheckpointed result.
14. Added visible dirty/current/saving state in the top bar and updated review/export checkpoint
    controls without adding browser persistence, analytics, or network transfer of sequence data.
15. Advanced source/engine labelling to `0.7.0-session-7`. Project schema remains `v1alpha7`
    because no serialized data contract changed.

## Deployment instructions

Place the contents of `rdp-wasm/` at the repository root. In GitHub, choose **Settings → Pages →
Source: GitHub Actions**, then push to the default branch. The workflow can also be run manually
from Actions. It publishes only the verified `dist/` artifact; no separate hosting configuration
or deployment branch is required.

## Fidelity boundary

The original RDP DLL/VB/manual attachments were not present after workspace maintenance and were
not included in the restored project ZIP. No new native algorithm was inferred or ported without
those references. The session 6 boundary therefore remains explicit: active `OKSeq` 0–14 (with
10/11 source-zero) and 17/18 are diagnostic, while final-list membership 15, later `FinalTrim`
expansion/special pruning, and the remaining `ConsensusOK` rebuild/straggler/score-and-prune path
are unported and cannot alter the automatic two-of-three membership.

## Static inspection performed

Only source and metadata were inspected. The workflow structure, action inputs, relative hosting
path, npm lock presence, version literals, React prop contracts, checkpoint mutation paths, source
delimiters, and archive contents were reviewed statically. The Pages workflow is intentionally the
first place compilation and type checking will occur when the user chooses to run it.

## Next fidelity phase

Once the supplied reference archives are attached again, resume from the post-`OKSeq 6`
`FinalTrim` expansions and selected-role special pruning that define `OKSeq 15`; then port the
remaining `ConsensusOK` list reconstruction and acceptance/pruning stages before changing automatic
membership. Separately, use the first authorized Actions run to resolve any compiler/type contract
diagnostics before native golden comparison.
