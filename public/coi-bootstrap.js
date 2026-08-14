/*
 * Enables SharedArrayBuffer on static hosts such as GitHub Pages by letting a
 * same-origin service worker attach COOP/COEP headers to the next navigation.
 * If registration is unavailable or isolation still fails after one reload,
 * RDP Web simply keeps using its single-worker WASM module.
 */
(() => {
  "use strict";

  if (!("serviceWorker" in navigator) ||
      !window.isSecureContext ||
      location.protocol === "file:") {
    return;
  }

  const reloadKey = "rdp-web-coi-reload-v1";
  const reloadOnce = () => {
    if (window.crossOriginIsolated || sessionStorage.getItem(reloadKey) === "1") return;
    sessionStorage.setItem(reloadKey, "1");
    location.reload();
  };

  if (window.crossOriginIsolated) {
    sessionStorage.removeItem(reloadKey);
  }

  navigator.serviceWorker.register(
    new URL("coi-serviceworker.js", document.baseURI),
    { updateViaCache: "none" },
  ).then(async (registration) => {
    await registration.update().catch(() => undefined);
    await navigator.serviceWorker.ready;
    if (window.crossOriginIsolated) {
      sessionStorage.removeItem(reloadKey);
      return;
    }
    if (navigator.serviceWorker.controller) {
      reloadOnce();
      return;
    }
    navigator.serviceWorker.addEventListener("controllerchange", reloadOnce, { once: true });
  }).catch(() => {
    // A denied service worker is a supported configuration: the analysis
    // worker will load rdp-core.mjs and expose one CPU in Settings.
  });
})();
