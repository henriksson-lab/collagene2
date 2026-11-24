# Viewport performance measurements

The development-only `scripts/performance_browser_qa.py` generates independent synthetic GenBank fixtures, opens them in the actual release application, waits for browser persistence, then dispatches a rotating plasmid drag and scrolls across the entire sequence viewport. Each interaction is sampled for 120 animation frames. It writes raw environment, import/persistence timing, animation-frame intervals, dispatch timing, long-task observations and rendered-geometry counts as JSON. It invokes no scientific software and is not part of the application runtime.

## Local checkpoint, 2026-11-24

Environment: Linux; Intel Xeon Gold 6138 at 2 GHz; browser-reported 40 hardware threads; headless Chrome 104.0.5112.79; 1600 × 1000 viewport, device scale 1; ANGLE Vulkan **SwiftShader software renderer**. Host load and container CPU allocation were not controlled. These are short local measurements, not acceptance evidence for every browser, GPU, operating system or Tauri webview.

| Fixture | Interaction | Before p95 frame interval | After p95 frame interval |
| --- | --- | --- | --- |
| 20 kb, 200 features | Rotation | 16.666 ms | 16.666 ms |
| 20 kb, 200 features | Full-reference scrolling | 16.666 ms | 16.666 ms |
| 1 Mb, 10,000 features | Rotation | 16.666 ms | 16.666 ms |
| 1 Mb, 10,000 features | Full-reference scrolling | 33.332 ms | 16.666 ms |

No interaction long tasks over the observer threshold were reported in these runs. Import plus persistence took about 193/410 ms before and 203/464 ms after for the two fixtures; these single observations include startup/layout/persistence and should not be interpreted as isolated parser timings or a statistically established regression. The post-change renderer reports all 200/10,000 features included and no geometry-budget simplification. Before this change, map geometry was capped at an arbitrary first 2,000 features.

The scrolling change removes a full molecule/metadata clone on every repaint, builds a reusable bounded annotation interval index, and avoids resetting unchanged canvas dimensions. The map retains one document copy per revision and changes view uniforms during rotation. Its static geometry now has a 12 MiB vertex-data budget rather than a fixed feature prefix; a visible notice identifies simplification when that budget is reached. The interval index returns original-order visible tracks, including compound locations, with independent exhaustive query fixtures.

A final repeat after adding annotation pagination retained 16.666 ms p95 for all four interactions, included all 10,000 map features, and checked moving forward/backward through 500-entry annotation pages. Import/persistence was 110/489 ms in that repeat, illustrating why isolated startup observations should not be treated as a controlled parser comparison. [Raw observations](benchmarks/viewport-2026-11-24.json) retain all three runs and final source/release hashes. The baseline artifact digest was not captured; its missing source identity is explicit rather than inferred retroactively.

Frame intervals are animation-frame callback observations, not GPU timer queries or measured input-to-photon latency. Dispatch timings do not include all later reactive/raster work. Browser `performance.memory` is coarsened and does not establish WASM linear-memory usage; no memory-budget acceptance claim is made from it. This run does not cover continuous zoom, analysis-concurrent pointer latency, sustained GC behavior, the 10 Mb fixture, cold download/compile, hardware GPU execution, Firefox, Safari, or desktop webviews. Those gates remain open in the project plan.

To reproduce after starting the isolated preview/browser described in `scripts/browser_qa.py`:

```sh
python3 scripts/performance_browser_qa.py --output /tmp/collagene2-performance.json
```

The helper clears IndexedDB for its test origin; use a dedicated test profile. The fixture generator, measurement code and raw observations are retained in the repository for future comparisons. A CI runner should report its environment and measurements, not infer a universal performance guarantee from this table.

## Configurable cached-track checkpoint, 2026-11-24

After connecting the dedicated immutable track worker and configurable canvas tracks, the same local 120-frame fixture interactions measured 16.666 ms p95 for rotation and scrolling at both 20 kb/200 features and 1 Mb/10,000 features. All map features were retained without geometry simplification, and annotation pagination passed. Initial viewport track reports were complete. No interaction long tasks were observed. [Raw observations and source/release hashes](benchmarks/viewport-tracks-2026-11-24.json) retain this separate checkpoint.

Time to import, persist **and prepare the initial track report** was 314 ms and 1,023 ms in these individual observations. This adds track readiness to the measured interval; it is not a controlled comparison with the older import/persistence-only metric. Scrolling teleports across the full reference each frame; cached scientific items outside the prior tile refresh after debounce/settling. These frame intervals therefore establish neither complete track computation at every animation frame nor scientific analysis throughput. The hardware/browser/memory/10 Mb gates above remain open.
