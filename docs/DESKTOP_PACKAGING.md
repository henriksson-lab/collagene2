# Desktop download artifacts

Automatic **Desktop build** CI produces a portable Windows executable, a macOS
DMG installer/application bundle, and a Linux Debian package. It builds the fresh
bounded Leptos frontend and scientific worker once, then calls the reusable
**Desktop packages (unsigned)** workflow with that exact same-run artifact. Each
native platform builds once; there is no duplicate native-only build.

The package workflow also supports manual dispatch: `development` is the default
for iteration, while `release` selects the repository's optimized native profile. Both
produce the same file formats. Neither publishes a GitHub release or uses signing
credentials.

| Runner | Download files | Actual payload check |
| --- | --- | --- |
| Ubuntu 22.04 x64 | `.deb` | Extract without installing; run its executable with `--mcp` |
| Windows hosted x64 | `collagene-desktop.exe` | Check native DLL imports, then run that exact executable with `--mcp` |
| macOS hosted architecture | `.dmg` and `.app.tar.gz` | Run the actual `.app/Contents/MacOS` executable with `--mcp` |

Windows uses `cargo tauri build --no-bundle --target x86_64-pc-windows-msvc`.
Target-scoped static-CRT flags and the pinned `tauri-build` implementation's
`STATIC_VCRUNTIME=true` avoid requiring an extra Visual C++ redistributable.
Explicitly selecting the target keeps those flags off host proc-macro builds.
CI inspects the resulting executable with the runner's MSVC `dumpbin`: dependency
imports must resolve to Windows system DLLs/API contracts, with versioned VC++
redistributable/debug CRT imports and an external `WebView2Loader.dll` rejected.
The target-specific output is checked and uploaded directly; no NSIS installer
is substituted for the portable executable.

The executable embeds the complete frontend and worker; it needs no adjacent
`dist` folder, development server or Node installation. It still requires the
Microsoft Edge WebView2 Runtime installed in Windows. Tauri's [official runtime
options](https://v2.tauri.app/distribute/windows-installer/#webview2-installation-options)
distinguish that OS runtime from application embedding. This plain executable
does not install or bundle a private browser runtime. Its native data root and
folder restrictions behave exactly as in the installed desktop application.

The macOS DMG contains the application for copying to Applications. The separate
tar archive preserves its executable bits and internal symlinks during artifact
upload. Native jobs select the runner's architecture; these are not universal
macOS or ARM Windows/Linux builds. The `.deb` uses Ubuntu 22.04 as its packaging
baseline and declares its operating-system runtime dependencies. See [Tauri's
Debian compatibility notes](https://v2.tauri.app/distribute/debian/#limitations).

Both automatic and manual builds use checksum-pinned upstream **cargo-tauri
2.11.3**, Trunk **0.21.14** and wasm-bindgen-cli **0.2.128** prebuilt tools, locked
Rust dependencies, two Cargo jobs, disabled incremental compilation and disabled
debug information. No missing CLI triggers a source-compilation fallback.
`tauri.conf.json` keeps ordinary bundling disabled; explicit platform bundle flags
select Debian or app/DMG output, while Windows uses `--no-bundle`. The packaging
overlay adds descriptions and the Education category. CI generates platform icon
layers from the existing provisional application PNG before compilation.

Download `desktop-packages-unsigned-<runner>-<profile>-<commit>` from the completed
**Desktop build** or manual package run. Artifacts expire after seven days. Each
contains exactly the requested platform files and `BUILD.json`: exact source SHA,
frontend artifact identity, profile, OS/architecture, runtime prerequisites and
per-file SHA-256 hashes. The workflow removes cached bundle output before each
build. Staging never merges a previous run's packages, and bounds output to 1 GiB
per file and 2 GiB per platform.

The headless MCP check imports synthetic FASTA and verifies a restriction digest.
It opens no GUI, browser, application HTTP listener, local repository or credential
vault. These checks establish build/artifact functionality; the user's review of
the final application remains the primary product acceptance step. Packages are
unsigned prerelease artifacts: no signing/notarization or updater publication is
claimed. Current workflow changes require a successful hosted run; no package
build or installation was performed locally.
