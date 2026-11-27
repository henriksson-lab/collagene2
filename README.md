# Collagene 2

A molecular cloning workspace in Rust, WebAssembly and Leptos, with WebGL2 maps and deterministic design verification. [MIT licensed](LICENSE). In active development.

![Sequence editor and circular plasmid map](docs/images/collagene-workspace.png)

![Vector assembly with selected backbone and insert](docs/images/collagene-assembly.png)

*Development screenshots, 27 November 2026. Assembly example uses synthetic sequences.*

## Install and start

With Rust installed:

```sh
rustup target add wasm32-unknown-unknown
cargo install trunk --locked
trunk serve
```

Open **http://127.0.0.1:8080**.

[Features, formats and platforms](docs/FEATURES.md) · [Application guide](docs/APPLICATION_GUIDE.md)

[Desktop installation/build](docs/TAURI_DESKTOP.md) · [Server setup](docs/SERVER_LIBRARY.md)

[Implementation status](docs/IMPLEMENTATION_STATUS.md) · [Product plan](docs/PROJECT_PLAN.md)
