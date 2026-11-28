pub(crate) mod ai;
#[cfg(target_arch = "wasm32")]
pub(crate) mod app;
pub(crate) mod assembly;
pub(crate) mod biology;
pub(crate) mod cloning;
pub(crate) mod genomes;
pub(crate) mod history;
pub(crate) mod io;
pub(crate) mod library;
pub(crate) mod platform;
pub(crate) mod primers;
pub(crate) mod protocol;
pub(crate) mod reads;
pub(crate) mod restriction;
pub(crate) mod sources;
pub(crate) mod viewer;
pub(crate) mod workspace;

#[cfg(target_arch = "wasm32")]
fn main() {
    leptos::mount::mount_to_body(app::App);
}

#[cfg(not(target_arch = "wasm32"))]
fn main() {
    println!("Build the web application with: trunk serve crates/web/index.html");
}
