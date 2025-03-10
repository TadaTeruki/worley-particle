extern crate flatc_rust;
use std::path::Path;

fn main() {
    flatc_rust::run(flatc_rust::Args {
        inputs: &[Path::new("schema/particlemap.fbs")],
        out_dir: Path::new("src/map/rw"),
        ..Default::default()
    })
    .expect("flatc failed to run");
}
