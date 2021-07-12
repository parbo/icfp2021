#![forbid(unsafe_code)]
#![cfg_attr(not(debug_assertions), deny(warnings))] // Forbid warnings in release builds
#![warn(clippy::all, rust_2018_idioms)]

use clap::{App, Arg, SubCommand};
use dotenv::dotenv;

// When compiling natively:
#[cfg(not(target_arch = "wasm32"))]
fn main() {
    dotenv().ok();
    let matches = App::new("ICFP2021 Polygon")
        .subcommand(
            SubCommand::with_name("console").arg(
                Arg::with_name("filename")
                    .short("f")
                    .long("filename")
                    .value_name("FILE")
                    .required(true)
                    .takes_value(true),
            ),
        )
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("console") {
        let filename = matches.value_of("filename").unwrap();
        let mut app = polygon::PolygonApp::default();
        app.solve_console(filename);
    } else {
        let app = polygon::PolygonApp::default();
        let native_options = eframe::NativeOptions::default();
        eframe::run_native(Box::new(app), native_options);
    }
}
