use eframe::{egui, epi};
//extern crate env_file;

use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub type Point = [i32; 2];

fn to_pos2(v: Point) -> egui::Pos2 {
    let x = v[0] as f32;
    let y = v[1] as f32;
    egui::Pos2 { x, y }
}

pub type Hole = Vec<Point>;

#[derive(Serialize, Deserialize, Debug)]
pub struct Figure {
    vertices: Vec<Point>,
    edges: Vec<(usize, usize)>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Problem {
    hole: Hole,
    figure: Figure,
    epsilon: i32,
}

fn read_problem_from_file<P: AsRef<Path>>(path: P) -> Result<Problem, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `Problem`.
    let u = serde_json::from_reader(reader)?;

    // Return the `User`.
    Ok(u)
}

/// We derive Deserialize/Serialize so we can persist app state on shutdown.
#[cfg_attr(feature = "persistence", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "persistence", serde(default))] // if we add new fields, give them default values when deserializing old state
pub struct TemplateApp {
    filename: String,
    grid: bool,

    #[cfg_attr(feature = "persistence", serde(skip))]
    problem: Option<Problem>,
}

impl Default for TemplateApp {
    fn default() -> Self {
        Self {
            // Example stuff:
            filename: "problems/1.problem".to_owned(),
            grid: false,
            problem: None,
        }
    }
}

impl epi::App for TemplateApp {
    fn name(&self) -> &str {
        "icfp 2021"
    }

    /// Called by the framework to load old app state (if any).
    #[cfg(feature = "persistence")]
    fn setup(
        &mut self,
        _ctx: &egui::CtxRef,
        _frame: &mut epi::Frame<'_>,
        storage: Option<&dyn epi::Storage>,
    ) {
        if let Some(storage) = storage {
            *self = epi::get_value(storage, epi::APP_KEY).unwrap_or_default()
        }
    }

    /// Called by the frame work to save state before shutdown.
    #[cfg(feature = "persistence")]
    fn save(&mut self, storage: &mut dyn epi::Storage) {
        epi::set_value(storage, epi::APP_KEY, self);
    }

    fn update(&mut self, ctx: &egui::CtxRef, frame: &mut epi::Frame<'_>) {
        let Self {
            filename,
            grid,
            problem,
        } = self;

        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            // The top panel is often a good place for a menu bar:
            egui::menu::bar(ui, |ui| {
                egui::menu::menu(ui, "File", |ui| {
                    if ui.button("Quit").clicked() {
                        frame.quit();
                    }
                });
            });
        });

        egui::SidePanel::left("side_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Filename: ");
                ui.text_edit_singleline(filename);
            });

            if ui.button("Load").clicked() {
                if let Ok(read_problem) = read_problem_from_file(&filename) {
                    *problem = Some(read_problem);
                    println!("Problem: {:?}", problem);
                }
            }
            ui.checkbox(grid, "Show Grid");
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::Frame::dark_canvas(ui.style()).show(ui, |ui| {
                ui.ctx().request_repaint();
                let desired_size = ui.available_size();
                if desired_size.x == 0.0 || desired_size.y == 0.0 {
                    return;
                }
                let (_id, rect) = ui.allocate_space(desired_size);
                ui.set_clip_rect(rect);

                let mut shapes = vec![];
                let hole_stroke = egui::Stroke::new(2.0, egui::Color32::LIGHT_BLUE);
                let pose_stroke = egui::Stroke::new(2.0, egui::Color32::RED);
                let grid_stroke = egui::Stroke::new(1.0, egui::Color32::LIGHT_GRAY);
                if let Some(problem) = problem {
                    // Transform so hole fits area
                    let min_x = problem.hole.iter().map(|x| x[0]).min().unwrap_or(0) as f32;
                    let mut max_x = problem.hole.iter().map(|x| x[0]).max().unwrap_or(0) as f32;
                    let min_y = problem.hole.iter().map(|x| x[1]).min().unwrap_or(0) as f32;
                    let mut max_y = problem.hole.iter().map(|x| x[1]).max().unwrap_or(0) as f32;
                    let width = max_x - min_x;
                    let height = max_y - min_y;
                    let x_ratio = desired_size.x / width;
                    let y_ratio = desired_size.y / height;
                    if y_ratio < x_ratio {
                        max_x = min_x + height * (desired_size.x / desired_size.y);
                    } else {
                        max_y = min_y + width * (desired_size.y / desired_size.x);
                    };
                    let to_screen = egui::emath::RectTransform::from_to(
                        egui::Rect::from_x_y_ranges(min_x..=max_x, min_y..=max_y),
                        rect,
                    );
                    if *grid {
                        // Draw grid
                        let mut x = min_x;
                        while x <= max_x {
                            let points = vec![
                                to_screen * egui::Pos2::from((x, min_y)),
                                to_screen * egui::Pos2::from((x, max_y)),
                            ];
                            shapes.push(egui::Shape::line(points, grid_stroke));
                            x += 10.0;
                        }
                        let mut y = min_y;
                        while y <= max_y {
                            let points = vec![
                                to_screen * egui::Pos2::from((min_x, y)),
                                to_screen * egui::Pos2::from((max_x, y)),
                            ];
                            shapes.push(egui::Shape::line(points, grid_stroke));
                            y += 10.0;
                        }
                    }
                    // Draw hole
                    let points = problem
                        .hole
                        .iter()
                        .map(|x| to_screen * to_pos2(*x))
                        .collect();
                    shapes.push(egui::Shape::closed_line(points, hole_stroke));
                    // Draw pose
                    for edge in &problem.figure.edges {
                        let points = vec![
                            to_screen * to_pos2(problem.figure.vertices[edge.0]),
                            to_screen * to_pos2(problem.figure.vertices[edge.1]),
                        ];
                        shapes.push(egui::Shape::line(points, pose_stroke));
                    }
                }
                ui.painter().extend(shapes);
            });
        });
    }
}
