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

fn intersects(a: (Point, Point), b: (Point, Point)) -> bool {
    let [x1, y1] = a.0;
    let [x2, y2] = a.1;
    let [x3, y3] = b.0;
    let [x4, y4] = b.1;
    let tq = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4);
    let td = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if td != 0 {
        let t = tq as f32 / td as f32;
        if 0.0 < t && t <= 1.0 {
            let uq = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3);
            let ud = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
            if ud != 0 {
                let u = uq as f32 / ud as f32;
                if 0.0 < u && u <= 1.0 {
                    return true;
                }
            }
        }
    }
    return false;
}

fn inside(poly: &[Point], p: Point) -> bool {
    if poly.len() == 0 {
        return false;
    }
    let mut i = 0;
    let mut j = poly.len() - 1;
    let mut c = false;
    while i < poly.len() {
        if (((poly[i][1] <= p[1]) && (p[1] < poly[j][1]))
            || ((poly[j][1] <= p[1]) && (p[1] < poly[i][1])))
            && (p[0]
                < (poly[j][0] - poly[i][0]) * (p[1] - poly[i][1]) / (poly[j][1] - poly[i][1])
                    + poly[i][0])
        {
            c = !c;
        }
        j = i;
        i += 1;
    }
    return c;
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
                let pose_stroke = egui::Stroke::new(2.0, egui::Color32::GREEN);
                let intersect_stroke = egui::Stroke::new(2.0, egui::Color32::RED);
                let outside_stroke = egui::Stroke::new(2.0, egui::Color32::YELLOW);
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
                        let p1 = problem.figure.vertices[edge.0];
                        let p2 = problem.figure.vertices[edge.1];
                        let points = vec![to_screen * to_pos2(p1), to_screen * to_pos2(p2)];
                        let mut stroke = if !inside(&problem.hole, p1) || !inside(&problem.hole, p2)
                        {
                            outside_stroke
                        } else {
                            pose_stroke
                        };
                        let mut i = 0;
                        let l = problem.hole.len();
                        while i < l {
                            let j = (i + 1) % l;
                            if intersects((p1, p2), (problem.hole[i], problem.hole[j])) {
                                stroke = intersect_stroke;
                                break;
                            }
                            i += 1;
                        }
                        shapes.push(egui::Shape::line(points, stroke));
                    }
                }
                ui.painter().extend(shapes);
            });
        });
    }
}
