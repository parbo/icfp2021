use eframe::{egui, epi};
//extern crate env_file;

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub type Point = egui::Pos2;

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

#[derive(Serialize, Deserialize, Debug)]
pub struct Pose {
    vertices: Vec<[i32; 2]>,
}

fn read_problem_from_file<P: AsRef<Path>>(path: P) -> Result<Problem, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `Problem`.
    let u = serde_json::from_reader(reader)?;

    // Return the `Problem`.
    Ok(u)
}

fn write_solution_to_file<P: AsRef<Path>>(path: P, verts: &[Point]) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let vertices: Vec<_> = verts.iter().map(|p| [p.x as i32, p.y as i32]).collect();
    let p = Pose { vertices };
    serde_json::to_writer(&file, &p)?;
    Ok(())
}

fn intersects(a: (Point, Point), b: (Point, Point)) -> bool {
    let Point { x: x1, y: y1 } = a.0;
    let Point { x: x2, y: y2 } = a.1;
    let Point { x: x3, y: y3 } = b.0;
    let Point { x: x4, y: y4 } = b.1;
    let tq = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4);
    let td = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if td != 0.0 {
        let t = tq as f32 / td as f32;
        if 0.0 < t && t < 1.0 {
            let uq = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3);
            let ud = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
            if ud != 0.0 {
                let u = uq as f32 / ud as f32;
                if 0.0 < u && u < 1.0 {
                    return true;
                }
            }
        }
    }
    false
}

fn inside(poly: &[egui::Pos2], p: egui::Pos2) -> bool {
    if poly.is_empty() {
        return false;
    }
    let mut i = 0;
    let mut j = poly.len() - 1;
    let mut c = false;
    while i < poly.len() {
        // If the point is in the polygon, then it is inside
        if p == poly[i] {
            return true;
        }
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
    c
}

fn place_vertice(
    problem: &Problem,
    memo: &mut HashMap<Vec<[i32; 2]>, Option<Vec<[i32; 2]>>>,
    verts: Vec<[i32; 2]>,
    min_x: i32,
    max_x: i32,
    min_y: i32,
    max_y: i32,
) -> Option<Vec<[i32; 2]>> {
    if let Some(p) = memo.get(&verts) {
        println!("memoized");
        return p.clone();
    }
    let ix = verts.len();
    if ix == problem.figure.vertices.len() {
        for (a, b) in &problem.figure.edges {
            let p1 = problem.figure.vertices[*a];
            let p2 = problem.figure.vertices[*b];
            let d = (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]);
            let pp1 = verts[*a];
            let pp2 = verts[*b];
            let dd = (pp1[0] - pp2[0]) * (pp1[0] - pp2[0]) + (pp1[1] - pp2[1]) * (pp1[1] - pp2[1]);
            let eps = 1000000.0 * ((d as f32 / dd as f32) - 1.0).abs();
            println!("{}, {}, {}, {}", d, dd, eps, problem.epsilon);
        }
        println!("all successfully placed! {:?}", verts);
        let v = verts.clone();
        memo.insert(verts, Some(v.clone()));
        return Some(v);
    }
    let hole = problem.hole.clone();
    for x in min_x..=max_x {
        for y in min_y..=max_y {
            let point = Point {
                x: x as f32,
                y: y as f32,
            };
            // Is the point inside?
            if !inside(&hole, point) {
                continue;
            }
            // Are all placed edges inside?
            // Are all placed edges within the constraint?
            let mut ok = true;
            for (a, b) in &problem.figure.edges {
                // Only check the new edges enabled by this point
                if (*a == ix || *b == ix) && *a <= ix && *b <= ix {
                    let p1 = problem.figure.vertices[*a];
                    let p2 = problem.figure.vertices[*b];
                    let d = (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]);
                    let pp1 = if *a == ix {
                        point
                    } else {
                        let [x, y] = verts[*a];
                        Point {
                            x: x as f32,
                            y: y as f32,
                        }
                    };
                    let pp2 = if *b == ix {
                        point
                    } else {
                        let [x, y] = verts[*b];
                        Point {
                            x: x as f32,
                            y: y as f32,
                        }
                    };
                    let dd = (pp1[0] - pp2[0]) * (pp1[0] - pp2[0])
                        + (pp1[1] - pp2[1]) * (pp1[1] - pp2[1]);
                    let eps = 1000000.0 * ((dd as f32 / d as f32) - 1.0).abs();
                    if eps as i32 > problem.epsilon {
                        ok = false;
                        break;
                    }
                    // check intersections
                    let mut i = 0;
                    let l = problem.hole.len();
                    while i < l {
                        let j = (i + 1) % l;
                        if intersects((pp1, pp2), (problem.hole[i], problem.hole[j])) {
                            ok = false;
                            break;
                        }
                        i += 1;
                    }
                    // check that the whole line is inside (there must be a fast and correct way to do this)
                    // This is not really 100% correct
                    for i in 0..1000 {
                        let dir = [pp2[0] - pp1[0], pp2[1] - pp1[1]];
                        let r = i as f32;
                        let pppp = Point {
                            x: pp1[0] + (r * dir[0]) / 1000.0,
                            y: pp1[1] + (r * dir[1]) / 1000.0,
                        };
                        if !inside(&problem.hole, pppp) {
                            ok = false;
                            break;
                        }
                    }
                }
            }
            if ok {
                let mut v = verts.clone();
                v.push([point.x as i32, point.y as i32]);
                if let Some(res) = place_vertice(problem, memo, v, min_x, max_x, min_y, max_y) {
                    memo.insert(verts, Some(res.clone()));
                    return Some(res);
                }
            }
        }
    }
    memo.insert(verts, None);
    None
}

fn solve(problem: &Problem) -> Option<Vec<[i32; 2]>> {
    let min_x = problem.hole.iter().map(|x| x[0] as i32).min().unwrap_or(0);
    let max_x = problem.hole.iter().map(|x| x[0] as i32).max().unwrap_or(0);
    let min_y = problem.hole.iter().map(|x| x[1] as i32).min().unwrap_or(0);
    let max_y = problem.hole.iter().map(|x| x[1] as i32).max().unwrap_or(0);
    let v = vec![];
    let mut m = HashMap::new();
    place_vertice(problem, &mut m, v, min_x, max_x, min_y, max_y)
}

/// We derive Deserialize/Serialize so we can persist app state on shutdown.
#[cfg_attr(feature = "persistence", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "persistence", serde(default))] // if we add new fields, give them default values when deserializing old state
pub struct PolygonApp {
    filename: String,
    grid: bool,

    #[cfg_attr(feature = "persistence", serde(skip))]
    problem: Option<Problem>,

    #[cfg_attr(feature = "persistence", serde(skip))]
    pose: Vec<Point>,

    #[cfg_attr(feature = "persistence", serde(skip))]
    selected: HashSet<usize>,
}

impl Default for PolygonApp {
    fn default() -> Self {
        Self {
            // Example stuff:
            filename: "problems/1.problem".to_owned(),
            grid: false,
            problem: None,
            pose: vec![],
            selected: HashSet::new(),
        }
    }
}

impl epi::App for PolygonApp {
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
            pose,
            selected,
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
                    *pose = read_problem.figure.vertices.clone();
                    *problem = Some(read_problem);
                    println!("Problem: {:?}", problem);
                }
            }
            ui.checkbox(grid, "Show Grid");
            ui.group(|ui| {
                ui.set_enabled(problem.is_some());
                if ui.button("Solve").clicked() {
                    if let Some(res) = solve(problem.as_ref().unwrap()) {
                        *pose = res
                            .into_iter()
                            .map(|p| Point {
                                x: p[0] as f32,
                                y: p[1] as f32,
                            })
                            .collect();
                    }
                }
                if ui.button("Save").clicked() {
                    let out_filename = filename.to_owned() + ".solution.json";
                    if write_solution_to_file(&out_filename, &pose).is_ok() {
                        println!("saved solution!");
                    }
                }
            });
            if let Some(p) = problem {
                if ui.button("Select All").clicked() {
                    *selected = HashSet::new();
                }
                egui::ScrollArea::auto_sized().show(ui, |ui| {
                    for i in 0..p.figure.edges.len() {
                        let checked = selected.contains(&i);
                        let (a, b) = p.figure.edges[i];
                        let p1 = p.figure.vertices[a];
                        let p2 = p.figure.vertices[b];
                        let d =
                            (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]);
                        let pp1 = pose[a];
                        let pp2 = pose[b];
                        let dd = (pp1[0] - pp2[0]) * (pp1[0] - pp2[0])
                            + (pp1[1] - pp2[1]) * (pp1[1] - pp2[1]);
                        let eps = 1000000.0 * ((dd as f32 / d as f32) - 1.0).abs();
                        let label = format!(
                            "{:?}, {}, {}, {}, {}",
                            p.figure.edges[i], d, dd, eps, p.epsilon
                        );
                        if ui.add(egui::SelectableLabel::new(checked, label)).clicked() {
                            if selected.contains(&i) {
                                selected.remove(&i);
                            } else {
                                selected.insert(i);
                            }
                        }
                    }
                });
            }
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
                    let min_x = problem.hole.iter().map(|x| x[0] as i32).min().unwrap_or(0) as f32;
                    let mut max_x =
                        problem.hole.iter().map(|x| x[0] as i32).max().unwrap_or(0) as f32;
                    let min_y = problem.hole.iter().map(|x| x[1] as i32).min().unwrap_or(0) as f32;
                    let mut max_y =
                        problem.hole.iter().map(|x| x[1] as i32).max().unwrap_or(0) as f32;
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
                    let points = problem.hole.iter().map(|x| to_screen * *x).collect();
                    shapes.push(egui::Shape::closed_line(points, hole_stroke));
                    // Draw pose
                    for ix in 0..problem.figure.edges.len() {
                        if !selected.is_empty() && !selected.contains(&ix) {
                            continue;
                        }
                        let edge = problem.figure.edges[ix];
                        let p1 = pose[edge.0];
                        let p2 = pose[edge.1];
                        let points = vec![to_screen * p1, to_screen * p2];
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
