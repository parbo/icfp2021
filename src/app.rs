use eframe::{egui, epi};

use itertools::Itertools;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;

pub type Point = egui::Pos2;

pub type Hole = Vec<Point>;

#[derive(Serialize, Deserialize, Debug, Clone, Eq, PartialEq)]
pub struct Figure {
    vertices: Vec<Point>,
    edges: Vec<(usize, usize)>,
}

#[derive(Serialize, Deserialize, Debug, Clone, Eq, PartialEq)]
pub struct Bonus {
    position: Point,
    bonus: String,
    problem: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone, Eq, PartialEq)]
pub struct Problem {
    hole: Hole,
    figure: Figure,
    epsilon: i32,
    bonuses: Vec<Bonus>,
}

#[derive(Serialize, Deserialize, Debug, Clone, Eq, PartialEq)]
pub struct Pose {
    vertices: Vec<[i32; 2]>,
    bonuses: Vec<Bonus>,
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
    let p = Pose {
        vertices,
        bonuses: vec![],
    };
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
        if (0.0..=1.0).contains(&t) {
            let uq = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3);
            let ud = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
            if ud != 0.0 {
                let u = uq as f32 / ud as f32;
                if (0.0..=1.0).contains(&u) {
                    // Check if only touching
                    return !(((u - 0.0).abs() < f32::EPSILON || (u - 1.0).abs() < f32::EPSILON)
                        && ((t - 0.0).abs() < f32::EPSILON || (t - 1.0).abs() < f32::EPSILON));
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

#[derive(Debug, Clone, Eq, PartialEq)]
struct Thing {
    num: usize,
    verts: HashMap<usize, [i32; 2]>,
    dislikes: i32,
    constrainedness: Vec<(usize, i32)>,
    edge_length: i32,
}

impl Ord for Thing {
    fn cmp(&self, other: &Thing) -> Ordering {
        self.verts
            .len()
            .cmp(&other.verts.len())
            .then_with(|| self.edge_length.cmp(&other.edge_length))
            .then_with(|| other.dislikes.cmp(&self.dislikes))
    }
}

// `PartialOrd` needs to be implemented as well.
impl PartialOrd for Thing {
    fn partial_cmp(&self, other: &Thing) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn calc_constrainedness(problem: &Problem, pose: &HashMap<usize, [i32; 2]>) -> Vec<(usize, i32)> {
    let mut cnt = HashMap::new();
    for (a, b) in &problem.figure.edges {
        if !pose.contains_key(a) {
            let inc = if pose.contains_key(b) { 10 } else { 1 };
            *cnt.entry(*a).or_insert(0) += inc;
        }
        if !pose.contains_key(b) {
            let inc = if pose.contains_key(a) { 10 } else { 1 };
            *cnt.entry(*b).or_insert(0) += inc;
        }
    }
    let mut v: Vec<_> = cnt.into_iter().collect();
    v.sort_by(|a, b| b.1.cmp(&a.1));
    v
}

fn calc_edge_length(problem: &Problem, pose: &HashMap<usize, [i32; 2]>) -> i32 {
    let mut tot = 0;
    for (a, b) in &problem.figure.edges {
        if pose.contains_key(a) && pose.contains_key(b) {
            let pa = problem.figure.vertices[*a];
            let pb = problem.figure.vertices[*b];
            let edge_length =
                ((pa.x - pb.x) * (pa.x - pb.x) + (pa.y - pb.y) * (pa.y - pb.y)) as i32;
            tot += edge_length;
        }
    }
    tot
}

fn calc_dislikes(problem: &Problem, pose: &HashMap<usize, [i32; 2]>) -> i32 {
    let mut d_sum = 0;
    for v in &problem.hole {
        let mut min_d = None;
        for p in pose.values() {
            let d = (v.x as i32 - p[0]) * (v.x as i32 - p[0])
                + (v.y as i32 - p[1]) * (v.y as i32 - p[1]);
            if let Some(md) = min_d {
                if d < md {
                    min_d = Some(d);
                }
            } else {
                min_d = Some(d);
            }
        }
        d_sum += min_d.unwrap_or(0);
    }
    d_sum
}

fn find_cycles(problem: &Problem, ix: usize, graph: Vec<usize>) -> Vec<Vec<usize>> {
    //println!("find cycle: {}, {:?}", ix, graph);
    // Only find smallish cycles
    if graph.len() > 6 {
        return vec![];
    }
    let expand = graph.last().unwrap();
    if graph.len() > 1 && *expand == ix {
        return vec![graph];
    }
    let connected: Vec<_> = problem
        .figure
        .edges
        .iter()
        .filter_map(|(a, b)| {
            if a == expand {
                Some(b)
            } else if b == expand {
                Some(a)
            } else {
                None
            }
        })
        .collect();
    let mut cycles = vec![];
    for i in connected {
        if !graph[1..].contains(&i) {
            let mut g = graph.clone();
            g.push(*i);
            cycles.extend(find_cycles(problem, ix, g));
        }
    }
    cycles
}

struct Solution {
    queue: BinaryHeap<Thing>,
    iterations: usize,
    solved: bool,
    seen: HashSet<Vec<[i32; 2]>>,
    lowest_dislikes: Option<i32>,
    candidates: Vec<((i32, i32), usize)>,
    longest_edge_in_cycle: HashMap<usize, i32>,
}

/// We derive Deserialize/Serialize so we can persist app state on shutdown.
#[cfg_attr(feature = "persistence", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "persistence", serde(default))] // if we add new fields, give them default values when deserializing old state
pub struct PolygonApp {
    filename: String,
    grid: bool,
    try_keep: bool,

    #[cfg_attr(feature = "persistence", serde(skip))]
    problem: Option<Problem>,

    #[cfg_attr(feature = "persistence", serde(skip))]
    solution: Solution,

    #[cfg_attr(feature = "persistence", serde(skip))]
    pose: HashMap<usize, [i32; 2]>,

    #[cfg_attr(feature = "persistence", serde(skip))]
    best_pose: HashMap<usize, [i32; 2]>,

    #[cfg_attr(feature = "persistence", serde(skip))]
    selected: HashSet<usize>,

    #[cfg_attr(feature = "persistence", serde(skip))]
    solving: bool,
}

impl PolygonApp {
    fn init_solution(&mut self) {
        let problem = self.problem.as_ref().unwrap();
        let min_x = problem.hole.iter().map(|x| x[0] as i32).min().unwrap_or(0);
        let max_x = problem.hole.iter().map(|x| x[0] as i32).max().unwrap_or(0);
        let min_y = problem.hole.iter().map(|x| x[1] as i32).min().unwrap_or(0);
        let max_y = problem.hole.iter().map(|x| x[1] as i32).max().unwrap_or(0);
        let x_range: Vec<_> = (min_x..=max_x).collect();
        let y_range: Vec<_> = (min_y..=max_y).collect();
        let mut positions: Vec<_> = x_range
            .into_iter()
            .cartesian_product(y_range.into_iter())
            .filter_map(|(x, y)| {
                // Filter out points outside the hole
                let point = Point {
                    x: x as f32,
                    y: y as f32,
                };
                // Is the point inside?
                if inside(&problem.hole, point) {
                    Some((x, y))
                } else {
                    None
                }
            })
            .collect();
        // Sort the positions on dislikes
        positions.sort_by(|a, b| {
            let mut posea = HashMap::new();
            let mut poseb = HashMap::new();
            posea.insert(0, [a.0, a.1]);
            poseb.insert(0, [b.0, b.1]);
            let da = calc_dislikes(&problem, &posea);
            let db = calc_dislikes(&problem, &poseb);
            da.cmp(&db)
        });
        let candidates: Vec<_> = positions
            .into_iter()
            .cartesian_product(0..problem.figure.vertices.len())
            .collect();
        self.solution = Solution {
            queue: BinaryHeap::new(),
            iterations: 0,
            solved: false,
            seen: HashSet::new(),
            lowest_dislikes: None,
            candidates,
            longest_edge_in_cycle: HashMap::new(),
        };
        // First try with keeping the vertices with no invalid edges
        if self.try_keep {
            let mut valid = HashMap::new();
            for i in 0..problem.figure.vertices.len() {
                let vertex = problem.figure.vertices[i];
                valid.insert(i, [vertex.x as i32, vertex.y as i32]);
            }
            for ix in 0..problem.figure.edges.len() {
                let edge = problem.figure.edges[ix];
                let p1 = problem.figure.vertices[edge.0];
                let p2 = problem.figure.vertices[edge.1];
                let mut ok = inside(&problem.hole, p1) && inside(&problem.hole, p2);
                if ok {
                    let mut ix = 0;
                    let len = problem.hole.len();
                    while ix < len {
                        let jx = (ix + 1) % len;
                        if intersects((p1, p2), (problem.hole[ix], problem.hole[jx])) {
                            ok = false;
                            break;
                        }
                        ix += 1;
                    }
                }
                if !ok {
                    valid.remove(&edge.0);
                    valid.remove(&edge.1);
                }
            }
            let dislikes = calc_dislikes(&problem, &valid);
            let constrainedness = calc_constrainedness(&problem, &valid);
            let edge_length = calc_edge_length(&problem, &valid);
            self.solution.queue.push(Thing {
                num: problem.figure.vertices.len(),
                verts: valid,
                dislikes,
                constrainedness,
                edge_length,
            });
        }

        // Also add one where nothing is placed yet
        let valid = HashMap::new();
        let constrainedness = calc_constrainedness(&problem, &valid);
        let edge_length = calc_edge_length(&problem, &valid);
        self.solution.queue.push(Thing {
            num: problem.figure.vertices.len(),
            verts: valid,
            dislikes: 0,
            constrainedness,
            edge_length,
        });
        println!("find cycles for {} vertices", problem.figure.vertices.len());
        for i in 0..problem.figure.vertices.len() {
            println!("find cycles for vertex {}", i);
            let g = vec![i];
            let cycles = find_cycles(problem, i, g);
            println!("found {} cycles", cycles.len());
            let mut smallest = None;
            for cycle in &cycles {
                if cycle.len() > 3 {
                    if let Some(s) = smallest {
                        if cycle.len() < s {
                            smallest = Some(cycle.len());
                        }
                    } else {
                        smallest = Some(cycle.len());
                    }
                }
            }
            let mut longest_edge = 0;
            let mut unique: HashSet<Vec<usize>> = HashSet::new();
            for cycle in &cycles {
                let cl = cycle.len();
                if Some(cl) == smallest {
                    let c = cycle.to_vec();
                    let mut rev = cycle.to_vec();
                    rev.reverse();
                    if !unique.contains(&c) && !unique.contains(&rev) {
                        unique.insert(c);
                        for i in 0..(cycle.len() - 1) {
                            let pa = problem.figure.vertices[i];
                            let pb =
                                problem.figure.vertices[(i + 1) % problem.figure.vertices.len()];
                            let edge_length = ((pa.x - pb.x) * (pa.x - pb.x)
                                + (pa.y - pb.y) * (pa.y - pb.y))
                                as i32;
                            if edge_length > longest_edge {
                                longest_edge = edge_length;
                            }
                        }
                    }
                }
            }
            self.solution.longest_edge_in_cycle.insert(i, longest_edge);
            for cycle in unique {
                println!("vertex {} has smallest cycle {:?}", i, cycle);
            }
            println!("vertex {} has longest edge {}", i, longest_edge);
        }
    }

    fn solve(&mut self, iterations: usize) {
        // println!("solving: {} iterations", iterations);
        let mut iter = 0;
        let problem = self.problem.as_ref().unwrap();
        while let Some(t) = self.solution.queue.pop() {
            let Thing {
                num: _wanted_num,
                verts,
                dislikes,
                constrainedness: _,
                edge_length: _,
            } = t.clone();
            // println!("{}, {:?}", dislikes, verts);
            let num = verts.len();
            if num == problem.figure.vertices.len() {
                if dislikes == 0 {
                    self.solving = false;
                    self.solution.solved = true;
                    self.best_pose = verts.clone();
                    self.pose = verts;
                    println!("solved! {:?}", self.best_pose);
                    return;
                }
                if let Some(lowest) = self.solution.lowest_dislikes {
                    if dislikes < lowest {
                        println!("New best! {:?}, {}", verts, dislikes);
                        self.solution.lowest_dislikes = Some(dislikes);
                        self.best_pose = verts;
                        self.save_best();
                        return;
                    }
                } else {
                    println!("First solution! {:?}, {}", verts, dislikes);
                    self.solution.lowest_dislikes = Some(dislikes);
                    self.best_pose = verts;
                    self.save_best();
                    return;
                }
            }
            let max_el = self
                .solution
                .longest_edge_in_cycle
                .iter()
                .filter_map(|(ix, el)| {
                    if verts.contains_key(ix) {
                        None
                    } else {
                        Some(el)
                    }
                })
                .max()
                .unwrap_or(&0);
            //	    println!("max el: {}", max_el);
            let cand: Vec<_> = self
                .solution
                .candidates
                .iter()
                .filter_map(|((x, y), ix)| {
                    if let Some(le) = self.solution.longest_edge_in_cycle.get(&ix) {
                        //			println!("le for {} is {}", ix, le);
                        if le == max_el {
                            Some(((x, y), ix))
                        } else {
                            None
                        }
                    } else {
                        Some(((x, y), ix))
                    }
                })
                .collect();
            //	    println!("cand: {:?}", cand);
            // If we didn't use all the nighbors, put it back in
            if cand.len() < num {
                self.solution.queue.push(t);
            }
            let ok_points: Vec<_> = cand
                .into_par_iter()
                .filter_map(|((x, y), ix)| {
                    let point = Point {
                        x: *x as f32,
                        y: *y as f32,
                    };
                    // Are all placed edges inside?
                    // Are all placed edges within the constraint?
                    let mut ok = true;
                    for (a, b) in &problem.figure.edges {
                        // Only check the new edges enabled by this point
                        if (*a == *ix || verts.contains_key(a))
                            && (*b == *ix || verts.contains_key(b))
                        {
                            let p1 = problem.figure.vertices[*a];
                            let p2 = problem.figure.vertices[*b];
                            let distance = (p1[0] - p2[0]) * (p1[0] - p2[0])
                                + (p1[1] - p2[1]) * (p1[1] - p2[1]);
                            let pp1 = if *a == *ix {
                                point
                            } else {
                                let [x, y] = verts[a];
                                Point {
                                    x: x as f32,
                                    y: y as f32,
                                }
                            };
                            let pp2 = if *b == *ix {
                                point
                            } else {
                                let [x, y] = verts[b];
                                Point {
                                    x: x as f32,
                                    y: y as f32,
                                }
                            };
                            let new_distance = (pp1[0] - pp2[0]) * (pp1[0] - pp2[0])
                                + (pp1[1] - pp2[1]) * (pp1[1] - pp2[1]);
                            let eps =
                                1000000.0 * ((new_distance as f32 / distance as f32) - 1.0).abs();
                            if eps as i32 > problem.epsilon {
                                ok = false;
                                // println!("edge too long/short: {}, {}, {}, {}", distance, new_distance, eps, problem.epsilon);
                                break;
                            }
                            // check intersections
                            let mut ix = 0;
                            let len = problem.hole.len();
                            while ix < len {
                                let jx = (ix + 1) % len;
                                if intersects((pp1, pp2), (problem.hole[ix], problem.hole[jx])) {
                                    // println!("edge intersects with hole");
                                    ok = false;
                                    break;
                                }
                                ix += 1;
                            }
                        }
                    }
                    if ok {
                        Some((ix, point))
                    } else {
                        None
                    }
                })
                .collect();
            for (ix, point) in ok_points {
                let mut points = vec![];
                for i in 0..problem.figure.vertices.len() {
                    if i == *ix {
                        points.push([point.x as i32, point.y as i32]);
                    } else {
                        let vertex = problem.figure.vertices[i];
                        let pos = [vertex.x as i32, vertex.y as i32];
                        let p = verts.get(&i).unwrap_or(&pos);
                        points.push(*p);
                    }
                }
                if self.solution.seen.insert(points) {
                    let mut v = verts.clone();
                    v.insert(*ix, [point.x as i32, point.y as i32]);
                    let new_dislikes = calc_dislikes(&problem, &v);
                    let new_constrainedness = calc_constrainedness(&problem, &v);
                    let new_edge_length = calc_edge_length(&problem, &v);
                    // println!("insert new thing: {:?}", v);
                    let t = Thing {
                        num: problem.figure.vertices.len(),
                        verts: v,
                        dislikes: new_dislikes,
                        constrainedness: new_constrainedness,
                        edge_length: new_edge_length,
                    };
                    self.solution.queue.push(t);
                }
            }
            iter += 1;
            self.solution.iterations += 1;
            if iter > iterations {
                self.pose = verts;
                // println!("solved {} iterations, {:?}", iterations, self.pose);
                return;
            }
        }
        println!("no more in queue");
        self.solving = false;
    }

    fn save_best(&self) {
        let out_filename = self.filename.to_owned() + ".solution.json";
        let mut pose = vec![];
        let problem = self.problem.as_ref().unwrap();
        for i in 0..problem.figure.vertices.len() {
            if let Some(p) = self.best_pose.get(&i) {
                pose.push(Point {
                    x: p[0] as f32,
                    y: p[1] as f32,
                });
            } else {
                pose.push(problem.figure.vertices[i])
            }
        }
        if write_solution_to_file(&out_filename, &pose).is_ok() {
            println!("saved solution!");
        }
    }

    fn load(&mut self) {
        if let Ok(read_problem) = read_problem_from_file(&self.filename) {
            let mut pose = HashMap::new();
            for i in 0..read_problem.figure.vertices.len() {
                let p = read_problem.figure.vertices[i];
                pose.insert(i, [p.x as i32, p.y as i32]);
            }
            self.pose = pose;
            self.best_pose = HashMap::new();
            self.problem = Some(read_problem);
            self.init_solution();
            println!("Problem: {:?}", self.problem);
        }
    }

    pub fn solve_console(&mut self, filename: &str) {
        self.filename = filename.to_owned();
        self.load();
        self.solving = true;
        while self.solving {
            self.solve(10);
            print!(".");
            std::io::stdout().flush().unwrap();
        }
    }
}

impl Default for PolygonApp {
    fn default() -> Self {
        Self {
            // Example stuff:
            filename: "problems/1.problem".to_owned(),
            grid: false,
            try_keep: false,
            problem: None,
            solution: Solution {
                queue: BinaryHeap::new(),
                iterations: 0,
                solved: false,
                seen: HashSet::new(),
                lowest_dislikes: None,
                candidates: vec![],
                longest_edge_in_cycle: HashMap::new(),
            },
            pose: HashMap::new(),
            selected: HashSet::new(),
            solving: false,
            best_pose: HashMap::new(),
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
                ui.text_edit_singleline(&mut self.filename);
            });

            if ui.button("Load").clicked() {
                self.load();
            }
            ui.checkbox(&mut self.grid, "Show Grid");
            ui.checkbox(&mut self.try_keep, "Try to keep valid vertices");
            ui.group(|ui| {
                ui.set_enabled(self.problem.is_some());
                ui.checkbox(&mut self.solving, "Solve");
                ui.label(format!("Iterations: {}", self.solution.iterations));
                ui.label(format!(
                    "Dislikes: {}",
                    self.solution.lowest_dislikes.unwrap_or(-1)
                ));
                ui.label(format!("Solved: {}", self.solution.solved));
                if ui.button("Save").clicked() {
                    self.save_best();
                }
            });
            let mut selected = self.selected.clone();
            if let Some(problem) = &self.problem {
                if ui.button("Select All").clicked() {
                    selected = HashSet::new();
                }
                egui::ScrollArea::auto_sized().show(ui, |ui| {
                    for i in 0..problem.figure.edges.len() {
                        let checked = selected.contains(&i);
                        let (a, b) = problem.figure.edges[i];
                        let p1 = problem.figure.vertices[a];
                        let p2 = problem.figure.vertices[b];
                        let d =
                            (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]);
                        let edge = problem.figure.edges[i];
                        let op1 = problem.figure.vertices[edge.0];
                        let opp1 = [op1.x as i32, op1.y as i32];
                        let op2 = problem.figure.vertices[edge.1];
                        let opp2 = [op2.x as i32, op2.y as i32];
                        let p1 = self.best_pose.get(&edge.0).unwrap_or(&opp1);
                        let p2 = self.best_pose.get(&edge.1).unwrap_or(&opp2);
                        let edge_valid =
                            self.pose.contains_key(&edge.0) && self.pose.contains_key(&edge.1);
                        let pp1 = Point {
                            x: p1[0] as f32,
                            y: p1[1] as f32,
                        };
                        let pp2 = Point {
                            x: p2[0] as f32,
                            y: p2[1] as f32,
                        };
                        let dd = (pp1[0] - pp2[0]) * (pp1[0] - pp2[0])
                            + (pp1[1] - pp2[1]) * (pp1[1] - pp2[1]);
                        let eps = 1000000.0 * ((dd as f32 / d as f32) - 1.0).abs();
                        let label = format!(
                            "{:?}, {}, {}, {}, {}, {}",
                            problem.figure.edges[i], edge_valid, d, dd, eps, problem.epsilon
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
            self.selected = selected;

            if self.solving {
                self.solve(10);
                ui.ctx().request_repaint();
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
                let invalid_stroke = egui::Stroke::new(1.0, egui::Color32::GOLD);
                let best_stroke = egui::Stroke::new(2.0, egui::Color32::BLUE);
                let grid_stroke = egui::Stroke::new(1.0, egui::Color32::LIGHT_GRAY);
                if let Some(problem) = &self.problem {
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
                    if self.grid {
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
                        if !self.selected.is_empty() && !self.selected.contains(&ix) {
                            continue;
                        }
                        let edge = problem.figure.edges[ix];
                        let op1 = problem.figure.vertices[edge.0];
                        let opp1 = [op1.x as i32, op1.y as i32];
                        let op2 = problem.figure.vertices[edge.1];
                        let opp2 = [op2.x as i32, op2.y as i32];
                        let p1 = self.pose.get(&edge.0).unwrap_or(&opp1);
                        let p2 = self.pose.get(&edge.1).unwrap_or(&opp2);
                        let pp1 = Point {
                            x: p1[0] as f32,
                            y: p1[1] as f32,
                        };
                        let pp2 = Point {
                            x: p2[0] as f32,
                            y: p2[1] as f32,
                        };
                        let points = vec![to_screen * pp1, to_screen * pp2];
                        let mut stroke =
                            if !inside(&problem.hole, pp1) || !inside(&problem.hole, pp2) {
                                outside_stroke
                            } else {
                                pose_stroke
                            };
                        let mut i = 0;
                        let l = problem.hole.len();
                        while i < l {
                            let j = (i + 1) % l;
                            if intersects((pp1, pp2), (problem.hole[i], problem.hole[j])) {
                                stroke = intersect_stroke;
                                break;
                            }
                            i += 1;
                        }
                        let edge_valid =
                            self.pose.contains_key(&edge.0) && self.pose.contains_key(&edge.1);
                        if !edge_valid {
                            stroke = invalid_stroke;
                        }
                        shapes.push(egui::Shape::line(points, stroke));
                    }
                    // Draw best pose
                    for ix in 0..problem.figure.edges.len() {
                        if !self.selected.is_empty() && !self.selected.contains(&ix) {
                            continue;
                        }
                        let edge = problem.figure.edges[ix];
                        if !self.best_pose.contains_key(&edge.0)
                            || !self.best_pose.contains_key(&edge.1)
                        {
                            continue;
                        }
                        let p1 = self.best_pose.get(&edge.0).unwrap();
                        let p2 = self.best_pose.get(&edge.1).unwrap();
                        let pp1 = Point {
                            x: p1[0] as f32,
                            y: p1[1] as f32,
                        };
                        let pp2 = Point {
                            x: p2[0] as f32,
                            y: p2[1] as f32,
                        };
                        let points = vec![to_screen * pp1, to_screen * pp2];
                        shapes.push(egui::Shape::line(points, best_stroke));
                    }
                }
                ui.painter().extend(shapes);
            });
        });
    }
}
