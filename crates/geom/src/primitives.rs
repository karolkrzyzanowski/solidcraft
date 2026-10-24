//! Finding planes and cylinders in a scanned mesh: random sample consensus (Fischler and Bolles,
//! 1981) on the mesh's vertices and their normals, in the manner Schnabel, Wahl and Klein
//! describe for primitives (2007). A shape is guessed from one vertex (a plane) or two nearby
//! ones (a cylinder, from where their normals meet), the vertices that lie on it and face the way
//! it does are counted, and the best guess is refitted to its largest connected patch of the mesh.
//! Those vertices are taken out and the search goes on. Everything is arithmetic on finite
//! numbers; input that is not finite or not a mesh gives no shapes.

use std::collections::HashMap;

use crate::Vec3;
use crate::registration::jacobi;

/// A mesh as vertices with unit normals and the triangles joining them.
#[derive(Clone, Debug, Default)]
pub struct Cloud {
    pub points: Vec<Vec3>,
    pub normals: Vec<Vec3>,
    pub triangles: Vec<[u32; 3]>,
}

impl Cloud {
    /// A cloud from a welded mesh, each vertex's normal the mean of the faces around it weighted
    /// by their area. `None` for a vertex that is not finite or a triangle that names a vertex
    /// that is not there.
    pub fn from_mesh(points: Vec<Vec3>, triangles: Vec<[u32; 3]>) -> Option<Cloud> {
        if points.iter().any(|p| !p.is_finite()) || triangles.iter().flatten().any(|i| *i as usize >= points.len()) {
            return None;
        }
        let mut sum = vec![Vec3::ZERO; points.len()];
        for t in &triangles {
            let (a, b, c) = (*points.get(t[0] as usize)?, *points.get(t[1] as usize)?, *points.get(t[2] as usize)?);
            // The cross product's length is twice the area: the weight comes with it.
            let n = (b - a).cross(c - a);
            for i in t {
                if let Some(s) = sum.get_mut(*i as usize) {
                    *s += n;
                }
            }
        }
        let normals = sum.into_iter().map(|n| n.normalized().unwrap_or(Vec3::Z)).collect();
        Some(Cloud { points, normals, triangles })
    }
}

/// What was found.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Shape {
    /// The plane through `origin` with unit `normal` (turned the way most of its vertices face).
    Plane { origin: Vec3, normal: Vec3 },
    /// The cylinder about the line through `origin` along unit `axis`, from `min` to `max` along it
    /// (measured from `origin`).
    Cylinder { origin: Vec3, axis: Vec3, radius: f64, min: f64, max: f64 },
}

/// A shape and the vertices on it.
#[derive(Clone, Debug, PartialEq)]
pub struct Found {
    pub shape: Shape,
    /// Indices of the vertices on the shape (one connected patch of the mesh).
    pub inliers: Vec<u32>,
    /// RMS distance of those vertices from the shape (mm).
    pub rms: f64,
}

/// How [`detect`] searches.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub planes: bool,
    pub cylinders: bool,
    /// How far a vertex may lie from a shape (mm).
    pub tolerance: f64,
    /// How far a vertex's normal may turn from the shape's (radians).
    pub angle: f64,
    /// Fewest vertices a shape must have.
    pub min_points: usize,
    /// Most shapes to find.
    pub max_shapes: usize,
    /// Guesses per shape.
    pub trials: usize,
    /// Seed of the generator: the same seed finds the same shapes.
    pub seed: u64,
}

impl Default for Options {
    fn default() -> Self {
        Options { planes: true, cylinders: true, tolerance: 0.2, angle: 15f64.to_radians(), min_points: 200, max_shapes: 20, trials: 400, seed: 1 }
    }
}

/// Most vertices a guess is scored against (an even spread of the rest).
const SCORE_CAP: usize = 40_000;
/// A cylinder's radius may not be less than this times the tolerance, nor more than this times the
/// size of the mesh (it is a plane then).
const MIN_RADIUS_TOLERANCES: f64 = 3.0;
const MAX_RADIUS_SIZES: f64 = 5.0;
/// A cylinder takes a patch from a plane only with this many times its vertices.
const CYLINDER_BIAS: f64 = 1.15;
/// A cylinder whose arc bulges less than this many tolerances is a plane.
const FLAT_SAGITTAS: f64 = 2.0;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn below(&mut self, n: usize) -> usize {
        if n == 0 { 0 } else { (self.next() % n as u64) as usize }
    }
}

/// For each vertex, those it shares a triangle side with.
struct Adjacency {
    start: Vec<u32>,
    list: Vec<u32>,
}

impl Adjacency {
    fn new(n: usize, triangles: &[[u32; 3]]) -> Adjacency {
        let mut degree = vec![0u32; n + 1];
        for t in triangles {
            for i in t {
                if let Some(d) = degree.get_mut(*i as usize) {
                    *d += 2;
                }
            }
        }
        let mut start = vec![0u32; n + 1];
        for i in 0..n {
            start[i + 1] = start[i] + degree[i];
        }
        let mut fill = start.clone();
        let mut list = vec![0u32; start[n] as usize];
        for t in triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                for (from, to) in [(a, b), (b, a)] {
                    if let Some(f) = fill.get_mut(from as usize)
                        && let Some(slot) = list.get_mut(*f as usize)
                    {
                        *slot = to;
                        *f += 1;
                    }
                }
            }
        }
        Adjacency { start, list }
    }

    fn of(&self, i: u32) -> &[u32] {
        let (a, b) = (self.start.get(i as usize).copied().unwrap_or(0) as usize, self.start.get(i as usize + 1).copied().unwrap_or(0) as usize);
        self.list.get(a..b).unwrap_or(&[])
    }
}

/// How far `p` (facing `n`) is from the shape, `None` when it faces another way.
fn residual(shape: &Shape, p: Vec3, n: Vec3, cos: f64) -> Option<f64> {
    match shape {
        Shape::Plane { origin, normal } => {
            if n.dot(*normal).abs() < cos {
                return None;
            }
            Some((p - *origin).dot(*normal).abs())
        }
        Shape::Cylinder { origin, axis, radius, .. } => {
            let v = p - *origin;
            let radial = v - *axis * v.dot(*axis);
            let len = radial.len();
            if len < 1e-12 || n.dot(radial / len).abs() < cos {
                return None;
            }
            Some((len - radius).abs())
        }
    }
}

fn plane_through(p: Vec3, n: Vec3) -> Shape {
    Shape::Plane { origin: p, normal: n }
}

/// The cylinder two vertices with their normals suggest: its axis is square to both normals, and
/// in the plane square to the axis the normals' lines meet at its centre.
fn cylinder_through(pa: Vec3, na: Vec3, pb: Vec3, nb: Vec3) -> Option<Shape> {
    let axis = na.cross(nb);
    let len = axis.len();
    // Normals that are nearly parallel give no axis.
    if len < 8f64.to_radians().sin() {
        return None;
    }
    let axis = axis / len;
    let u = axis.any_perp();
    let v = axis.cross(u);
    let flat = |p: Vec3| (p.dot(u), p.dot(v));
    let (a, b) = (flat(pa), flat(pb));
    let (da, db) = (flat(na), flat(nb));
    let (dal, dbl) = ((da.0 * da.0 + da.1 * da.1).sqrt(), (db.0 * db.0 + db.1 * db.1).sqrt());
    if dal < 1e-9 || dbl < 1e-9 {
        return None;
    }
    let (da, db) = ((da.0 / dal, da.1 / dal), (db.0 / dbl, db.1 / dbl));
    // a + s·da = b + t·db
    let det = db.0 * da.1 - da.0 * db.1;
    if det.abs() < 1e-9 {
        return None;
    }
    let (rx, ry) = (b.0 - a.0, b.1 - a.1);
    let s = (-rx * db.1 + db.0 * ry) / det;
    let centre = (a.0 + s * da.0, a.1 + s * da.1);
    let radius = s.abs();
    if !(radius.is_finite() && radius > 0.0) {
        return None;
    }
    Some(Shape::Cylinder { origin: u * centre.0 + v * centre.1 + axis * pa.dot(axis), axis, radius, min: 0.0, max: 0.0 })
}

fn refit_plane(cloud: &Cloud, inliers: &[u32]) -> Option<Shape> {
    let pts: Vec<Vec3> = inliers.iter().filter_map(|i| cloud.points.get(*i as usize).copied()).collect();
    let c = crate::registration::centroid(&pts)?;
    let mut m = [[0.0; 3]; 3];
    for p in &pts {
        let d = [p.x - c.x, p.y - c.y, p.z - c.z];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell += d[i] * d[j];
            }
        }
    }
    let (values, vectors) = jacobi(m);
    let k = (0..3).fold(0, |best, i| if values[i] < values[best] { i } else { best });
    let mut normal = Vec3::new(vectors[0][k], vectors[1][k], vectors[2][k]).normalized()?;
    let facing: Vec3 = inliers.iter().filter_map(|i| cloud.normals.get(*i as usize).copied()).fold(Vec3::ZERO, |s, n| s + n);
    if normal.dot(facing) < 0.0 {
        normal = -normal;
    }
    Some(Shape::Plane { origin: c, normal })
}

fn refit_cylinder(cloud: &Cloud, inliers: &[u32]) -> Option<Shape> {
    let pts: Vec<Vec3> = inliers.iter().filter_map(|i| cloud.points.get(*i as usize).copied()).collect();
    let c = crate::registration::centroid(&pts)?;
    // The axis is the direction the normals never point along.
    let mut m = [[0.0; 3]; 3];
    for i in inliers {
        let n = cloud.normals.get(*i as usize).copied()?;
        let d = [n.x, n.y, n.z];
        for (a, row) in m.iter_mut().enumerate() {
            for (b, cell) in row.iter_mut().enumerate() {
                *cell += d[a] * d[b];
            }
        }
    }
    let (values, vectors) = jacobi(m);
    let mut order = [0usize, 1, 2];
    order.sort_by(|a, b| values[*a].total_cmp(&values[*b]));
    // Normals in one direction only (a plane) leave the axis open.
    if values[order[1]] < 0.05 * inliers.len() as f64 {
        return None;
    }
    let k = order[0];
    let axis = Vec3::new(vectors[0][k], vectors[1][k], vectors[2][k]).normalized()?;
    let u = axis.any_perp();
    let v = axis.cross(u);
    // Circle through the vertices seen along the axis (Kåsa): x² + y² + D·x + E·y + F = 0.
    let (mut sxx, mut sxy, mut syy, mut sx, mut sy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut sxz, mut syz, mut sz) = (0.0, 0.0, 0.0);
    for p in &pts {
        let d = *p - c;
        let (x, y) = (d.dot(u), d.dot(v));
        let z = x * x + y * y;
        sxx += x * x;
        sxy += x * y;
        syy += y * y;
        sx += x;
        sy += y;
        sxz += x * z;
        syz += y * z;
        sz += z;
    }
    let n = pts.len() as f64;
    let (a, b, cc) = ([sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n]);
    let rhs = [-sxz, -syz, -sz];
    let det3 = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let det = det3([a, b, cc]);
    if det.abs() < 1e-12 * (sxx * syy * n).abs().max(1e-300) {
        return None;
    }
    let col = |k: usize| {
        let mut rows = [a, b, cc];
        for (r, row) in rows.iter_mut().enumerate() {
            row[k] = rhs[r];
        }
        det3(rows) / det
    };
    let (d, e, f) = (col(0), col(1), col(2));
    let (cx, cy) = (-d / 2.0, -e / 2.0);
    let r2 = cx * cx + cy * cy - f;
    if !(r2.is_finite() && r2 > 0.0) {
        return None;
    }
    let origin = c + u * cx + v * cy;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in &pts {
        let t = (*p - origin).dot(axis);
        lo = lo.min(t);
        hi = hi.max(t);
    }
    Some(Shape::Cylinder { origin, axis, radius: r2.sqrt(), min: lo, max: hi })
}

/// Scratch space for finding connected patches.
struct Work {
    flag: Vec<u32>,
    seen: Vec<u32>,
    stamp: u32,
}

impl Work {
    fn new(n: usize) -> Work {
        Work { flag: vec![0; n], seen: vec![0; n], stamp: 0 }
    }

    /// The largest set of `inliers` joined by triangle sides.
    fn largest_patch(&mut self, inliers: &[u32], adj: &Adjacency) -> Vec<u32> {
        self.stamp += 1;
        let stamp = self.stamp;
        for i in inliers {
            if let Some(f) = self.flag.get_mut(*i as usize) {
                *f = stamp;
            }
        }
        let mut best: Vec<u32> = Vec::new();
        for &start in inliers {
            if self.seen.get(start as usize).copied() == Some(stamp) {
                continue;
            }
            let mut patch = vec![start];
            if let Some(s) = self.seen.get_mut(start as usize) {
                *s = stamp;
            }
            let mut k = 0;
            while k < patch.len() {
                let v = patch[k];
                k += 1;
                for &w in adj.of(v) {
                    if self.flag.get(w as usize).copied() == Some(stamp)
                        && let Some(s) = self.seen.get_mut(w as usize)
                        && *s != stamp
                    {
                        *s = stamp;
                        patch.push(w);
                    }
                }
            }
            if patch.len() > best.len() {
                best = patch;
            }
        }
        best
    }
}

fn on_shape(cloud: &Cloud, shape: &Shape, alive: &[u32], tol: f64, cos: f64) -> Vec<u32> {
    alive
        .iter()
        .copied()
        .filter(|i| match (cloud.points.get(*i as usize), cloud.normals.get(*i as usize)) {
            (Some(p), Some(n)) => residual(shape, *p, *n, cos).is_some_and(|d| d <= tol),
            _ => false,
        })
        .collect()
}

/// A guess made exact: the vertices on it, the largest patch of them, the shape fitted to the patch,
/// again a few times.
fn refine(cloud: &Cloud, adj: &Adjacency, work: &mut Work, guess: Shape, alive: &[u32], opts: &Options) -> Option<Found> {
    let cos = opts.angle.cos();
    let mut shape = guess;
    for _ in 0..4 {
        let on = on_shape(cloud, &shape, alive, opts.tolerance, cos);
        let patch = work.largest_patch(&on, adj);
        if patch.len() < 3 {
            return None;
        }
        shape = match shape {
            Shape::Plane { .. } => refit_plane(cloud, &patch)?,
            Shape::Cylinder { .. } => refit_cylinder(cloud, &patch)?,
        };
    }
    // What the last fit holds, and how well.
    let on = on_shape(cloud, &shape, alive, opts.tolerance, cos);
    let patch = work.largest_patch(&on, adj);
    if patch.len() < 3 {
        return None;
    }
    shape = match shape {
        Shape::Plane { .. } => refit_plane(cloud, &patch)?,
        Shape::Cylinder { .. } => refit_cylinder(cloud, &patch)?,
    };
    let sum: f64 =
        patch.iter().filter_map(|i| Some(residual(&shape, *cloud.points.get(*i as usize)?, *cloud.normals.get(*i as usize)?, -1.0)?.powi(2))).sum();
    let rms = (sum / patch.len() as f64).sqrt();
    rms.is_finite().then_some(Found { shape, inliers: patch, rms })
}

/// How far the middle of a cylinder's arc stands out from the chord across its ends (mm): small
/// means the patch is flat to within the scan's noise and a plane describes it as well.
fn sagitta(cloud: &Cloud, f: &Found) -> f64 {
    let Shape::Cylinder { origin, axis, radius, .. } = f.shape else { return f64::INFINITY };
    let u = axis.any_perp();
    let v = axis.cross(u);
    let angles: Vec<f64> = f
        .inliers
        .iter()
        .filter_map(|i| {
            let d = *cloud.points.get(*i as usize)? - origin;
            Some(d.dot(v).atan2(d.dot(u)))
        })
        .collect();
    if angles.is_empty() {
        return f64::INFINITY;
    }
    // The mean direction, then the widest turn away from it either way.
    let (s, c) = angles.iter().fold((0.0, 0.0), |(s, c), a| (s + a.sin(), c + a.cos()));
    let mean = s.atan2(c);
    let wrap = |a: f64| {
        let d = (a - mean).rem_euclid(std::f64::consts::TAU);
        if d > std::f64::consts::PI { d - std::f64::consts::TAU } else { d }
    };
    let (lo, hi) = angles.iter().map(|a| wrap(*a)).fold((0.0f64, 0.0f64), |(lo, hi), d| (lo.min(d), hi.max(d)));
    let half = ((hi - lo) / 2.0).min(std::f64::consts::PI);
    radius * (1.0 - half.cos())
}

fn rms_of(cloud: &Cloud, shape: &Shape, inliers: &[u32]) -> f64 {
    let sum: f64 =
        inliers.iter().filter_map(|i| Some(residual(shape, *cloud.points.get(*i as usize)?, *cloud.normals.get(*i as usize)?, -1.0)?.powi(2))).sum();
    (sum / inliers.len().max(1) as f64).sqrt()
}

/// Whether two shapes are one: planes with the same facing and place, cylinders with the same axis
/// line and radius.
fn same_shape(a: &Shape, b: &Shape, tol: f64) -> bool {
    match (a, b) {
        (Shape::Plane { origin: oa, normal: na }, Shape::Plane { origin: ob, normal: nb }) => {
            na.dot(*nb).abs() > 5f64.to_radians().cos() && (*ob - *oa).dot(*na).abs() < tol && (*oa - *ob).dot(*nb).abs() < tol
        }
        (Shape::Cylinder { origin: oa, axis: xa, radius: ra, .. }, Shape::Cylinder { origin: ob, axis: xb, radius: rb, .. }) => {
            let d = *ob - *oa;
            let off = d - *xa * d.dot(*xa);
            xa.dot(*xb).abs() > 3f64.to_radians().cos() && off.len() < tol && (ra - rb).abs() < tol
        }
        _ => false,
    }
}

/// Found shapes that are the same shape joined into one, refitted to all their vertices.
fn merge(cloud: &Cloud, found: Vec<Found>, tol: f64) -> Vec<Found> {
    let mut out: Vec<Found> = Vec::new();
    for f in found {
        let Some(k) = out.iter().position(|o| same_shape(&o.shape, &f.shape, tol)) else {
            out.push(f);
            continue;
        };
        let mut all = out[k].inliers.clone();
        all.extend_from_slice(&f.inliers);
        let refit = match f.shape {
            Shape::Plane { .. } => refit_plane(cloud, &all),
            Shape::Cylinder { .. } => refit_cylinder(cloud, &all),
        };
        match refit {
            Some(shape) => {
                let rms = rms_of(cloud, &shape, &all);
                // Only if the joined fit is still as good as the two were.
                if rms <= out[k].rms.max(f.rms) * 1.5 + tol * 0.2 {
                    out[k] = Found { shape, inliers: all, rms };
                } else {
                    out.push(f);
                }
            }
            None => out.push(f),
        }
    }
    out.sort_by_key(|f| std::cmp::Reverse(f.inliers.len()));
    out
}

/// Planes and cylinders on the mesh, the one with the most vertices first.
pub fn detect(cloud: &Cloud, opts: &Options) -> Vec<Found> {
    let n = cloud.points.len();
    let finite = |x: f64| x.is_finite() && x > 0.0;
    if n < 3
        || cloud.normals.len() != n
        || !(finite(opts.tolerance) && finite(opts.angle) && opts.angle < std::f64::consts::FRAC_PI_2)
        || opts.min_points < 3
        || opts.trials == 0
        || opts.max_shapes == 0
        || !(opts.planes || opts.cylinders)
        || cloud.triangles.iter().flatten().any(|i| *i as usize >= n)
    {
        return Vec::new();
    }
    let size = crate::registration::bounds(&cloud.points).diagonal();
    if !finite(size) {
        return Vec::new();
    }
    let cos = opts.angle.cos();
    let adj = Adjacency::new(n, &cloud.triangles);
    let mut work = Work::new(n);
    let mut rng = Rng(opts.seed ^ 0x9E37_79B9_7F4A_7C15);
    let mut alive = vec![true; n];
    let mut out: Vec<Found> = Vec::new();
    for _ in 0..opts.max_shapes {
        let live: Vec<u32> = (0..n as u32).filter(|i| alive.get(*i as usize).copied().unwrap_or(false)).collect();
        if live.len() < opts.min_points {
            break;
        }
        let step = live.len().div_ceil(SCORE_CAP);
        let sample: Vec<u32> = live.iter().step_by(step.max(1)).copied().collect();
        let mut best_plane: Option<(usize, Shape)> = None;
        let mut best_cyl: Option<(usize, Shape)> = None;
        for _ in 0..opts.trials {
            let a = live[rng.below(live.len())];
            let (Some(pa), Some(na)) = (cloud.points.get(a as usize).copied(), cloud.normals.get(a as usize).copied()) else { continue };
            if opts.planes {
                let guess = plane_through(pa, na);
                let score = on_shape(cloud, &guess, &sample, opts.tolerance, cos).len();
                if best_plane.as_ref().is_none_or(|(s, _)| score > *s) {
                    best_plane = Some((score, guess));
                }
            }
            if opts.cylinders {
                // A second vertex a few steps along the mesh from the first.
                let mut b = a;
                for _ in 0..(3 + rng.below(12)) {
                    let next = adj.of(b);
                    if next.is_empty() {
                        break;
                    }
                    b = next[rng.below(next.len())];
                }
                if b == a || !alive.get(b as usize).copied().unwrap_or(false) {
                    continue;
                }
                let (Some(pb), Some(nb)) = (cloud.points.get(b as usize).copied(), cloud.normals.get(b as usize).copied()) else { continue };
                let Some(guess) = cylinder_through(pa, na, pb, nb) else { continue };
                if let Shape::Cylinder { radius, .. } = guess
                    && !(MIN_RADIUS_TOLERANCES * opts.tolerance..=MAX_RADIUS_SIZES * size).contains(&radius)
                {
                    continue;
                }
                let score = on_shape(cloud, &guess, &sample, opts.tolerance, cos).len();
                if best_cyl.as_ref().is_none_or(|(s, _)| score > *s) {
                    best_cyl = Some((score, guess));
                }
            }
        }
        let plane = best_plane.and_then(|(_, g)| refine(cloud, &adj, &mut work, g, &live, opts));
        let cylinder = best_cyl
            .and_then(|(_, g)| refine(cloud, &adj, &mut work, g, &live, opts))
            .filter(|f| matches!(f.shape, Shape::Cylinder { radius, .. } if (MIN_RADIUS_TOLERANCES * opts.tolerance..=MAX_RADIUS_SIZES * size).contains(&radius)))
            .filter(|f| sagitta(cloud, f) > FLAT_SAGITTAS * opts.tolerance);
        let chosen = match (plane, cylinder) {
            (Some(p), Some(c)) => {
                if c.inliers.len() as f64 > CYLINDER_BIAS * p.inliers.len() as f64 {
                    c
                } else {
                    p
                }
            }
            (Some(f), None) | (None, Some(f)) => f,
            (None, None) => break,
        };
        if chosen.inliers.len() < opts.min_points {
            break;
        }
        for i in &chosen.inliers {
            if let Some(a) = alive.get_mut(*i as usize) {
                *a = false;
            }
        }
        out.push(chosen);
    }
    merge(cloud, out, opts.tolerance)
}

/// A coarser mesh: the vertices that fall in one cell of a grid `cell` mm across become one at their
/// mean, triangles that collapse go. `None` for a cell that is not a positive finite number, a vertex
/// that is not finite or a triangle that names a vertex that is not there.
pub fn cluster_decimate(points: &[Vec3], triangles: &[[u32; 3]], cell: f64) -> Option<(Vec<Vec3>, Vec<[u32; 3]>)> {
    if !(cell.is_finite() && cell > 0.0) || triangles.iter().flatten().any(|i| *i as usize >= points.len()) {
        return None;
    }
    let mut index: HashMap<[i64; 3], u32> = HashMap::new();
    let mut sums: Vec<(Vec3, f64)> = Vec::new();
    let mut map: Vec<u32> = Vec::with_capacity(points.len());
    for p in points {
        let k = |x: f64| {
            let q = (x / cell).floor();
            (q.is_finite() && q.abs() < 1e15).then_some(q as i64)
        };
        let key = [k(p.x)?, k(p.y)?, k(p.z)?];
        let next = u32::try_from(sums.len()).ok()?;
        let id = *index.entry(key).or_insert(next);
        if id == next {
            sums.push((Vec3::ZERO, 0.0));
        }
        let slot = sums.get_mut(id as usize)?;
        slot.0 += *p;
        slot.1 += 1.0;
        map.push(id);
    }
    let merged: Vec<Vec3> = sums.iter().map(|(s, c)| *s / *c).collect();
    let tris = triangles
        .iter()
        .filter_map(|t| {
            let (a, b, c) = (*map.get(t[0] as usize)?, *map.get(t[1] as usize)?, *map.get(t[2] as usize)?);
            (a != b && b != c && c != a).then_some([a, b, c])
        })
        .collect();
    Some((merged, tris))
}

#[cfg(test)]
#[path = "primitives_tests.rs"]
mod tests;
