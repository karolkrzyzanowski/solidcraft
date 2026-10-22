//! Placing a scanned mesh: its principal axes, rigid best fit of one point set onto another
//! (Horn's closed form), iterative closest point against a triangle mesh, and the distances of
//! points from a reference surface. All of it is plain arithmetic on finite numbers; input that
//! is not finite, or too small to define a placement, gives `None`.

use crate::mesh_distance::TriangleSet;
use crate::{Aabb3, Vec3};

fn arr(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// A rigid motion: `p' = r · p + t`, with `r` a rotation (rows).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rigid {
    pub r: [[f64; 3]; 3],
    pub t: Vec3,
}

impl Rigid {
    pub const IDENTITY: Rigid = Rigid { r: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], t: Vec3::new(0.0, 0.0, 0.0) };

    pub fn apply(&self, p: Vec3) -> Vec3 {
        let [a, b, c] = self.r;
        Vec3::new(a[0] * p.x + a[1] * p.y + a[2] * p.z, b[0] * p.x + b[1] * p.y + b[2] * p.z, c[0] * p.x + c[1] * p.y + c[2] * p.z) + self.t
    }

    /// `self` after `first`: `p ↦ self(first(p))`.
    pub fn after(&self, first: &Rigid) -> Rigid {
        let mut r = [[0.0; 3]; 3];
        for (i, row) in r.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = (0..3).map(|k| self.r[i][k] * first.r[k][j]).sum();
            }
        }
        Rigid { r, t: self.apply(first.t) }
    }

    /// The rotation as a unit axis and an angle in radians (0 to π).
    pub fn axis_angle(&self) -> (Vec3, f64) {
        let r = &self.r;
        let v = Vec3::new(r[2][1] - r[1][2], r[0][2] - r[2][0], r[1][0] - r[0][1]) * 0.5;
        let cos = ((r[0][0] + r[1][1] + r[2][2]) - 1.0) * 0.5;
        let sin = v.len();
        let angle = sin.atan2(cos);
        if sin > 1e-9 {
            return (v / sin, angle);
        }
        if cos > 0.0 {
            return (Vec3::Z, 0.0);
        }
        // Half a turn: the axis is the column of r + I with the largest diagonal.
        let k = (0..3).fold(0, |best, i| if r[i][i] > r[best][best] { i } else { best });
        let col =
            Vec3::new(r[0][k] + if k == 0 { 1.0 } else { 0.0 }, r[1][k] + if k == 1 { 1.0 } else { 0.0 }, r[2][k] + if k == 2 { 1.0 } else { 0.0 });
        (col.normalized().unwrap_or(Vec3::Z), std::f64::consts::PI)
    }
}

/// Eigenvalues and eigenvectors (as columns) of a symmetric matrix, by Jacobi rotations.
pub(crate) fn jacobi<const N: usize>(mut a: [[f64; N]; N]) -> ([f64; N], [[f64; N]; N]) {
    let mut v = [[0.0; N]; N];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let scale: f64 = a.iter().flatten().map(|x| x * x).sum::<f64>().max(1e-300);
    for _ in 0..64 {
        let off: f64 = (0..N).flat_map(|i| (0..N).filter(move |j| *j != i).map(move |j| (i, j))).map(|(i, j)| a[i][j] * a[i][j]).sum();
        if off <= 1e-28 * scale {
            break;
        }
        for p in 0..N {
            for q in (p + 1)..N {
                let apq = a[p][q];
                if apq.abs() <= 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for row in a.iter_mut() {
                    let (x, y) = (row[p], row[q]);
                    row[p] = c * x - s * y;
                    row[q] = s * x + c * y;
                }
                let (rp, rq) = (a[p], a[q]);
                for (k, (x, y)) in rp.iter().zip(&rq).enumerate() {
                    a[p][k] = c * x - s * y;
                    a[q][k] = s * x + c * y;
                }
                for row in v.iter_mut() {
                    let (x, y) = (row[p], row[q]);
                    row[p] = c * x - s * y;
                    row[q] = s * x + c * y;
                }
            }
        }
    }
    let mut values = [0.0; N];
    for (i, x) in values.iter_mut().enumerate() {
        *x = a[i][i];
    }
    (values, v)
}

/// The mean of the points, `None` for none or a point that is not finite.
pub fn centroid(points: &[Vec3]) -> Option<Vec3> {
    if points.is_empty() || points.iter().any(|p| !p.is_finite()) {
        return None;
    }
    Some(points.iter().fold(Vec3::ZERO, |s, p| s + *p) / points.len() as f64)
}

/// The principal frame of a point set: its centroid and three unit axes (a right-handed frame),
/// longest spread first. The sense of each axis is the one the points' skew favours, so a part
/// that is not symmetric gets the same frame however it was scanned. `None` for fewer than
/// three points or points on a line (a plane still gives a frame; its third axis is the normal).
pub fn principal_frame(points: &[Vec3]) -> Option<(Vec3, [Vec3; 3])> {
    let c = centroid(points)?;
    if points.len() < 3 {
        return None;
    }
    let mut m = [[0.0; 3]; 3];
    for p in points {
        let d = arr(*p - c);
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += d[i] * d[j];
            }
        }
    }
    let (values, vectors) = jacobi(m);
    let mut order = [0usize, 1, 2];
    order.sort_by(|a, b| values[*b].total_cmp(&values[*a]));
    let top = values[order[0]];
    if !(top > 0.0) || values[order[1]] <= 1e-12 * top {
        return None;
    }
    let col = |k: usize| Vec3::new(vectors[0][k], vectors[1][k], vectors[2][k]);
    let mut axes = [col(order[0]), col(order[1]), col(order[2])];
    // Sense: the axis points to where the points are skewed (third moment positive).
    for a in axes.iter_mut().take(2) {
        let skew: f64 = points.iter().map(|p| (*p - c).dot(*a).powi(3)).sum();
        if skew < 0.0 {
            *a = -*a;
        }
    }
    axes[2] = axes[0].cross(axes[1]);
    Some((c, axes))
}

/// The motion that puts a point set's principal frame on the world axes (longest spread along X,
/// shortest along Z) with its centroid at `center`.
pub fn align_to_axes(points: &[Vec3], center: Vec3) -> Option<Rigid> {
    let (c, axes) = principal_frame(points)?;
    if !center.is_finite() {
        return None;
    }
    // The rows of r are the axes: r · axis_i = e_i.
    let r = [arr(axes[0]), arr(axes[1]), arr(axes[2])];
    let rig = Rigid { r, t: Vec3::ZERO };
    Some(Rigid { r, t: center - rig.apply(c) })
}

/// The rigid motion that best maps `from[i]` onto `to[i]` in least squares (Horn 1987).
/// `None` for fewer than three pairs, unequal lengths or points that are not finite.
pub fn fit_rigid(from: &[Vec3], to: &[Vec3]) -> Option<Rigid> {
    if from.len() != to.len() || from.len() < 3 {
        return None;
    }
    let (cf, ct) = (centroid(from)?, centroid(to)?);
    let mut s = [[0.0; 3]; 3];
    for (a, b) in from.iter().zip(to) {
        let (a, b) = (arr(*a - cf), arr(*b - ct));
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] += a[i] * b[j];
            }
        }
    }
    let (sxx, sxy, sxz, syx, syy, syz, szx, szy, szz) = (s[0][0], s[0][1], s[0][2], s[1][0], s[1][1], s[1][2], s[2][0], s[2][1], s[2][2]);
    let n = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];
    if n.iter().flatten().any(|x| !x.is_finite()) {
        return None;
    }
    let (values, vectors) = jacobi(n);
    let k = (0..4).fold(0, |best, i| if values[i] > values[best] { i } else { best });
    let (w, x, y, z) = (vectors[0][k], vectors[1][k], vectors[2][k], vectors[3][k]);
    let norm = (w * w + x * x + y * y + z * z).sqrt();
    if !(norm > 0.0) {
        return None;
    }
    let (w, x, y, z) = (w / norm, x / norm, y / norm, z / norm);
    let r = [
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - w * z), 2.0 * (x * z + w * y)],
        [2.0 * (x * y + w * z), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - w * x)],
        [2.0 * (x * z - w * y), 2.0 * (y * z + w * x), 1.0 - 2.0 * (x * x + y * y)],
    ];
    let rig = Rigid { r, t: Vec3::ZERO };
    Some(Rigid { r, t: ct - rig.apply(cf) })
}

/// How [`icp`] runs.
#[derive(Clone, Copy, Debug)]
pub struct IcpOptions {
    /// Most rounds.
    pub iterations: usize,
    /// Stops when a round lowers the RMS distance by less than this share of it.
    pub tolerance: f64,
    /// The share of the points (the nearest) that take part in each round, so scan noise and
    /// parts the model lacks do not pull the fit (0.5 to 1).
    pub keep: f64,
}

impl Default for IcpOptions {
    fn default() -> Self {
        IcpOptions { iterations: 50, tolerance: 1e-6, keep: 0.9 }
    }
}

/// What [`icp`] found.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IcpResult {
    /// The motion of the points from where they were given.
    pub rigid: Rigid,
    /// RMS distance of the kept points from the surface after the last round (mm).
    pub rms: f64,
    /// RMS distance of all the points, the farthest included (mm): a fit that matches nine
    /// tenths of the points and leaves the rest hanging off the model is told apart by it.
    pub rms_all: f64,
    pub iterations: usize,
}

/// Iterative closest point: moves `points` (rigidly, starting from `start`) onto the surface of
/// `target`. Each round pairs every point with its nearest point on the surface, drops the
/// farthest `1 - keep` of the pairs and fits the rest. `None` when there is nothing to fit.
pub fn icp(points: &[Vec3], target: &TriangleSet, start: Rigid, options: &IcpOptions) -> Option<IcpResult> {
    if points.len() < 3 || target.is_empty() || options.iterations == 0 {
        return None;
    }
    let keep = if options.keep.is_finite() { options.keep.clamp(0.5, 1.0) } else { 0.9 };
    let mut rigid = start;
    let mut last = f64::INFINITY;
    let mut rounds = 0;
    let mut rms = f64::INFINITY;
    for round in 0..options.iterations {
        rounds = round + 1;
        let mut pairs: Vec<(f64, Vec3, Vec3)> = Vec::with_capacity(points.len());
        for p in points {
            let q = rigid.apply(*p);
            if let Some((d, hit, _)) = target.closest_point(q, f64::INFINITY) {
                pairs.push((d, *p, hit));
            }
        }
        if pairs.len() < 3 {
            return None;
        }
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        let n = ((pairs.len() as f64 * keep).ceil() as usize).clamp(3, pairs.len());
        let kept = pairs.get(..n)?;
        rms = (kept.iter().map(|x| x.0 * x.0).sum::<f64>() / n as f64).sqrt();
        let from: Vec<Vec3> = kept.iter().map(|x| x.1).collect();
        let to: Vec<Vec3> = kept.iter().map(|x| x.2).collect();
        rigid = fit_rigid(&from, &to)?;
        if last.is_finite() && (last - rms).abs() <= options.tolerance * last.max(1e-12) {
            break;
        }
        last = rms;
    }
    // The RMS of the motion just fitted, over the kept share, for the caller to judge it by.
    let mut d2: Vec<f64> = points.iter().filter_map(|p| target.closest_point(rigid.apply(*p), f64::INFINITY)).map(|x| x.0 * x.0).collect();
    d2.sort_by(f64::total_cmp);
    let n = ((d2.len() as f64 * keep).ceil() as usize).clamp(1, d2.len().max(1));
    if let Some(kept) = d2.get(..n) {
        rms = (kept.iter().sum::<f64>() / n as f64).sqrt();
    }
    let rms_all = if d2.is_empty() { f64::INFINITY } else { (d2.iter().sum::<f64>() / d2.len() as f64).sqrt() };
    Some(IcpResult { rigid, rms, rms_all, iterations: rounds })
}

/// The signed distance of a point from a reference surface: positive outside it (on the side
/// its normal faces), negative inside.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Deviation {
    pub signed: f64,
    pub nearest: Vec3,
}

/// Signed distance of `p` from `target`'s surface, `None` when it is not finite or the surface
/// is empty.
pub fn deviation(p: Vec3, target: &TriangleSet) -> Option<Deviation> {
    let (d, q, t) = target.closest_point(p, f64::INFINITY)?;
    let n = target.normal(t)?;
    let side = (p - q).dot(n);
    // On the surface itself there is no side: the distance is zero.
    let signed = if d <= 1e-12 {
        0.0
    } else if side < 0.0 {
        -d
    } else {
        d
    };
    Some(Deviation { signed, nearest: q })
}

/// Distribution of signed deviations.
#[derive(Clone, Debug, PartialEq)]
pub struct DeviationStats {
    pub count: usize,
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    /// Root mean square of the signed values.
    pub rms: f64,
    pub std_dev: f64,
    /// Share (0–1) of the points whose distance is within the tolerance, either side.
    pub within_tolerance: f64,
    /// The median, and the 95th percentile of the absolute values.
    pub median: f64,
    pub p95_abs: f64,
    /// Points per band between `-range` and `range`, then the two outside: `bins.len()` bands
    /// of equal width, first the below-range count, last the above-range one.
    pub histogram: Vec<usize>,
}

/// Statistics of `values` (signed deviations, mm) with `bins` equal bands over `-range..range`.
/// `None` for no values or any that is not finite.
pub fn deviation_stats(values: &[f64], tolerance: f64, range: f64, bins: usize) -> Option<DeviationStats> {
    if values.is_empty() || values.iter().any(|v| !v.is_finite()) || !(range.is_finite() && range > 0.0) || bins == 0 || bins > 10_000 {
        return None;
    }
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let rms = (values.iter().map(|v| v * v).sum::<f64>() / n).sqrt();
    let std_dev = (values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n).sqrt();
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let at = |q: f64| sorted.get(((sorted.len() - 1) as f64 * q).round() as usize).copied().unwrap_or(0.0);
    let mut abs: Vec<f64> = values.iter().map(|v| v.abs()).collect();
    abs.sort_by(f64::total_cmp);
    let p95_abs = abs.get(((abs.len() - 1) as f64 * 0.95).round() as usize).copied().unwrap_or(0.0);
    let mut histogram = vec![0usize; bins + 2];
    let width = 2.0 * range / bins as f64;
    for v in values {
        let slot = if *v < -range {
            0
        } else if *v >= range {
            bins + 1
        } else {
            1 + (((*v + range) / width) as usize).min(bins - 1)
        };
        if let Some(c) = histogram.get_mut(slot) {
            *c += 1;
        }
    }
    Some(DeviationStats {
        count: values.len(),
        min: at(0.0),
        max: at(1.0),
        mean,
        rms,
        std_dev,
        within_tolerance: values.iter().filter(|v| v.abs() <= tolerance).count() as f64 / n,
        median: at(0.5),
        p95_abs,
        histogram,
    })
}

/// A colour for a deviation: blue at `-range` and below, cyan, green at 0, yellow, red at
/// `range` and above. Not finite: grey.
pub fn heat_color(value: f64, range: f64) -> [u8; 3] {
    if !(value.is_finite() && range.is_finite() && range > 0.0) {
        return [128, 128, 128];
    }
    // Position on the scale 0 (blue) to 1 (red), in four equal steps.
    let x = (value / range).clamp(-1.0, 1.0) * 0.5 + 0.5;
    let (r, g, b) = match x {
        x if x < 0.25 => (0.0, x * 4.0, 1.0),
        x if x < 0.5 => (0.0, 1.0, 1.0 - (x - 0.25) * 4.0),
        x if x < 0.75 => ((x - 0.5) * 4.0, 1.0, 0.0),
        x => (1.0, 1.0 - (x - 0.75) * 4.0, 0.0),
    };
    let byte = |c: f64| (c.clamp(0.0, 1.0) * 255.0).round() as u8;
    [byte(r), byte(g), byte(b)]
}

/// A bounding box of `points`.
pub fn bounds(points: &[Vec3]) -> Aabb3 {
    let mut b = Aabb3::EMPTY;
    for p in points {
        b.add(*p);
    }
    b
}

#[cfg(test)]
#[path = "registration_tests.rs"]
mod tests;
