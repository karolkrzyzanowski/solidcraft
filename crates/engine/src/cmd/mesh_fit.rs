//! Reverse engineering from a scan: place a scanned mesh (on the world axes by its principal
//! frame, or onto a model by iterative closest point) and measure how far it lies from a
//! reference body, as numbers and as a colour map saved with the scan's vertices.

use std::collections::HashMap;

use serde_json::{Value, json};
use solidcraft_doc::{Dir, FeatureKind};
use solidcraft_geom::mesh_distance::TriangleSet;
use solidcraft_geom::registration::{self, IcpOptions, Rigid};
use solidcraft_geom::{Mesh, Vec3};

use super::{CommandSpec, features::add_feature};
use crate::params::{bad, bool_, str_, vec3};
use crate::{Result, Sel, Session};

pub static COMMANDS: &[CommandSpec] = &[
    CommandSpec::new("mesh.align", "Align Mesh", align).icon("move").params(
        "body?: mesh body (default: the selected mesh body); to?: reference body (default: the other selected body): fit the scan onto it (iterative closest point, \
         started from the best of the scan's principal frames and its current place); without `to` the scan's principal axes are put on the world axes \
         (longest along X, shortest along Z) with its centre of mass at `center?` ([x,y,z] mm, default the origin); initial?: auto | current (fit only: \
         auto tries the principal frames, current starts where the scan is); iterations? (default 50, up to 500); keep?: share of the closest points \
         that take part in each round (0.5 to 1, default 0.9: noise and features the model lacks do not pull the fit); max_points? (default 20000, up to \
         200000); tolerance?: mm (0.001 to 1, default 0.02), how finely a curved reference is meshed; name? → a Move feature (timeline step, undoable), \
         rotation {axis, degrees}, translation, rms_mm (fit), candidate (fit), warnings",
    ),
    CommandSpec::new("mesh.deviation", "Mesh Deviation", deviation).icon("compare").noundo().params(
        "body?: scan, a mesh body (default: the selected mesh body); to?: reference body (default: the other selected body); tolerance?: mm (default 0.1), \
         the band counted as good; range?: mm (default 1), where the colour scale and the histogram end; bins? (2 to 200, default 20); chord?: mm (default \
         0.02), how finely a curved reference is meshed; max_points? (default 500000, up to 3000000) vertices measured; save?: path of a .ply file with the \
         scan coloured by its signed deviation (blue inside the reference, green on it, red outside; every vertex is measured) → count, min_mm, max_mm, \
         mean_mm, rms_mm, std_dev_mm, median_mm, p95_abs_mm, within_tolerance (share), worst_outside {mm, at}, worst_inside {mm, at}, histogram {range_mm, \
         bins (counts below, in each band, above)}, verdict, ply {path, bytes}. Positive is outside the reference",
    ),
];

const MAX_ITERATIONS: usize = 500;

fn number(cmd: &str, p: &Value, key: &str, default: f64, range: std::ops::RangeInclusive<f64>) -> Result<f64> {
    match p.get(key) {
        None | Some(Value::Null) => Ok(default),
        Some(v) => v
            .as_f64()
            .filter(|x| x.is_finite() && range.contains(x))
            .ok_or_else(|| bad(cmd, format!("`{key}` must be a number from {} to {}", range.start(), range.end()))),
    }
}

fn name_param(cmd: &str, p: &Value, key: &str) -> Result<Option<String>> {
    match p.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(n)) if !n.trim().is_empty() => Ok(Some(n.clone())),
        Some(_) => Err(bad(cmd, format!("`{key}` must be a body name"))),
    }
}

/// The scan and the reference: `body` and `to`, else the selected bodies (the mesh body is the
/// scan; with two selected, the other is the reference).
fn pair(s: &Session, p: &Value, cmd: &str, need_reference: bool) -> Result<(String, Option<String>)> {
    let st = s.world_state();
    let selected: Vec<String> = s.selection.iter().filter_map(|x| if let Sel::Body { name } = x { Some(name.clone()) } else { None }).collect();
    let is_mesh = |n: &str| st.body(n).is_some_and(|b| b.body.is_mesh());
    let scan = match name_param(cmd, p, "body")? {
        Some(n) => n,
        None => {
            let meshes: Vec<&String> = selected.iter().filter(|n| is_mesh(n)).collect();
            match meshes.as_slice() {
                [one] => (*one).clone(),
                [] => return Err(bad(cmd, "select a mesh body (or give `body`)")),
                [first, ..] if need_reference && selected.len() == 2 => (*first).clone(),
                _ => return Err(bad(cmd, "select one mesh body")),
            }
        }
    };
    let reference = match name_param(cmd, p, "to")? {
        Some(n) => Some(n),
        None => {
            let others: Vec<&String> = selected.iter().filter(|n| **n != scan).collect();
            match others.as_slice() {
                [one] => Some((*one).clone()),
                [] if !need_reference => None,
                [] => return Err(bad(cmd, "select the reference body too (or give `to`)")),
                _ => return Err(bad(cmd, "select the scan and one reference body")),
            }
        }
    };
    let Some(b) = st.body(&scan) else { return Err(bad(cmd, format!("no body `{scan}`"))) };
    if !b.body.is_mesh() {
        return Err(bad(cmd, format!("`{scan}` is not a mesh body")));
    }
    if let Some(r) = &reference {
        if *r == scan {
            return Err(bad(cmd, "the scan and the reference must be two different bodies"));
        }
        if st.body(r).is_none() {
            return Err(bad(cmd, format!("no body `{r}`")));
        }
    }
    Ok((scan, reference))
}

/// A body's surface as a mesh, in world placement.
fn surface(s: &Session, name: &str, chord: f64, cmd: &str) -> Result<Mesh> {
    let st = s.world_state();
    let b = st.body(name).ok_or_else(|| bad(cmd, format!("no body `{name}`")))?;
    let mesh = b.body.tessellate_with(chord, true).map_err(|e| bad(cmd, format!("`{name}`: {e}")))?;
    if mesh.triangles.is_empty() {
        return Err(bad(cmd, format!("`{name}` has no surface")));
    }
    Ok(mesh)
}

/// The distinct vertices of a mesh whose corners are split per triangle, with the triangles
/// rewritten to them.
fn weld(mesh: &Mesh) -> (Vec<Vec3>, Vec<[u32; 3]>) {
    let mut index: HashMap<[u64; 3], u32> = HashMap::new();
    let mut points: Vec<Vec3> = Vec::new();
    let mut map = |p: Vec3| -> Option<u32> {
        let key = [(p.x + 0.0).to_bits(), (p.y + 0.0).to_bits(), (p.z + 0.0).to_bits()];
        if let Some(i) = index.get(&key) {
            return Some(*i);
        }
        let i = u32::try_from(points.len()).ok()?;
        points.push(p);
        index.insert(key, i);
        Some(i)
    };
    let mut triangles = Vec::with_capacity(mesh.triangles.len());
    for t in &mesh.triangles {
        let corners: Option<Vec<u32>> = t.iter().map(|i| mesh.positions.get(*i as usize).copied().and_then(&mut map)).collect();
        if let Some([a, b, c]) = corners.as_deref() {
            triangles.push([*a, *b, *c]);
        }
    }
    (points, triangles)
}

/// At most `n` of the points, evenly spread through the list.
fn sample(points: &[Vec3], n: usize) -> Vec<Vec3> {
    if points.len() <= n || n == 0 {
        return points.to_vec();
    }
    let step = points.len() as f64 / n as f64;
    (0..n).filter_map(|i| points.get((i as f64 * step) as usize).copied()).collect()
}

fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// The starting motions for a fit: the scan's centre on the reference's, turned each of the four
/// proper ways that put its principal axes on the reference's.
fn starts(scan: &[Vec3], reference: &[Vec3]) -> Vec<Rigid> {
    let (Some(cs), Some(cr)) = (registration::centroid(scan), registration::centroid(reference)) else { return Vec::new() };
    let shift = Rigid { t: cr - cs, ..Rigid::IDENTITY };
    let mut out = vec![shift];
    if let (Some((_, a)), Some((_, b))) = (registration::principal_frame(scan), registration::principal_frame(reference)) {
        for signs in [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]] {
            let mut r = [[0.0; 3]; 3];
            for (i, row) in r.iter_mut().enumerate() {
                for (j, cell) in row.iter_mut().enumerate() {
                    *cell = (0..3).map(|k| array(b[k])[i] * signs[k] * array(a[k])[j]).sum();
                }
            }
            let turn = Rigid { r, t: Vec3::ZERO };
            out.push(Rigid { r, t: cr - turn.apply(cs) });
        }
    }
    out
}

fn motion_json(m: &Rigid) -> Value {
    let (axis, angle) = m.axis_angle();
    json!({"rotation": {"axis": array(axis), "degrees": angle.to_degrees()}, "translation": array(m.t)})
}

fn align(s: &mut Session, p: &Value) -> Result<Value> {
    let cmd = "mesh.align";
    let (scan, reference) = pair(s, p, cmd, false)?;
    let st = s.model.state();
    let source = st.body(&scan).ok_or_else(|| bad(cmd, format!("no body `{scan}`")))?;
    if s.doc.model_body_component(source) != s.active_component {
        return Err(bad(cmd, "activate the mesh body's component before aligning it"));
    }
    let iterations = number(cmd, p, "iterations", 50.0, 1.0..=MAX_ITERATIONS as f64)? as usize;
    let keep = number(cmd, p, "keep", 0.9, 0.5..=1.0)?;
    let max_points = number(cmd, p, "max_points", 20_000.0, 100.0..=200_000.0)? as usize;
    let chord = number(cmd, p, "tolerance", 0.02, 0.001..=1.0)?;
    let initial = str_(p, "initial").unwrap_or("auto");
    if !matches!(initial, "auto" | "current") {
        return Err(bad(cmd, "`initial` must be auto or current"));
    }
    let center = match p.get("center") {
        None | Some(Value::Null) => Vec3::ZERO,
        Some(v) => vec3(v).filter(|c| c.is_finite()).ok_or_else(|| bad(cmd, "`center` must be [x, y, z] in millimetres"))?,
    };
    let points = {
        let (verts, _) = weld(&surface(s, &scan, 0.05, cmd)?);
        sample(&verts, max_points)
    };
    let mut warnings: Vec<String> = Vec::new();
    let (motion, extra) = match &reference {
        None => {
            let m = registration::align_to_axes(&points, center)
                .ok_or_else(|| bad(cmd, "the mesh has no principal axes (too few points, or all on a line)"))?;
            (m, json!({"mode": "axes"}))
        }
        Some(r) => {
            let model = surface(s, r, chord, cmd)?;
            let target = TriangleSet::new(&model);
            let (ref_points, _) = weld(&model);
            let candidates = if initial == "current" { vec![Rigid::IDENTITY] } else { starts(&points, &ref_points) };
            if candidates.is_empty() {
                return Err(bad(cmd, "nothing to fit"));
            }
            // Every start is tried on a thin sample; the best is then refined on the whole one.
            let coarse = sample(&points, 2000);
            let rough = IcpOptions { iterations: iterations.min(30), keep, ..IcpOptions::default() };
            let mut best: Option<(usize, registration::IcpResult)> = None;
            for (i, start) in candidates.iter().enumerate() {
                if let Some(found) = registration::icp(&coarse, &target, *start, &rough)
                    && best.as_ref().is_none_or(|(_, b)| found.rms < b.rms)
                {
                    best = Some((i, found));
                }
            }
            let (index, rough_fit) = best.ok_or_else(|| bad(cmd, "the scan could not be fitted onto the reference"))?;
            let fine = IcpOptions { iterations, keep, ..IcpOptions::default() };
            let found = registration::icp(&points, &target, rough_fit.rigid, &fine).ok_or_else(|| bad(cmd, "the scan could not be fitted onto the reference"))?;
            let diag = registration::bounds(&ref_points).diagonal();
            if found.rms > 0.05 * diag {
                warnings.push(format!(
                    "the fit is poor (RMS {:.3} mm against a reference {:.1} mm across): the scan may be another part, or too different in shape; try `initial: \"current\"` after placing it by hand",
                    found.rms, diag
                ));
            }
            (found.rigid, json!({"mode": "fit", "reference": r, "rms_mm": found.rms, "candidate": index, "candidates": candidates.len(), "iterations": found.iterations}))
        }
    };
    if !motion.t.is_finite() || motion.r.iter().flatten().any(|x| !x.is_finite()) {
        return Err(bad(cmd, "the placement is not finite"));
    }
    let (axis, angle) = motion.axis_angle();
    let num = |x: f64| format!("{x:.12}");
    let mut feature = add_feature(
        s,
        p,
        FeatureKind::Move {
            bodies: vec![scan.clone().into()],
            translate: [format!("{} mm", num(motion.t.x)), format!("{} mm", num(motion.t.y)), format!("{} mm", num(motion.t.z))],
            rotate_axis: Some(Dir(axis)),
            angle: Some(format!("{} rad", num(angle))),
            copy: false,
        },
    )?;
    if let (Some(out), Some(info)) = (feature.as_object_mut(), extra.as_object()) {
        out.extend(motion_json(&motion).as_object().cloned().unwrap_or_default());
        out.extend(info.clone());
        out.insert("warnings".into(), json!(warnings));
    }
    Ok(feature)
}

fn deviation(s: &mut Session, p: &Value) -> Result<Value> {
    let cmd = "mesh.deviation";
    let (scan, reference) = pair(s, p, cmd, true)?;
    let Some(reference) = reference else { return Err(bad(cmd, "give the reference body in `to`")) };
    let tolerance = number(cmd, p, "tolerance", 0.1, 1e-6..=1e6)?;
    let range = number(cmd, p, "range", 1.0, 1e-6..=1e6)?;
    let bins = number(cmd, p, "bins", 20.0, 2.0..=200.0)? as usize;
    let chord = number(cmd, p, "chord", 0.02, 0.001..=1.0)?;
    let max_points = number(cmd, p, "max_points", 500_000.0, 100.0..=3_000_000.0)? as usize;
    let save = match p.get("save") {
        None | Some(Value::Null) => None,
        Some(Value::String(path)) if !path.trim().is_empty() && path.len() < 4096 && path.to_ascii_lowercase().ends_with(".ply") => Some(path.clone()),
        Some(_) => return Err(bad(cmd, "`save` must be the path of a .ply file")),
    };
    let model = TriangleSet::new(&surface(s, &reference, chord, cmd)?);
    let (verts, triangles) = weld(&surface(s, &scan, 0.05, cmd)?);
    if save.is_some() && verts.len() > max_points {
        return Err(bad(cmd, format!("the scan has {} vertices; raise `max_points` to at least that to save its colour map", verts.len())));
    }
    let measured = sample(&verts, max_points);
    let devs: Vec<Option<registration::Deviation>> = if save.is_some() {
        measured.iter().map(|v| registration::deviation(*v, &model)).collect()
    } else {
        // Large scans are measured on several threads.
        let chunks: Vec<&[Vec3]> = measured.chunks(20_000).collect();
        solidcraft_kernel::map_on_threads(&chunks, 16, |c| c.iter().map(|v| registration::deviation(*v, &model)).collect::<Vec<_>>())
            .into_iter()
            .flat_map(|r| r.unwrap_or_default())
            .collect()
    };
    let values: Vec<f64> = devs.iter().flatten().map(|d| d.signed).collect();
    let stats = registration::deviation_stats(&values, tolerance, range, bins).ok_or_else(|| bad(cmd, "nothing could be measured"))?;
    let worst = |sign: f64| {
        devs.iter()
            .zip(&measured)
            .filter_map(|(d, v)| d.map(|d| (d.signed * sign, *v)))
            .max_by(|a, b| a.0.total_cmp(&b.0))
            .map(|(d, at)| json!({"mm": d * sign, "at": array(at)}))
    };
    let verdict = format!(
        "{} points: {:.1}% within ±{} mm of `{reference}`; RMS {:.4} mm, mean {:+.4} mm, range {:+.4} to {:+.4} mm",
        stats.count,
        stats.within_tolerance * 100.0,
        tolerance,
        stats.rms,
        stats.mean,
        stats.min,
        stats.max
    );
    let mut out = json!({
        "scan": scan,
        "reference": reference,
        "count": stats.count,
        "min_mm": stats.min,
        "max_mm": stats.max,
        "mean_mm": stats.mean,
        "rms_mm": stats.rms,
        "std_dev_mm": stats.std_dev,
        "median_mm": stats.median,
        "p95_abs_mm": stats.p95_abs,
        "tolerance_mm": tolerance,
        "within_tolerance": stats.within_tolerance,
        "worst_outside": worst(1.0),
        "worst_inside": worst(-1.0),
        "histogram": {"range_mm": range, "bins": stats.histogram},
        "sampled": measured.len() < verts.len(),
        "verdict": verdict,
    });
    if let Some(path) = save {
        let colors: Vec<[u8; 3]> = devs.iter().map(|d| registration::heat_color(d.map_or(f64::NAN, |d| d.signed), range)).collect();
        let bytes = solidcraft_io::write_ply_colored(&measured, &colors, &triangles).map_err(|e| bad(cmd, e))?;
        solidcraft_io::vfs::write(&path, &bytes).map_err(|e| crate::EngineError::Other(format!("{path}: {e}")))?;
        out["ply"] = json!({"path": path, "bytes": bytes.len()});
    }
    if bool_(p, "echo").unwrap_or(true) {
        s.echo(verdict_line(&out));
    }
    Ok(out)
}

fn verdict_line(out: &Value) -> String {
    out.get("verdict").and_then(Value::as_str).unwrap_or_default().to_owned()
}

#[cfg(test)]
#[path = "mesh_fit_tests.rs"]
mod tests;
