//! `mesh.detect`: planes and cylinders in a scanned mesh body.

use serde_json::{Value, json};
use solidcraft_geom::primitives::{self, Cloud, Options, Shape};
use solidcraft_geom::{Vec3, registration};

use super::CommandSpec;
use super::mesh_fit::{array, number, pair, surface, weld};
use crate::params::{bad, bool_};
use crate::{Result, Session};

pub static COMMANDS: &[CommandSpec] = &[CommandSpec::new("mesh.detect", "Detect Shapes", detect).icon("compare").noundo().params(
    "body?: a mesh body (default: the selected mesh body); planes? (default true); cylinders? (default true); tolerance?: mm (0.001 to 100, default 0.2), \
     how far a vertex may lie from a shape; angle?: degrees (1 to 60, default 15), how far its normal may turn; min_points? (default 200); \
     max_shapes? (1 to 100, default 20); trials? (10 to 5000, default 400) guesses per shape; seed?; save?: path of a .ply file with each shape in its own \
     colour (grey where nothing was found) → shapes [{kind: plane | cylinder, points, share, rms_mm, plane: {origin, normal} | cylinder: {origin, axis, \
     radius_mm, from_mm, to_mm}}], covered (share of the vertices on a shape), vertices (used), verdict, ply {path, bytes}",
)];

/// Most vertices searched; a finer scan is merged cell by cell first.
const MAX_VERTICES: usize = 150_000;

const PALETTE: [[u8; 3]; 8] =
    [[230, 80, 80], [80, 160, 230], [90, 200, 110], [240, 190, 60], [170, 100, 220], [60, 200, 200], [240, 130, 40], [200, 100, 160]];

fn detect(s: &mut Session, p: &Value) -> Result<Value> {
    let cmd = "mesh.detect";
    let (scan, _) = pair(s, p, cmd, false)?;
    let planes = bool_(p, "planes").unwrap_or(true);
    let cylinders = bool_(p, "cylinders").unwrap_or(true);
    if !(planes || cylinders) {
        return Err(bad(cmd, "`planes` and `cylinders` are both off"));
    }
    let tolerance = number(cmd, p, "tolerance", 0.2, 0.001..=100.0)?;
    let angle = number(cmd, p, "angle", 15.0, 1.0..=60.0)?.to_radians();
    let min_points = number(cmd, p, "min_points", 200.0, 3.0..=1e7)? as usize;
    let max_shapes = number(cmd, p, "max_shapes", 20.0, 1.0..=100.0)? as usize;
    let trials = number(cmd, p, "trials", 400.0, 10.0..=5000.0)? as usize;
    let seed = number(cmd, p, "seed", 1.0, 0.0..=1e12)? as u64;
    let save = match p.get("save") {
        None | Some(Value::Null) => None,
        Some(Value::String(path)) if !path.trim().is_empty() && path.len() < 4096 && path.to_ascii_lowercase().ends_with(".ply") => {
            Some(path.clone())
        }
        Some(_) => return Err(bad(cmd, "`save` must be the path of a .ply file")),
    };
    let (mut points, mut triangles) = weld(&surface(s, &scan, 0.05, cmd)?);
    let original = points.len();
    if points.len() > MAX_VERTICES {
        let diag = registration::bounds(&points).diagonal();
        let mut cell = diag / 400.0;
        for _ in 0..12 {
            let (q, t) = primitives::cluster_decimate(&points, &triangles, cell).ok_or_else(|| bad(cmd, "the scan could not be thinned out"))?;
            if q.len() <= MAX_VERTICES || cell <= 0.0 {
                points = q;
                triangles = t;
                break;
            }
            cell *= 1.5;
        }
    }
    let used = points.len();
    let cloud = Cloud::from_mesh(points, triangles).ok_or_else(|| bad(cmd, "the scan has vertices that are not finite"))?;
    let found = primitives::detect(&cloud, &Options { planes, cylinders, tolerance, angle, min_points, max_shapes, trials, seed });
    let on: usize = found.iter().map(|f| f.inliers.len()).sum();
    let shapes: Vec<Value> = found
        .iter()
        .map(|f| {
            let mut v = json!({
                "kind": match f.shape { Shape::Plane { .. } => "plane", Shape::Cylinder { .. } => "cylinder" },
                "points": f.inliers.len(),
                "share": f.inliers.len() as f64 / used.max(1) as f64,
                "rms_mm": f.rms,
            });
            match f.shape {
                Shape::Plane { origin, normal } => v["plane"] = json!({"origin": array(origin), "normal": array(normal)}),
                Shape::Cylinder { origin, axis, radius, min, max } => {
                    v["cylinder"] = json!({"origin": array(origin), "axis": array(axis), "radius_mm": radius, "from_mm": min, "to_mm": max});
                }
            }
            v
        })
        .collect();
    let verdict = if found.is_empty() {
        format!("`{scan}`: no shapes found (tolerance {tolerance} mm; try a larger one or fewer min_points)")
    } else {
        let planes_n = found.iter().filter(|f| matches!(f.shape, Shape::Plane { .. })).count();
        format!(
            "`{scan}`: {} planes, {} cylinders; {:.1}% of {} vertices lie on them",
            planes_n,
            found.len() - planes_n,
            100.0 * on as f64 / used.max(1) as f64,
            used
        )
    };
    let mut out = json!({
        "scan": scan,
        "shapes": shapes,
        "covered": on as f64 / used.max(1) as f64,
        "vertices": used,
        "thinned": used < original,
        "verdict": verdict,
    });
    if let Some(path) = save {
        let mut colors = vec![[150u8, 150, 150]; cloud.points.len()];
        for (k, f) in found.iter().enumerate() {
            let c = PALETTE[k % PALETTE.len()];
            for i in &f.inliers {
                if let Some(slot) = colors.get_mut(*i as usize) {
                    *slot = c;
                }
            }
        }
        let positions: Vec<Vec3> = cloud.points.clone();
        let bytes = solidcraft_io::write_ply_colored(&positions, &colors, &cloud.triangles).map_err(|e| bad(cmd, e))?;
        solidcraft_io::vfs::write(&path, &bytes).map_err(|e| crate::EngineError::Other(format!("{path}: {e}")))?;
        out["ply"] = json!({"path": path, "bytes": bytes.len()});
    }
    if bool_(p, "echo").unwrap_or(true) {
        s.echo(verdict);
    }
    Ok(out)
}
