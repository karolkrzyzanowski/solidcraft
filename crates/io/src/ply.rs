//! PLY (Polygon File Format, Stanford) output: a triangle mesh with a colour on every vertex,
//! binary little-endian. Written from the public format description.

use solidcraft_geom::Vec3;

/// Most vertices or triangles one file takes.
const MAX_ITEMS: usize = 50_000_000;

/// A binary PLY file of the triangles `triangles` over `positions` (mm), vertex `i` coloured
/// `colors[i]`. Errors for unequal lengths, a vertex that is not finite or a triangle that names
/// a vertex that does not exist.
pub fn write_ply_colored(positions: &[Vec3], colors: &[[u8; 3]], triangles: &[[u32; 3]]) -> Result<Vec<u8>, String> {
    if positions.len() != colors.len() {
        return Err("a colour is needed for every vertex".into());
    }
    if positions.len() > MAX_ITEMS || triangles.len() > MAX_ITEMS {
        return Err("mesh too large for a PLY file".into());
    }
    if positions.iter().any(|p| !p.is_finite()) {
        return Err("a vertex is not finite".into());
    }
    let n = u32::try_from(positions.len()).map_err(|_| "too many vertices")?;
    if triangles.iter().flatten().any(|i| *i >= n) {
        return Err("a triangle names a vertex that does not exist".into());
    }
    let header = format!(
        "ply\nformat binary_little_endian 1.0\ncomment SolidCraft deviation map\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\n\
         property uchar red\nproperty uchar green\nproperty uchar blue\nelement face {}\nproperty list uchar uint vertex_indices\nend_header\n",
        positions.len(),
        triangles.len()
    );
    let mut out = Vec::with_capacity(header.len() + positions.len() * 15 + triangles.len() * 13);
    out.extend_from_slice(header.as_bytes());
    for (p, c) in positions.iter().zip(colors) {
        for v in [p.x, p.y, p.z] {
            out.extend_from_slice(&(v as f32).to_le_bytes());
        }
        out.extend_from_slice(c);
    }
    for t in triangles {
        out.push(3);
        for i in t {
            out.extend_from_slice(&i.to_le_bytes());
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_triangle_is_written_and_bad_input_refused() {
        let pts = [Vec3::ZERO, Vec3::X, Vec3::Z];
        let cols = [[255, 0, 0], [0, 255, 0], [0, 0, 255]];
        let bytes = write_ply_colored(&pts, &cols, &[[0, 1, 2]]).expect("ply");
        let text = String::from_utf8_lossy(&bytes[..bytes.len().min(300)]).into_owned();
        assert!(text.starts_with("ply\nformat binary_little_endian 1.0\n"));
        assert!(text.contains("element vertex 3") && text.contains("element face 1"));
        let header_end = bytes.windows(11).position(|w| w == b"end_header\n").map(|i| i + 11).unwrap_or(0);
        assert_eq!(bytes.len(), header_end + 3 * 15 + 13);
        assert!(write_ply_colored(&pts, &cols[..2], &[]).is_err());
        assert!(write_ply_colored(&pts, &cols, &[[0, 1, 3]]).is_err());
        assert!(write_ply_colored(&[Vec3::new(f64::NAN, 0.0, 0.0)], &cols[..1], &[]).is_err());
    }
}
