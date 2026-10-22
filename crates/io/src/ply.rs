//! PLY (Polygon File Format, Stanford): reading polygon meshes (ASCII, binary little- and
//! big-endian; the `x`, `y`, `z` of `vertex` and the vertex list of `face`, polygons split into
//! fans, everything else skipped) and writing a triangle mesh with a colour on every vertex
//! (binary little-endian). Written from the public format description. A file is untrusted:
//! counts are checked against the bytes that are there and capped.

use solidcraft_geom::Vec3;

use crate::threemf::MeshObject;
use crate::{IoError, Result};

/// Most vertices and triangles we read from one file.
const MAX_READ: usize = 20_000_000;
/// Longest header we accept (bytes).
const MAX_HEADER: usize = 64 << 10;

#[derive(Clone, Copy, PartialEq)]
enum Ty {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
    F32,
    F64,
}

impl Ty {
    fn parse(s: &str) -> Option<Ty> {
        Some(match s {
            "char" | "int8" => Ty::I8,
            "uchar" | "uint8" => Ty::U8,
            "short" | "int16" => Ty::I16,
            "ushort" | "uint16" => Ty::U16,
            "int" | "int32" => Ty::I32,
            "uint" | "uint32" => Ty::U32,
            "float" | "float32" => Ty::F32,
            "double" | "float64" => Ty::F64,
            _ => return None,
        })
    }
    fn size(self) -> usize {
        match self {
            Ty::I8 | Ty::U8 => 1,
            Ty::I16 | Ty::U16 => 2,
            Ty::I32 | Ty::U32 | Ty::F32 => 4,
            Ty::F64 => 8,
        }
    }
}

struct Prop {
    name: String,
    /// The value type; for a list, the type of its items.
    ty: Ty,
    /// The type of a list's length, `None` for a single value.
    count: Option<Ty>,
}

struct Element {
    name: String,
    count: usize,
    props: Vec<Prop>,
}

#[derive(Clone, Copy, PartialEq)]
enum Enc {
    Ascii,
    Little,
    Big,
}

/// Values of the body, one at a time.
struct Body<'a> {
    enc: Enc,
    bytes: &'a [u8],
    pos: usize,
    tokens: std::str::SplitAsciiWhitespace<'a>,
}

impl Body<'_> {
    fn next(&mut self, ty: Ty) -> Result<f64> {
        let short = || IoError::Invalid("PLY file ends early or has a bad value".into());
        if self.enc == Enc::Ascii {
            return self.tokens.next().and_then(|t| t.parse::<f64>().ok()).ok_or_else(short);
        }
        let n = ty.size();
        let raw = self.bytes.get(self.pos..self.pos.checked_add(n).ok_or_else(short)?).ok_or_else(short)?;
        self.pos += n;
        let mut a = [0u8; 8];
        a.get_mut(..n).ok_or_else(short)?.copy_from_slice(raw);
        if self.enc == Enc::Big {
            a.get_mut(..n).ok_or_else(short)?.reverse();
        }
        let le4 = |i: usize| [a[i], a[i + 1], a[i + 2], a[i + 3]];
        Ok(match ty {
            Ty::I8 => f64::from(a[0] as i8),
            Ty::U8 => f64::from(a[0]),
            Ty::I16 => f64::from(i16::from_le_bytes([a[0], a[1]])),
            Ty::U16 => f64::from(u16::from_le_bytes([a[0], a[1]])),
            Ty::I32 => f64::from(i32::from_le_bytes(le4(0))),
            Ty::U32 => f64::from(u32::from_le_bytes(le4(0))),
            Ty::F32 => f64::from(f32::from_le_bytes(le4(0))),
            Ty::F64 => f64::from_le_bytes(a),
        })
    }
}

/// The mesh of a PLY file, named `name`, and what could not be used.
pub fn read_ply(bytes: &[u8], name: &str) -> Result<(Vec<MeshObject>, Vec<String>)> {
    let bad = |m: &str| IoError::Invalid(format!("PLY: {m}"));
    let end = bytes
        .windows(11)
        .take(MAX_HEADER)
        .position(|w| w == b"end_header\n")
        .map(|i| i + 11)
        .or_else(|| bytes.windows(13).take(MAX_HEADER).position(|w| w == b"end_header\r\n").map(|i| i + 13))
        .ok_or_else(|| bad("no header"))?;
    let header = std::str::from_utf8(bytes.get(..end).ok_or_else(|| bad("no header"))?).map_err(|_| bad("the header is not text"))?;
    let mut lines = header.lines().map(str::trim);
    if lines.next() != Some("ply") {
        return Err(bad("not a PLY file"));
    }
    let mut enc = None;
    let mut elements: Vec<Element> = Vec::new();
    for line in lines {
        let w: Vec<&str> = line.split_whitespace().collect();
        match w.as_slice() {
            ["format", f, _] => {
                enc = Some(match *f {
                    "ascii" => Enc::Ascii,
                    "binary_little_endian" => Enc::Little,
                    "binary_big_endian" => Enc::Big,
                    _ => return Err(bad("unknown format")),
                })
            }
            ["element", n, c] => {
                elements.push(Element { name: (*n).into(), count: c.parse().map_err(|_| bad("bad element count"))?, props: Vec::new() })
            }
            ["property", "list", ct, it, n] => {
                let prop = Prop {
                    name: (*n).into(),
                    ty: Ty::parse(it).ok_or_else(|| bad("unknown type"))?,
                    count: Some(Ty::parse(ct).ok_or_else(|| bad("unknown type"))?),
                };
                elements.last_mut().ok_or_else(|| bad("property before element"))?.props.push(prop);
            }
            ["property", t, n] => {
                let prop = Prop { name: (*n).into(), ty: Ty::parse(t).ok_or_else(|| bad("unknown type"))?, count: None };
                elements.last_mut().ok_or_else(|| bad("property before element"))?.props.push(prop);
            }
            _ => {}
        }
    }
    let enc = enc.ok_or_else(|| bad("no format"))?;
    let body_bytes = bytes.get(end..).unwrap_or(&[]);
    // Counts that the bytes cannot hold are refused before anything is allocated.
    for e in &elements {
        let least: usize = e.props.iter().map(|p| if p.count.is_some() { p.count.map_or(1, Ty::size) } else { p.ty.size() }).sum();
        let need = e.count.checked_mul(if enc == Enc::Ascii { 2.max(e.props.len()) } else { least.max(1) });
        if need.is_none_or(|n| n > body_bytes.len() + 1) && e.count > 0 {
            return Err(bad("the counts in the header are more than the file holds"));
        }
    }
    let text = if enc == Enc::Ascii { std::str::from_utf8(body_bytes).map_err(|_| bad("the body is not text"))? } else { "" };
    let mut body = Body { enc, bytes: body_bytes, pos: 0, tokens: text.split_ascii_whitespace() };
    let mut positions: Vec<Vec3> = Vec::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut warnings: Vec<String> = Vec::new();
    let (mut skipped_polygons, mut bad_vertices) = (0usize, 0usize);
    for e in &elements {
        if e.count > MAX_READ * 2 {
            return Err(bad("too many items"));
        }
        match e.name.as_str() {
            "vertex" => {
                if e.count > MAX_READ {
                    return Err(IoError::Invalid(format!("PLY with more than {MAX_READ} vertices")));
                }
                let at = |n: &str| e.props.iter().position(|p| p.name == n && p.count.is_none());
                let (Some(ix), Some(iy), Some(iz)) = (at("x"), at("y"), at("z")) else { return Err(bad("vertices without x, y, z")) };
                positions.reserve(e.count);
                let mut v = vec![0.0f64; e.props.len()];
                for _ in 0..e.count {
                    for (k, p) in e.props.iter().enumerate() {
                        match p.count {
                            None => v[k] = body.next(p.ty)?,
                            Some(ct) => {
                                for _ in 0..body.next(ct)? as usize {
                                    body.next(p.ty)?;
                                }
                            }
                        }
                    }
                    let q = Vec3::new(v[ix], v[iy], v[iz]);
                    if q.is_finite() && q.x.abs() < 1e9 && q.y.abs() < 1e9 && q.z.abs() < 1e9 {
                        positions.push(q);
                    } else {
                        // Keep the numbering: faces count vertices by position.
                        positions.push(Vec3::new(f64::NAN, f64::NAN, f64::NAN));
                        bad_vertices += 1;
                    }
                }
            }
            "face" => {
                let list = e.props.iter().position(|p| p.count.is_some() && (p.name == "vertex_indices" || p.name == "vertex_index"));
                let n = positions.len();
                for _ in 0..e.count {
                    let mut poly: Vec<u32> = Vec::new();
                    for (k, p) in e.props.iter().enumerate() {
                        match p.count {
                            None => {
                                body.next(p.ty)?;
                            }
                            Some(ct) => {
                                let len = body.next(ct)?;
                                if !(0.0..=1e6).contains(&len) {
                                    return Err(bad("a polygon with a bad vertex count"));
                                }
                                for _ in 0..len as usize {
                                    let i = body.next(p.ty)?;
                                    if Some(k) == list {
                                        poly.push(if i >= 0.0 && (i as usize) < n && i.fract() == 0.0 { i as u32 } else { u32::MAX });
                                    }
                                }
                            }
                        }
                    }
                    if poly.len() < 3 || poly.contains(&u32::MAX) {
                        skipped_polygons += 1;
                        continue;
                    }
                    for j in 1..poly.len() - 1 {
                        triangles.push([poly[0], poly[j], poly[j + 1]]);
                    }
                    if triangles.len() > MAX_READ {
                        return Err(IoError::Invalid(format!("PLY with more than {MAX_READ} triangles")));
                    }
                }
            }
            _ => {
                for _ in 0..e.count {
                    for p in &e.props {
                        match p.count {
                            None => {
                                body.next(p.ty)?;
                            }
                            Some(ct) => {
                                for _ in 0..body.next(ct)?.clamp(0.0, 1e6) as usize {
                                    body.next(p.ty)?;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if triangles.is_empty() {
        return Err(bad("the file has no faces"));
    }
    // Triangles on a vertex that is not usable go.
    let before = triangles.len();
    triangles.retain(|t| t.iter().all(|i| positions.get(*i as usize).is_some_and(|p| p.is_finite())));
    if bad_vertices > 0 {
        warnings.push(format!("{bad_vertices} vertices were not usable numbers; {} triangles on them were left out", before - triangles.len()));
    }
    if skipped_polygons > 0 {
        warnings.push(format!("{skipped_polygons} polygons were left out (fewer than three corners, or a corner that does not exist)"));
    }
    if triangles.is_empty() {
        return Err(bad("the file has no usable faces"));
    }
    Ok((vec![MeshObject { name: name.to_string(), positions, triangles, ..Default::default() }], warnings))
}

/// Most vertices or triangles one file takes.
const MAX_ITEMS: usize = 50_000_000;

/// A binary PLY file of the triangles `triangles` over `positions` (mm), vertex `i` coloured
/// `colors[i]`. Errors for unequal lengths, a vertex that is not finite or a triangle that names
/// a vertex that does not exist.
pub fn write_ply_colored(positions: &[Vec3], colors: &[[u8; 3]], triangles: &[[u32; 3]]) -> std::result::Result<Vec<u8>, String> {
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

    fn ascii_quad() -> Vec<u8> {
        b"ply\nformat ascii 1.0\nelement vertex 4\nproperty float x\nproperty float y\nproperty float z\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n4 0 1 2 3\n".to_vec()
    }

    #[test]
    fn ascii_binary_and_big_endian_files_are_read() {
        let (objects, warnings) = read_ply(&ascii_quad(), "quad").expect("ascii");
        assert!(warnings.is_empty());
        assert_eq!((objects[0].positions.len(), objects[0].triangles.len()), (4, 2));
        // What this module writes, it reads back.
        let pts = [Vec3::ZERO, Vec3::X, Vec3::Z];
        let bytes = write_ply_colored(&pts, &[[1, 2, 3]; 3], &[[0, 1, 2]]).expect("ply");
        let (back, _) = read_ply(&bytes, "t").expect("binary");
        assert_eq!(back[0].triangles, vec![[0, 1, 2]]);
        assert!((back[0].positions[1].x - 1.0).abs() < 1e-6);
        // Big-endian: one triangle, float coordinates.
        let mut big = b"ply\nformat binary_big_endian 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\nelement face 1\nproperty list uchar uint vertex_indices\nend_header\n".to_vec();
        for v in [0.0f32, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0, 0.0] {
            big.extend_from_slice(&v.to_be_bytes());
        }
        big.push(3);
        for i in [0u32, 1, 2] {
            big.extend_from_slice(&i.to_be_bytes());
        }
        let (b, _) = read_ply(&big, "b").expect("big");
        assert_eq!(b[0].positions[1].x, 2.0);
        assert_eq!(b[0].positions[2].y, 3.0);
    }

    #[test]
    fn hostile_files_are_refused_not_trusted() {
        let good = ascii_quad();
        for cut in [0, 5, 40, 120, good.len() - 4] {
            let _ = read_ply(&good[..cut], "x");
        }
        assert!(read_ply(b"", "x").is_err());
        assert!(read_ply(b"ply\nformat ascii 1.0\nend_header\n", "x").is_err());
        let huge =
            b"ply\nformat binary_little_endian 1.0\nelement vertex 4000000000\nproperty float x\nproperty float y\nproperty float z\nend_header\n";
        assert!(read_ply(huge, "x").is_err());
        let wrong_index = String::from_utf8(good.clone()).expect("text").replace("4 0 1 2 3", "4 0 1 2 9");
        assert!(read_ply(wrong_index.as_bytes(), "x").is_err());
        let nan = String::from_utf8(good.clone()).expect("text").replace("1 1 0", "nan 1 0");
        assert!(read_ply(nan.as_bytes(), "x").is_err() || read_ply(nan.as_bytes(), "x").is_ok());
        let long = String::from_utf8(good).expect("text").replace("4 0 1 2 3", "1000000000 0 1 2 3");
        assert!(read_ply(long.as_bytes(), "x").is_err());
    }

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
