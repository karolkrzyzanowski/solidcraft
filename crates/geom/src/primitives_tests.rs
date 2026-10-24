use super::*;

/// A grid-meshed rectangle in the z = 0 plane, `n` × `n` vertices.
fn plane_grid(n: usize, size: f64) -> (Vec<Vec3>, Vec<[u32; 3]>) {
    let mut pts = Vec::new();
    for i in 0..n {
        for j in 0..n {
            pts.push(Vec3::new(i as f64 * size / n as f64, j as f64 * size / n as f64, 0.0));
        }
    }
    let mut tris = Vec::new();
    for i in 0..n - 1 {
        for j in 0..n - 1 {
            let a = (i * n + j) as u32;
            let (b, c, d) = (a + 1, a + n as u32, a + n as u32 + 1);
            tris.push([a, c, b]);
            tris.push([b, c, d]);
        }
    }
    (pts, tris)
}

/// Half a cylinder of radius `r` about the x axis, `n` × `n` vertices, 20 mm long.
fn cylinder_grid(n: usize, r: f64, offset: u32) -> (Vec<Vec3>, Vec<[u32; 3]>) {
    let mut pts = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let a = std::f64::consts::PI * j as f64 / (n - 1) as f64;
            pts.push(Vec3::new(i as f64 * 20.0 / n as f64, r * a.cos(), r * a.sin() + 5.0));
        }
    }
    let mut tris = Vec::new();
    for i in 0..n - 1 {
        for j in 0..n - 1 {
            let a = (i * n + j) as u32 + offset;
            let (b, c, d) = (a + 1, a + n as u32, a + n as u32 + 1);
            tris.push([a, c, b]);
            tris.push([b, c, d]);
        }
    }
    (pts, tris)
}

#[test]
fn finds_a_plane() {
    let (p, t) = plane_grid(40, 40.0);
    let cloud = Cloud::from_mesh(p, t).unwrap();
    let found = detect(&cloud, &Options { cylinders: false, min_points: 50, ..Options::default() });
    assert_eq!(found.len(), 1);
    let Shape::Plane { normal, .. } = found[0].shape else { panic!("not a plane") };
    assert!(normal.z.abs() > 0.999);
    assert!(found[0].rms < 1e-6);
    assert_eq!(found[0].inliers.len(), 1600);
}

#[test]
fn finds_a_cylinder_and_a_plane_apart() {
    let (mut p, mut t) = plane_grid(30, 30.0);
    for v in &mut p {
        v.z -= 20.0;
    }
    let offset = p.len() as u32;
    let (cp, ct) = cylinder_grid(30, 8.0, offset);
    p.extend(cp);
    t.extend(ct);
    let cloud = Cloud::from_mesh(p, t).unwrap();
    let found = detect(&cloud, &Options { min_points: 100, ..Options::default() });
    let cyl = found.iter().find(|f| matches!(f.shape, Shape::Cylinder { .. })).expect("cylinder");
    let Shape::Cylinder { radius, axis, .. } = cyl.shape else { panic!() };
    assert!((radius - 8.0).abs() < 0.05, "radius {radius}");
    assert!(axis.x.abs() > 0.999);
    assert!(found.iter().any(|f| matches!(f.shape, Shape::Plane { .. })));
}

#[test]
fn hostile_input_gives_nothing() {
    let (p, t) = plane_grid(10, 10.0);
    let cloud = Cloud::from_mesh(p, t).unwrap();
    for opts in [
        Options { tolerance: f64::NAN, ..Options::default() },
        Options { tolerance: -1.0, ..Options::default() },
        Options { angle: 4.0, ..Options::default() },
        Options { trials: 0, ..Options::default() },
        Options { min_points: 0, ..Options::default() },
        Options { planes: false, cylinders: false, ..Options::default() },
    ] {
        assert!(detect(&cloud, &opts).is_empty());
    }
    assert!(Cloud::from_mesh(vec![Vec3::new(f64::NAN, 0.0, 0.0)], vec![]).is_none());
    assert!(Cloud::from_mesh(vec![Vec3::ZERO], vec![[0, 1, 2]]).is_none());
    assert!(detect(&Cloud::default(), &Options::default()).is_empty());
    let bad = Cloud { points: vec![Vec3::ZERO; 5], normals: vec![Vec3::Z; 2], triangles: vec![] };
    assert!(detect(&bad, &Options::default()).is_empty());
}

#[test]
fn decimation_merges_cells_and_refuses_bad_input() {
    let (p, t) = plane_grid(20, 20.0);
    let (q, u) = cluster_decimate(&p, &t, 4.0).unwrap();
    assert!(q.len() < p.len() && !u.is_empty());
    assert!(cluster_decimate(&p, &t, 0.0).is_none());
    assert!(cluster_decimate(&p, &t, f64::NAN).is_none());
    assert!(cluster_decimate(&p, &[[0, 1, 9999]], 1.0).is_none());
    assert!(cluster_decimate(&[Vec3::new(f64::INFINITY, 0.0, 0.0)], &[], 1.0).is_none());
}

#[test]
fn a_plane_split_in_two_is_one_plane_and_a_nearly_flat_arc_is_not_a_cylinder() {
    let (p, t) = plane_grid(30, 30.0);
    let cloud = Cloud::from_mesh(p, t).unwrap();
    let mk = |idx: Vec<u32>| Found { shape: refit_plane(&cloud, &idx).unwrap(), rms: 0.0, inliers: idx };
    let left = mk((0..450).collect());
    let right = mk((450..900).collect());
    let merged = merge(&cloud, vec![left, right], 0.1);
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0].inliers.len(), 900);
    // A 30 mm wide strip of a cylinder 500 mm in radius bulges 0.225 mm.
    let flat = Found {
        shape: Shape::Cylinder { origin: Vec3::new(0.0, 0.0, -500.0), axis: Vec3::new(1.0, 0.0, 0.0), radius: 500.0, min: 0.0, max: 30.0 },
        inliers: (0..900).collect(),
        rms: 0.0,
    };
    assert!((sagitta(&cloud, &flat) - 0.21).abs() < 0.02, "{}", sagitta(&cloud, &flat));
}

#[test]
fn a_face_is_solid_and_scattered_points_are_not() {
    let (p, t) = plane_grid(30, 30.0);
    let adj = Adjacency::new(p.len(), &t);
    let face: Vec<u32> = (0..p.len() as u32).filter(|i| (i / 30) < 15).collect();
    assert!(solidity(&adj, p.len(), &face) > 0.9);
    // Every other row: each vertex has only its own row's neighbours in the patch.
    let rows: Vec<u32> = (0..p.len() as u32).filter(|i| (i / 30) % 2 == 0).collect();
    assert!(solidity(&adj, p.len(), &rows) < 0.6, "{}", solidity(&adj, p.len(), &rows));
}
