use super::*;
use crate::Mesh;

/// A box from the origin to `(x, y, z)`, triangles counter-clockwise seen from outside.
fn box_mesh(x: f64, y: f64, z: f64) -> Mesh {
    let positions = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(x, 0.0, 0.0),
        Vec3::new(x, y, 0.0),
        Vec3::new(0.0, y, 0.0),
        Vec3::new(0.0, 0.0, z),
        Vec3::new(x, 0.0, z),
        Vec3::new(x, y, z),
        Vec3::new(0.0, y, z),
    ];
    let triangles = vec![
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ];
    let tri_face = vec![0; triangles.len()];
    Mesh { positions, triangles, tri_face, ..Mesh::default() }
}

/// Points over the box's surface, denser toward one corner so the cloud is not symmetric.
fn surface_points(x: f64, y: f64, z: f64) -> Vec<Vec3> {
    let mut out = Vec::new();
    let steps = 14;
    for i in 0..=steps {
        for j in 0..=steps {
            let (u, v) = (f64::from(i) / f64::from(steps), f64::from(j) / f64::from(steps));
            let (u, v) = (u * u, v * v);
            out.push(Vec3::new(u * x, v * y, 0.0));
            out.push(Vec3::new(u * x, v * y, z));
            out.push(Vec3::new(u * x, 0.0, v * z));
            out.push(Vec3::new(u * x, y, v * z));
            out.push(Vec3::new(0.0, u * y, v * z));
            out.push(Vec3::new(x, u * y, v * z));
        }
    }
    out
}

fn rotation(axis: Vec3, angle: f64) -> Rigid {
    let a = axis.normalized().unwrap_or(Vec3::Z);
    let (c, s) = (angle.cos(), angle.sin());
    let t = 1.0 - c;
    Rigid {
        r: [
            [t * a.x * a.x + c, t * a.x * a.y - s * a.z, t * a.x * a.z + s * a.y],
            [t * a.x * a.y + s * a.z, t * a.y * a.y + c, t * a.y * a.z - s * a.x],
            [t * a.x * a.z - s * a.y, t * a.y * a.z + s * a.x, t * a.z * a.z + c],
        ],
        t: Vec3::ZERO,
    }
}

fn close(a: Vec3, b: Vec3, tol: f64) -> bool {
    a.dist(b) <= tol
}

#[test]
fn fit_rigid_recovers_a_known_motion() {
    let pts = surface_points(40.0, 30.0, 20.0);
    let mut motion = rotation(Vec3::new(1.0, 2.0, 3.0), 0.7);
    motion.t = Vec3::new(5.0, -7.0, 11.0);
    let moved: Vec<Vec3> = pts.iter().map(|p| motion.apply(*p)).collect();
    let fit = fit_rigid(&pts, &moved).expect("fit");
    assert!(pts.iter().zip(&moved).all(|(p, q)| close(fit.apply(*p), *q, 1e-9)));
}

#[test]
fn align_to_axes_puts_the_long_side_on_x_and_the_centroid_at_the_centre() {
    let pts = surface_points(40.0, 30.0, 20.0);
    let scrambled = {
        let mut m = rotation(Vec3::new(0.3, -1.0, 0.8), 2.1);
        m.t = Vec3::new(100.0, 50.0, -20.0);
        m
    };
    let scan: Vec<Vec3> = pts.iter().map(|p| scrambled.apply(*p)).collect();
    let motion = align_to_axes(&scan, Vec3::ZERO).expect("frame");
    let aligned: Vec<Vec3> = scan.iter().map(|p| motion.apply(*p)).collect();
    assert!(close(centroid(&aligned).expect("centroid"), Vec3::ZERO, 1e-9));
    // The spread is now along the axes, longest on X and shortest on Z.
    let mut cov = [[0.0f64; 3]; 3];
    for p in &aligned {
        let d = [p.x, p.y, p.z];
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += d[i] * d[j];
            }
        }
    }
    let scale = cov[0][0];
    assert!(cov[0][0] > cov[1][1] && cov[1][1] > cov[2][2]);
    assert!(cov[0][1].abs() < 1e-9 * scale && cov[0][2].abs() < 1e-9 * scale && cov[1][2].abs() < 1e-9 * scale, "{cov:?}");
    // The same part scanned in another pose lands in the same place.
    let other = rotation(Vec3::new(-1.0, 0.2, 0.4), -1.3);
    let scan2: Vec<Vec3> = pts.iter().map(|p| other.apply(*p)).collect();
    let motion2 = align_to_axes(&scan2, Vec3::ZERO).expect("frame");
    let again: Vec<Vec3> = scan2.iter().map(|p| motion2.apply(*p)).collect();
    // Both are the original cloud (about its centroid) in the principal frame.
    let reference: Vec<Vec3> = scan.iter().map(|p| motion.apply(*p)).collect();
    assert!(again.iter().zip(&reference).all(|(a, b)| close(*a, *b, 1e-6)));
}

#[test]
fn icp_pulls_a_displaced_scan_back_onto_the_model() {
    let model = TriangleSet::new(&box_mesh(40.0, 30.0, 20.0));
    let pts = surface_points(40.0, 30.0, 20.0);
    let mut off = rotation(Vec3::new(0.2, 0.5, 1.0), 0.05);
    off.t = Vec3::new(0.8, -0.6, 0.4);
    let scan: Vec<Vec3> = pts.iter().map(|p| off.apply(*p)).collect();
    let found = icp(&scan, &model, Rigid::IDENTITY, &IcpOptions { iterations: 100, ..IcpOptions::default() }).expect("icp");
    assert!(found.rms < 1e-3, "rms {}", found.rms);
    assert!(scan.iter().all(|p| model.closest_point(found.rigid.apply(*p), f64::INFINITY).is_some_and(|c| c.0 < 1e-2)));
}

#[test]
fn deviation_is_signed_by_the_surface_side() {
    let model = TriangleSet::new(&box_mesh(10.0, 10.0, 10.0));
    let out = deviation(Vec3::new(5.0, 5.0, 12.5), &model).expect("outside");
    let inside = deviation(Vec3::new(5.0, 5.0, 8.0), &model).expect("inside");
    let on = deviation(Vec3::new(5.0, 5.0, 10.0), &model).expect("on");
    assert!((out.signed - 2.5).abs() < 1e-9, "{out:?}");
    assert!((inside.signed + 2.0).abs() < 1e-9, "{inside:?}");
    assert_eq!(on.signed, 0.0);
}

#[test]
fn stats_count_each_value_once() {
    let values = [-0.3, -0.1, 0.0, 0.1, 0.2, 0.9];
    let s = deviation_stats(&values, 0.15, 0.5, 5).expect("stats");
    assert_eq!(s.count, 6);
    assert_eq!(s.histogram.iter().sum::<usize>(), 6);
    assert_eq!(s.histogram.first(), Some(&0));
    assert_eq!(s.histogram.last(), Some(&1));
    assert!((s.within_tolerance - 3.0 / 6.0).abs() < 1e-12);
    assert!((s.min + 0.3).abs() < 1e-12 && (s.max - 0.9).abs() < 1e-12);
}

#[test]
fn hostile_input_gives_none_never_a_panic() {
    let nan = Vec3::new(f64::NAN, 0.0, 0.0);
    let set = TriangleSet::new(&box_mesh(1.0, 1.0, 1.0));
    assert!(principal_frame(&[nan, Vec3::X, Vec3::Y]).is_none());
    assert!(principal_frame(&[Vec3::ZERO, Vec3::X]).is_none());
    // Points on a line define no frame.
    assert!(principal_frame(&[Vec3::ZERO, Vec3::X, Vec3::X * 2.0, Vec3::X * 3.0]).is_none());
    assert!(fit_rigid(&[Vec3::ZERO], &[Vec3::ZERO]).is_none());
    assert!(fit_rigid(&[Vec3::ZERO, Vec3::X, Vec3::Y], &[nan, Vec3::X, Vec3::Y]).is_none());
    assert!(icp(&[Vec3::ZERO; 2], &set, Rigid::IDENTITY, &IcpOptions::default()).is_none());
    assert!(icp(&[Vec3::ZERO, Vec3::X, Vec3::Y], &TriangleSet::default(), Rigid::IDENTITY, &IcpOptions::default()).is_none());
    assert!(deviation(nan, &set).is_none());
    assert!(deviation_stats(&[f64::NAN], 1.0, 1.0, 4).is_none());
    assert!(deviation_stats(&[0.0], 1.0, 0.0, 4).is_none());
    assert!(deviation_stats(&[0.0], 1.0, 1.0, 0).is_none());
    assert_eq!(heat_color(f64::NAN, 1.0), [128, 128, 128]);
}

#[test]
fn half_turn_has_an_axis() {
    let m = rotation(Vec3::new(0.0, 1.0, 0.0), std::f64::consts::PI);
    let (axis, angle) = m.axis_angle();
    assert!((angle - std::f64::consts::PI).abs() < 1e-9);
    assert!(axis.dot(Vec3::Y).abs() > 1.0 - 1e-9);
}
