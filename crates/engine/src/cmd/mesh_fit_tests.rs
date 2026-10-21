use serde_json::{Value, json};
use solidcraft_doc::{Document, FeatureKind, MeshData};

use crate::Session;

fn run(s: &mut Session, id: &str, p: Value) -> Value {
    match s.execute(id, &p) {
        Ok(v) => v,
        Err(e) => panic!("{id} {p}: {e}"),
    }
}

/// An L-shaped bracket with a bore (a solid named `Bracket`) and, beside it in the design, a
/// mesh body `Scan` of the same part that has been turned and moved away.
fn scanned_bracket() -> Session {
    let mut s = Session::default();
    run(&mut s, "sketch.create", json!({"plane":"XY","name":"Section"}));
    run(&mut s, "sketch.line", json!({"points":[[0,0],[40,0],[40,20],[20,20],[20,30],[0,30]],"closed":true}));
    run(&mut s, "sketch.circle.center", json!({"center":[10,10],"radius":4}));
    run(&mut s, "sketch.finish", json!({}));
    run(&mut s, "solid.extrude", json!({"sketch":"Section","profiles":[{"point":[30,10]}],"distance":15,"body_name":"Bracket"}));
    let mesh = s.model.state().body("Bracket").unwrap().body.tessellate(0.01).unwrap();
    s.doc_mut()
        .add_feature(
            FeatureKind::MeshImport {
                file: "public-safe-bracket.stl".into(),
                meshes: vec![MeshData {
                    name: "Scan".into(),
                    positions: mesh.positions.iter().flat_map(|p| [p.x, p.y, p.z]).collect(),
                    triangles: mesh.triangles.iter().flatten().copied().collect(),
                    color: None,
                }],
            },
            Some("Imported scan"),
        )
        .unwrap();
    s.model = Session::new(s.doc.clone()).model;
    run(&mut s, "solid.move", json!({"bodies":["Scan"],"translate":[30,-12,7],"axis":[1,2,3],"angle":"25 deg"}));
    s
}

#[test]
fn a_displaced_scan_is_fitted_onto_its_model_and_measured() {
    let mut s = scanned_bracket();
    let before = run(&mut s, "mesh.deviation", json!({"body":"Scan","to":"Bracket"}));
    assert!(before["rms_mm"].as_f64().unwrap() > 5.0, "{before}");
    let fit = run(&mut s, "mesh.align", json!({"body":"Scan","to":"Bracket","name":"Fit scan"}));
    assert_eq!(fit["mode"], "fit");
    assert!(fit["rms_mm"].as_f64().unwrap() < 0.05, "{fit}");
    assert!(s.model.results.iter().all(|r| r.error.is_none()), "{:?}", s.model.results);
    let after = run(&mut s, "mesh.deviation", json!({"body":"Scan","to":"Bracket","tolerance":0.1}));
    assert!(after["rms_mm"].as_f64().unwrap() < 0.05, "{after}");
    assert!(after["within_tolerance"].as_f64().unwrap() > 0.99, "{after}");
    let bins: Vec<u64> = after["histogram"]["bins"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap()).collect();
    assert_eq!(bins.iter().sum::<u64>(), after["count"].as_u64().unwrap());
    // The fit is a timeline step: undo puts the scan back where it was.
    run(&mut s, "edit.undo", json!({}));
    let undone = run(&mut s, "mesh.deviation", json!({"body":"Scan","to":"Bracket"}));
    assert!((undone["rms_mm"].as_f64().unwrap() - before["rms_mm"].as_f64().unwrap()).abs() < 1e-9);
}

#[test]
fn aligning_to_the_axes_is_stable() {
    let mut s = scanned_bracket();
    let first = run(&mut s, "mesh.align", json!({"body":"Scan","center":[0,0,0]}));
    assert_eq!(first["mode"], "axes");
    let second = run(&mut s, "mesh.align", json!({"body":"Scan"}));
    assert!(second["rotation"]["degrees"].as_f64().unwrap() < 1e-4, "{second}");
    assert!(second["translation"].as_array().unwrap().iter().all(|x| x.as_f64().unwrap().abs() < 1e-6), "{second}");
}

#[test]
fn the_selected_bodies_are_the_scan_and_the_reference() {
    let mut s = scanned_bracket();
    s.selection = vec![crate::Sel::Body { name: "Bracket".into() }, crate::Sel::Body { name: "Scan".into() }];
    let fit = run(&mut s, "mesh.align", json!({}));
    assert_eq!(fit["reference"], "Bracket");
    assert!(fit["rms_mm"].as_f64().unwrap() < 0.05, "{fit}");
    let dev = run(&mut s, "mesh.deviation", json!({}));
    assert_eq!(dev["scan"], "Scan");
}

#[test]
fn the_colour_map_is_saved_as_a_ply_with_a_colour_for_every_vertex() {
    let mut s = scanned_bracket();
    run(&mut s, "mesh.align", json!({"body":"Scan","to":"Bracket"}));
    let path = std::env::temp_dir().join(format!("solidcraft-deviation-{}.ply", std::process::id()));
    let r = run(&mut s, "mesh.deviation", json!({"body":"Scan","to":"Bracket","save":path,"range":0.5}));
    let bytes = std::fs::read(&path).unwrap();
    assert_eq!(r["ply"]["bytes"].as_u64().unwrap() as usize, bytes.len());
    let head = String::from_utf8_lossy(&bytes[..bytes.len().min(400)]).into_owned();
    assert!(head.starts_with("ply\nformat binary_little_endian 1.0\n") && head.contains(&format!("element vertex {}", r["count"])), "{head}");
    std::fs::remove_file(path).unwrap();
}

#[test]
fn bad_parameters_are_refused_and_change_nothing() {
    let mut s = scanned_bracket();
    let before = s.doc.to_json();
    for (id, p) in [
        ("mesh.align", json!({})),
        ("mesh.align", json!({"body":"Missing"})),
        ("mesh.align", json!({"body":"Bracket"})),
        ("mesh.align", json!({"body":"Scan","to":"Scan"})),
        ("mesh.align", json!({"body":"Scan","to":"Nothing"})),
        ("mesh.align", json!({"body":"Scan","to":"Bracket","iterations":0})),
        ("mesh.align", json!({"body":"Scan","to":"Bracket","iterations":1e9})),
        ("mesh.align", json!({"body":"Scan","to":"Bracket","keep":2})),
        ("mesh.align", json!({"body":"Scan","to":"Bracket","keep":"all"})),
        ("mesh.align", json!({"body":"Scan","to":"Bracket","initial":"random"})),
        ("mesh.align", json!({"body":"Scan","to":"Bracket","tolerance":0})),
        ("mesh.align", json!({"body":"Scan","center":[1e308,0,0]})),
        ("mesh.align", json!({"body":"Scan","center":[0,0]})),
        ("mesh.deviation", json!({})),
        ("mesh.deviation", json!({"body":"Scan"})),
        ("mesh.deviation", json!({"body":"Scan","to":"Bracket","bins":1})),
        ("mesh.deviation", json!({"body":"Scan","to":"Bracket","range":0})),
        ("mesh.deviation", json!({"body":"Scan","to":"Bracket","tolerance":-1})),
        ("mesh.deviation", json!({"body":"Scan","to":"Bracket","max_points":1e12})),
        ("mesh.deviation", json!({"body":"Scan","to":"Bracket","save":"scan.txt"})),
        ("mesh.deviation", json!({"body":"Scan","to":"Bracket","save":7})),
        ("mesh.deviation", json!({"body":"Scan","to":"Bracket","save":"/tmp/solidcraft-too-small.ply","max_points":100})),
    ] {
        assert!(s.execute(id, &p).is_err(), "{id} {p} should be refused");
        assert_eq!(s.doc.to_json(), before, "{id} {p}");
    }
}

#[test]
fn an_empty_design_has_nothing_to_align() {
    let mut s = Session::new(Document::default());
    assert!(s.execute("mesh.align", &json!({"body":"Scan"})).is_err());
    assert!(s.execute("mesh.deviation", &json!({"body":"Scan","to":"Other"})).is_err());
}
