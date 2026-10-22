# Reverse engineering from scans — handoff notes (fork karolkrzyzanowski/solidcraft)

Owner: Karol (CELUST / BANDIT OEM; Bambu Lab X1C, 3D scanner). Reply to him in **Polish**, short and ready to use.
Any agent (Claude, ChatGPT/Codex, other) continuing this work must first read `AGENTS.md` (no panics/unwrap,
Result errors, clean-room, layering, `cargo xtask ci`, tests in existing files, new labels in every complete catalog).

## Goal
QuickSurface-style tools (mesh → CAD) as SolidCraft **commands** (`mesh.*`), usable from GUI, CLI, MCP.
Order decided with Karol: 1) alignment + deviation (this branch) → 2) primitive detection (RANSAC planes/cylinders/
spheres/cones/tori) → 3) primitives → solid → 4) freeform surfaces. Also wanted later: Bambu-style 3MF (needs a
reference .3mf from Bambu Studio for X1C).

## Done (branches on the fork)
- `pl-translation-completion` (a11803c): complete Polish UI. Running on Karol's Mac as "SolidCraft PL"
  (~/projects/solidcraft, build with `export PATH="$HOME/.cargo/bin:$PATH"; cargo build --release -p solidcraft`).
- `mesh-align-deviation` (WIP, **not yet verified**):
  - `crates/geom/src/registration.rs` (+ `registration_tests.rs`, 7 tests pass): principal frame, `align_to_axes`,
    Horn `fit_rigid`, trimmed `icp`, signed `deviation`, `deviation_stats`, `heat_color`.
  - `TriangleSet::closest_point` and `normal` in `crates/geom/src/mesh_distance.rs`.
  - `crates/io/src/ply.rs`: `write_ply_colored` (binary PLY, per-vertex colour).
  - `crates/engine/src/cmd/mesh_fit.rs` (+ `mesh_fit_tests.rs`): `mesh.align` (axes mode, or fit to reference body
    by ICP from best of 5 starts; emits a normal Move feature, so undoable) and `mesh.deviation` (stats, histogram,
    worst points, optional `save: x.ply` colour map; read-only, echoes verdict).
  - Catalog keys `cmd.mesh_align`, `cmd.mesh_deviation` in en/pl/ja/zh/de/ko; both listed in SURFACE → CONVERT
    (`crates/ui-egui/src/workspace.rs`). Without a dialog the toolbar click runs the command with the selection
    (mesh body = scan, other selected body = reference).

## Assumptions
- `inspect.compare` already compares two bodies (face colouring); it does not give per-vertex scan deviation, hence `mesh.deviation`.
- Move rotates about the world origin then translates, so `p' = R p + t` maps directly to Move params.
- Positive deviation = scan outside the reference.
- Colour map is only in the PLY file; no in-app heat map yet (needs vertex-colour support in render/gpu/ui).
- No real scan from Karol yet; tests use a synthetic bracket.

## Not done / next
1. Get `cargo test -p solidcraft-engine --lib mesh_fit` green (was compiling; slow in the sandbox, build on the Mac).
2. Run `cargo fmt`, `cargo clippy -D warnings`, `cargo xtask i18n`, `cargo xtask layers`, `cargo xtask wasm`, `cargo xtask ci`.
3. Check the toolbar/command-registry tests (icons, docs/features are generated on main — do not regenerate).
4. Rebuild on the Mac, try with a real STL; consider a results window for the numbers and an in-app colour map.
5. Then module 2 (primitive detection), new commands in a new file under `crates/engine/src/cmd/`.

## Backlog (added 2026-10-22)
- `mesh.align` / `mesh.deviation`: activate the scan's component automatically (today: "activate the mesh body's component before aligning it" after file.insert_mesh; same rule as mesh.convert).
- Colour map of deviations inside the app (vertex colours in render/gpu/ui) and a results window for the numbers; today only PLY export + status line.
- Done since: `.ply` import (`read_ply`), checked on an 8.7M-triangle Revopoint scan.
- Surface-from-four-curves panel tool (idea from the Blender add-on "Surface Mesh" by 3Dartmodeling, concept-car panels from four Bézier boundary curves, live update, shared corners, mirror, flip normals). **Licence: that add-on is GPL and sold on Superhive: never read or copy its source.** Build from public documentation of what it does and the textbook Coons patch (S. A. Coons, 1967); SolidCraft already has the maths in `crates/geom/src/network.rs` and commands `surface.boundary`, `surface.patch`, `surface.match`, `surface.blend`, `surface.fair`. Missing: a one-step four-curve panel command with auto-fixed curve directions, shared-corner curve helpers, mirrored panel. Lower priority than primitive detection.
- Engine error strings are English (not translated); Dock label shows "solidcraft".
- Unchecked: the 3 failing upstream tests on the Mac (prism_local_tests, robustness_goldens seed_36/seed_71): confirm they fail on a clean branch too.
- Upstream is read-only; push only to the fork. Mac device calls time out at 60 s: run long jobs with nohup and poll.
