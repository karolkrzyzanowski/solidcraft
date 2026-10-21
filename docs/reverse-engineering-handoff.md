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
- Upstream is read-only; push only to the fork. Mac device calls time out at 60 s: run long jobs with nohup and poll.
