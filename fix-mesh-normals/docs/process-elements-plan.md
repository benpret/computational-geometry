# Plan: Per-Element Mesh Processing

Goal: allow processing either whole `UsdGeomMesh` prims (current behavior) or each connected element within a mesh when `--process-elements` is set. An "element" is a connected set of faces. Each element is processed independently through the 11-step pipeline, then elements are recombined into a single `<original_name>_fixed` mesh (and a single `_voxelised` mesh if requested) placed back at the original world location.

## Target Behavior
- Default: process each `UsdGeomMesh` as today, but run the new pipeline (recenter → weld → fix winding/orientation → SDF oracle → voxelise → restore transform) before writing `_fixed` and optionally `_voxelised`.
- `--process-elements`: split each mesh into connected components of faces. Process each component independently through the same pipeline, then reassemble into a single combined `_fixed` (and optional `_voxelised`) mesh. Avoid touching prims that already end with `_fixed` or `_voxelised`.

## Pipeline per Item
1) Center pivot to the item: compute pivot (component bounds center) in object space.  
2) Store original location: record pivot offset so we can restore later.  
3) Move item to world origin: subtract pivot from all points for processing.  
4) Weld vertices: merge coincident points using a configurable tolerance (default 1 mm) and remap indices.  
5) FixWindingAndOrientation: run on welded geometry to get consistent winding and component orientation.  
6) Convert to `MeshData` for SDF creation.  
7) Convert topology to triangles/quads for VDB (fan-triangulate n-gons).  
8) Create SDF from winding-fixed mesh (oracle).  
9) Use SDF oracle to fix component face orientations.  
10) Voxelise the component (always computed).  
11) Restore location: add pivot back and merge components into `_fixed`; if voxel output is enabled, also merge into `_voxelised`, then write prims.

## Implementation Steps
- CLI plumbing (`cli.h/.cpp`, `main.cpp`):
  - Add `Options::processElements` with `--process-elements` flag; update usage text.
  - Add `--write-voxelised` flag to control writing `_voxelised` output (voxelisation still runs internally for orientation oracle). Default: on.
  - Thread the flags into `ProcessStageFile` signature and call sites.
- Mesh/element partitioning (`usd_io.cpp`, maybe `topology.*`):
  - Build a helper to split a mesh into connected components using face adjacency (reuse `UnifyWinding` components or a dedicated BFS). Each component should expose sublists of counts/indices and the referenced point indices.
  - For default mode, produce a single component covering the entire mesh so downstream code is shared.
- Pivot + translation handling:
  - For each component, compute bounds of its referenced points; pivot = bounds center (fallback to geometric centroid if bounds empty).
  - Subtract pivot from a working copy of the component’s points before welding; keep the pivot vector to reapply later.
- Vertex welding (new helper in `mesh_processing` or `usd_io`):
  - Merge points within a tolerance (default 1 mm; also accept user override).
  - Return welded points plus remapped face indices; drop degenerate faces created by welding.
- Winding/orientation stage:
  - Run `FixWindingAndOrientation` on the welded component; keep the component membership info to avoid re-splitting.
  - Build `MeshData` from the winding-fixed topology (triangulate n-gons) and capture updated bounds in the translated space.
- Oracle correction + voxelisation:
  - Create the SDF grid from `MeshData` (step 8) and call `FixOrientationWithVDBOracle` on the same component.
  - Generate `_fixed` geometry from oracle-corrected topology, translate points back by pivot, recompute extents, and merge all components into one `_fixed` mesh.
  - Feed the winding-fixed `MeshData` through `VoxelizeMesh` to get the VDB rebuild; translate result points back by pivot, and if voxel output is enabled, merge all components into one `_voxelised` mesh and write it.
- Output naming and placement:
  - Keep existing suffixes. Always emit a single `<mesh>_fixed`; optionally emit a single `<mesh>_voxelised` when requested. Preserve current skip/overwrite logic for existing prims.
- Validation:
  - Add small USD fixtures (multi-component mesh) to verify per-element processing and translation round-trips.
  - Exercise both modes: whole-mesh and `--process-elements`, with voxel output on/off.

## Open Decisions to Settle Before Coding
- None — decisions recorded above (use `--write-voxelised` to control `_voxelised` emission).
