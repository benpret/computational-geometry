# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FixMeshNormals is a C++ command-line utility for processing USD (Universal Scene Description) mesh files. It voxelizes meshes using OpenVDB to generate consistent, outward-facing normals for rendering. The tool converts USD meshes to signed distance fields and back, creating new "_voxelised" variants with correct face normals.

**Key Dependencies:**
- Pixar USD (Universal Scene Description)
- OpenVDB (sparse volumetric data structures)
- Intel TBB (Threading Building Blocks)
- Boost (via USD)
- MaterialX

## Build System

This is a **Visual Studio 2022** Windows project using MSBuild (`.vcxproj`).

### Environment Variables Required

The project relies on environment variables for dependency paths:

**Debug builds:**
- `USD_DEBUG_ROOT` - Path to USD debug libraries and headers
- `PYTHON_INSTALL_ROOT` - Python installation path

**Release builds:**
- `USD_RELEASE_ROOT` - Path to USD release libraries and headers
- `PYTHON_INSTALL_ROOT` - Python installation path

### Build Commands

```bash
# Build Debug x64 (recommended for development)
msbuild FixMeshNormals.sln /p:Configuration=Debug /p:Platform=x64

# Build Release x64
msbuild FixMeshNormals.sln /p:Configuration=Release /p:Platform=x64

# Clean
msbuild FixMeshNormals.sln /t:Clean /p:Configuration=Debug /p:Platform=x64
```

**Note:** The project is configured for x64 architecture only (Win32 configs exist but x64 is primary target).

## Running the Tool

```bash
# Process USD file in-place (modifies original)
./x64/Debug/FixMeshNormals.exe path/to/mesh.usda

# Process with explicit output
./x64/Debug/FixMeshNormals.exe --out output.usda input.usda

# Process multiple files (each written back in-place)
./x64/Debug/FixMeshNormals.exe file1.usda file2.usda file3.usda

# Show help
./x64/Debug/FixMeshNormals.exe --help
```

## Code Architecture

The codebase follows a **modular, layered architecture** with clear separation of concerns:

### Module Structure

**cli.h/cpp** - Command-line interface (Pure)
- `Options` struct - holds parsed arguments
- `ParseArgs()` - parses command-line arguments
- `PrintUsage()` - displays help text
- No dependencies on USD or OpenVDB

**mesh_processing.h/cpp** - Core algorithms (Pure Core)
- `MeshData` struct - intermediate representation for voxelization
- `ProcessedMesh` struct - voxelized output with normals
- `VoxelizeMesh()` - main processing pipeline (mesh → volume → mesh)
- `ComputeFaceNormals()` - geometric normal calculation
- `ComputeExtentFromPoints()` - bounding box computation
- Only depends on Pixar math types and OpenVDB

**usd_io.h/cpp** - USD file operations (Impure Edge)
- `ReadMeshFromUsd()` - extracts mesh data from USD prim
- `WriteMeshToUsd()` - writes processed mesh to USD
- `ProcessStageFile()` - orchestrates stage-level processing
- All USD I/O isolated here

**main.cpp** - Thin orchestration layer
- Initializes OpenVDB
- Validates inputs
- Loops over files and calls `ProcessStageFile()`
- Only 60 lines

### Processing Pipeline

1. **USD Mesh Reading** (`ReadMeshFromUsd` in usd_io.cpp) - Extracts points and topology from USD mesh prims
2. **Voxelization** (`VoxelizeMesh` in mesh_processing.cpp) - Converts mesh to OpenVDB signed distance field and back
3. **Normal Computation** (`ComputeFaceNormals` in mesh_processing.cpp) - Computes face normals from reconstructed geometry
4. **USD Writing** (`WriteMeshToUsd` in usd_io.cpp) - Creates new `_voxelised` mesh prim alongside original

### Important Implementation Details

**Mesh Handling:**
- Skips existing `_voxelised` meshes to allow idempotent re-runs (usd_io.cpp)
- Handles both `GfVec3f` and `GfVec3d` point types from USD
- Fan-triangulates n-gons for OpenVDB compatibility (usd_io.cpp)

**Voxelization Parameters:**
- Auto-computes voxel size as `max_dimension / 64` (mesh_processing.cpp)
- Uses 3-voxel exterior band and 1-voxel interior band
- Isovalue of 0.0 for surface extraction

**Normal Generation:**
- Computes face normals (uniform interpolation), not vertex normals
- Uses cross product of first triangle edges per face
- Zero-length safeguard for degenerate faces (mesh_processing.cpp)

**Path Handling:**
- Converts relative paths to absolute (Windows WSL compatibility)
- Uses `generic_string()` for forward slashes (USD prefers this on Windows) (usd_io.cpp)

## Design Principles

This project follows the **Agentic-Agnostic Software Design Guide** in `resources/Agnostic-Software-Design-Guide.md`. Key principles applied:

- **Single Responsibility:** Each function has one clear job (parse mesh, voxelize, compute normals, etc.)
- **Pure Core, Impure Edge:** Mesh processing logic is isolated; USD I/O only at boundaries
- **Fail Fast:** Validates mesh data and skips invalid geometry with clear error messages
- **Minimal Logging:** Only actionable messages (found mesh, created output, errors)

When extending this code:
- Keep processing functions pure (take data, return data)
- Add validation early in the pipeline
- Use self-documenting names over comments
- Maintain separation between USD/VDB conversions

## Common Modification Patterns

**Adjusting voxel resolution:**
Modify `defaultDivisions` in `ComputeVoxelSize()` (mesh_processing.cpp)

**Changing normal interpolation:**
Modify `SetNormalsInterpolation()` call in `WriteMeshToUsd()` (usd_io.cpp)
Options: `"uniform"` (face), `"vertex"`, `"faceVarying"`

**Adding mesh attributes:**
- Extend `ProcessedMesh` struct in mesh_processing.h
- Compute new attributes in `VoxelizeMesh()` (mesh_processing.cpp)
- Write attributes in `WriteMeshToUsd()` (usd_io.cpp)

**Supporting additional command-line options:**
Extend `Options` struct and `ParseArgs()` in cli.h/cpp

**Changing processing algorithm:**
Modify `VoxelizeMesh()` in mesh_processing.cpp (core logic is isolated here)
