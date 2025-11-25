#pragma once

#include "mesh_processing.h"

#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/mesh.h>

#include <filesystem>

namespace FixMeshNormals {

// Reads mesh data from a USD prim into our processing format
bool ReadMeshFromUsd(const pxr::UsdGeomMesh& mesh, MeshData& outData);

// Writes processed mesh to USD as a new prim with "_voxelised" suffix
bool WriteMeshToUsd(pxr::UsdStageRefPtr stage, const pxr::SdfPath& originalPath, const ProcessedMesh& mesh, bool overwrite);

// Processes all meshes in a USD stage file
std::size_t ProcessStageFile(const std::filesystem::path& inputPath,
							 const std::filesystem::path& outputPath,
							 bool overwrite,
							 bool processElements,
							 bool writeVoxelised,
							 double weldTolerance);

} // namespace FixMeshNormals
