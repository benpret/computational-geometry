#pragma once

#include <pxr/base/gf/range3d.h>
#include <pxr/base/gf/vec3f.h>
#include <pxr/base/vt/array.h>

#include <openvdb/openvdb.h>

#include <vector>

namespace FixMeshNormals {

// Mesh representation for voxelization
struct MeshData {
	std::vector<openvdb::Vec3s> points;
	std::vector<openvdb::Vec3I> triangles;
	std::vector<openvdb::Vec4I> quads;
	pxr::GfRange3d bounds;

	bool IsEmpty() const {
		return points.empty() || (triangles.empty() && quads.empty());
	}
};

// Output mesh from voxelization
struct ProcessedMesh {
	pxr::VtArray<pxr::GfVec3f> points;
	pxr::VtArray<int> faceVertexCounts;
	pxr::VtArray<int> faceVertexIndices;
	pxr::VtArray<pxr::GfVec3f> normals;
	pxr::VtArray<pxr::GfVec3f> extent;
};

// Fix face winding and orientation before voxelization
// Returns: corrected (counts, indices)
std::pair<pxr::VtArray<int>, pxr::VtArray<int>> FixWindingAndOrientation(
	const pxr::VtArray<pxr::GfVec3f>& points,
	const pxr::VtArray<int>& counts,
	const pxr::VtArray<int>& indices);

// Use VDB signed distance field as oracle to fix face orientations
// Returns: (corrected counts, corrected indices, flip count)
std::tuple<pxr::VtArray<int>, pxr::VtArray<int>, int> FixOrientationWithVDBOracle(
	const pxr::VtArray<pxr::GfVec3f>& points,
	const pxr::VtArray<int>& counts,
	const pxr::VtArray<int>& indices,
	openvdb::FloatGrid::Ptr sdfGrid);

// Core voxelization: converts mesh to volume and back with corrected normals
ProcessedMesh VoxelizeMesh(const MeshData& input);

// Create VDB signed distance field from mesh (exposed for oracle usage)
openvdb::FloatGrid::Ptr CreateSDFFromMesh(const MeshData& input);

// Utilities
pxr::VtArray<pxr::GfVec3f> ComputeFaceNormals(const pxr::VtArray<pxr::GfVec3f>& points,
											  const pxr::VtArray<int>& counts,
											  const pxr::VtArray<int>& indices);

pxr::VtArray<pxr::GfVec3f> ComputeExtentFromPoints(const pxr::VtArray<pxr::GfVec3f>& points);

} // namespace FixMeshNormals
