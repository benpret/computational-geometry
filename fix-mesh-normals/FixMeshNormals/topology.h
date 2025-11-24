#pragma once

#include <pxr/base/gf/vec3d.h>
#include <pxr/base/gf/vec3f.h>

#include <map>
#include <vector>

namespace FixMeshNormals {

using Face = std::vector<int>;
using EdgeKey = std::pair<int, int>;

// Edge -> [(faceIdx, v0, v1), ...]
using EdgeMap = std::map<EdgeKey, std::vector<std::tuple<int, int, int>>>;

// Build edge-to-face adjacency map
EdgeMap BuildEdgeMap(const std::vector<Face>& faces);

// Unify face winding using half-edge BFS
// Returns: (orientation per face, connected components, conflict count)
struct WindingResult {
	std::vector<int> orientation;        // +1 or -1 per face
	std::vector<std::vector<int>> components;  // connected components (face indices)
	int conflicts;                       // edges with conflicting winding
};

WindingResult UnifyWinding(const std::vector<Face>& faces);

// Apply orientation to faces (-1 = reverse, +1 = keep)
std::vector<Face> ApplyOrientation(const std::vector<Face>& faces, const std::vector<int>& orientation);

// Compute polygon area vector (for normal and volume calculations)
pxr::GfVec3d PolygonAreaVector(const std::vector<pxr::GfVec3f>& points, const Face& face);

// Compute signed volume of a set of faces
double SignedVolume(const std::vector<pxr::GfVec3f>& points,
					const std::vector<Face>& faces,
					const std::vector<int>& faceIndices);

// Compute centroid flux (for open meshes)
double CentroidFlux(const std::vector<pxr::GfVec3f>& points,
					const std::vector<Face>& faces,
					const std::vector<int>& faceIndices);

// Count faces pointing outward vs inward from component centroid
std::pair<int, int> FaceRadialSigns(const std::vector<pxr::GfVec3f>& points,
									 const std::vector<Face>& faces,
									 const std::vector<int>& faceIndices);

// Compute component centroid
pxr::GfVec3d ComponentCentroid(const std::vector<pxr::GfVec3f>& points,
							   const std::vector<Face>& faces,
							   const std::vector<int>& faceIndices);

} // namespace FixMeshNormals
