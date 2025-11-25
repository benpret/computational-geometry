#include "mesh_processing.h"
#include "topology.h"

#include <pxr/base/gf/vec3d.h>

#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace FixMeshNormals {

namespace {

double ComputeVoxelSize(const pxr::GfRange3d& bounds) {
	if (bounds.IsEmpty()) {
		return 0.01;
	}

	const pxr::GfVec3d size = bounds.GetSize();
	const double maxDim = std::max({size[0], size[1], size[2]});
	const double defaultDivisions = 512;
	double voxelSize = maxDim > 0.0 ? maxDim / defaultDivisions : 0.01;

	const double minVoxelSize = 1e-4;
	if (voxelSize < minVoxelSize) {
		voxelSize = minVoxelSize;
	}

	return voxelSize;
}

openvdb::FloatGrid::Ptr MeshToGrid(const MeshData& mesh) {
	if (mesh.IsEmpty()) {
		return nullptr;
	}

	const double voxelSize = ComputeVoxelSize(mesh.bounds);
	openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxelSize);

	// Band widths in world units
	const double exteriorBandWidth = 3.0 * voxelSize;
	// Use max double for interior band - tells OpenVDB to fill entire interior
	const double interiorBandWidth = std::numeric_limits<double>::max();

	return openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(
		*transform, mesh.points, mesh.triangles, mesh.quads, exteriorBandWidth, interiorBandWidth);
}

bool GridToMesh(const openvdb::FloatGrid& grid,
				pxr::VtArray<pxr::GfVec3f>& outPoints,
				pxr::VtArray<int>& outCounts,
				pxr::VtArray<int>& outIndices) {
	std::vector<openvdb::Vec3s> vdbPoints;
	std::vector<openvdb::Vec3I> triPolys;
	std::vector<openvdb::Vec4I> quadPolys;

	openvdb::tools::volumeToMesh(grid, vdbPoints, triPolys, quadPolys, /*isovalue=*/0.0, /*adaptivity=*/0.0);

	if (vdbPoints.empty() || (triPolys.empty() && quadPolys.empty())) {
		return false;
	}

	outPoints.reserve(vdbPoints.size());
	for (const auto& p : vdbPoints) {
		outPoints.emplace_back(p[0], p[1], p[2]);
	}

	outCounts.reserve(triPolys.size() + quadPolys.size());
	outIndices.reserve(triPolys.size() * 3 + quadPolys.size() * 4);

	// Reverse winding order to flip normals outward
	// OpenVDB's volumeToMesh produces inward-facing normals by default
	for (const auto& tri : triPolys) {
		outCounts.push_back(3);
		outIndices.push_back(tri[0]);
		outIndices.push_back(tri[2]);  // Reversed
		outIndices.push_back(tri[1]);  // Reversed
	}

	for (const auto& quad : quadPolys) {
		outCounts.push_back(4);
		outIndices.push_back(quad[0]);
		outIndices.push_back(quad[3]);  // Reversed
		outIndices.push_back(quad[2]);  // Reversed
		outIndices.push_back(quad[1]);  // Reversed
	}

	return true;
}

std::vector<Face> UnpackFaces(const pxr::VtArray<int>& counts, const pxr::VtArray<int>& indices) {
	std::vector<Face> faces;
	std::size_t cursor = 0;
	for (int count : counts) {
		Face face;
		for (int i = 0; i < count; ++i) {
			face.push_back(indices[cursor + i]);
		}
		faces.push_back(face);
		cursor += count;
	}
	return faces;
}

std::pair<pxr::VtArray<int>, pxr::VtArray<int>> FlattenFaces(const std::vector<Face>& faces) {
	pxr::VtArray<int> counts;
	pxr::VtArray<int> flatIndices;
	for (const Face& face : faces) {
		counts.push_back(static_cast<int>(face.size()));
		for (int v : face) {
			flatIndices.push_back(v);
		}
	}
	return {counts, flatIndices};
}

} // anonymous namespace

std::pair<pxr::VtArray<int>, pxr::VtArray<int>> FixWindingAndOrientation(
	const pxr::VtArray<pxr::GfVec3f>& points,
	const pxr::VtArray<int>& counts,
	const pxr::VtArray<int>& indices) {

	// Convert to Face representation
	std::vector<Face> faces = UnpackFaces(counts, indices);
	std::vector<pxr::GfVec3f> pointsVec(points.begin(), points.end());

	// Step 1: Unify winding
	WindingResult windingResult = UnifyWinding(faces);
	std::vector<Face> orientedFaces = ApplyOrientation(faces, windingResult.orientation);

	// Step 2: Fix component orientation
	std::vector<int> finalOrientation = windingResult.orientation;

	// Compute bounding box for tolerance
	pxr::GfRange3d bbox;
	for (const auto& p : pointsVec) {
		bbox.UnionWith(pxr::GfVec3d(p));
	}
	pxr::GfVec3d bboxSize = bbox.GetSize();
	double bboxDiag = std::sqrt(bboxSize[0] * bboxSize[0] + bboxSize[1] * bboxSize[1] + bboxSize[2] * bboxSize[2]);
	if (bboxDiag <= 0.0) {
		bboxDiag = 1.0;
	}
	double volumeTolerance = 1e-12 * bboxDiag * bboxDiag * bboxDiag;

	int componentFlips = 0;
	for (const auto& comp : windingResult.components) {
		pxr::GfVec3d compCentroid = ComponentCentroid(pointsVec, orientedFaces, comp);
		double vol = SignedVolume(pointsVec, orientedFaces, comp);
		double flux = CentroidFlux(pointsVec, orientedFaces, comp);
		auto [pos, neg] = FaceRadialSigns(pointsVec, orientedFaces, comp);

		bool flipNeeded = false;

		// Decision hierarchy (matches Python logic)
		if (std::abs(vol) > volumeTolerance) {
			flipNeeded = vol < 0;
		} else if (flux != 0) {
			flipNeeded = flux < 0;
		} else if (pos != neg) {
			flipNeeded = neg > pos;
		} else {
			// Tie-break: bias away from global centroid
			pxr::GfVec3d globalCentroid(0.0);
			for (const auto& p : pointsVec) {
				globalCentroid += pxr::GfVec3d(p);
			}
			globalCentroid /= static_cast<double>(pointsVec.size());

			double score = 0.0;
			for (int faceIdx : comp) {
				const Face& face = orientedFaces[faceIdx];
				if (face.size() < 3) {
					continue;
				}
				pxr::GfVec3d areaVec = PolygonAreaVector(pointsVec, face);
				pxr::GfVec3d faceCentroid(0.0);
				for (int v : face) {
					faceCentroid += pxr::GfVec3d(pointsVec[v]);
				}
				faceCentroid /= static_cast<double>(face.size());
				score += pxr::GfDot(areaVec, faceCentroid - globalCentroid);
			}
			flipNeeded = score < 0;
		}

		if (flipNeeded) {
			++componentFlips;
			for (int faceIdx : comp) {
				std::reverse(orientedFaces[faceIdx].begin(), orientedFaces[faceIdx].end());
				finalOrientation[faceIdx] *= -1;
			}
		}

		// Step 3: Flip outlier faces within component
		std::vector<std::pair<int, double>> radialDots;
		for (int faceIdx : comp) {
			const Face& face = orientedFaces[faceIdx];
			if (face.size() < 3) {
				continue;
			}
			pxr::GfVec3d areaVec = PolygonAreaVector(pointsVec, face);
			pxr::GfVec3d faceCentroid(0.0);
			for (int v : face) {
				faceCentroid += pxr::GfVec3d(pointsVec[v]);
			}
			faceCentroid /= static_cast<double>(face.size());
			radialDots.emplace_back(faceIdx, pxr::GfDot(areaVec, faceCentroid - compCentroid));
		}

		std::vector<int> negFaces;
		int posCount = 0;
		for (const auto& [faceIdx, dot] : radialDots) {
			if (dot < -1e-9) {
				negFaces.push_back(faceIdx);
			} else if (dot > 1e-9) {
				++posCount;
			}
		}

		int maxOutliers = std::max(3, static_cast<int>(0.2 * comp.size()));
		if (!negFaces.empty() && static_cast<int>(negFaces.size()) <= maxOutliers && posCount > static_cast<int>(negFaces.size())) {
			for (int faceIdx : negFaces) {
				std::reverse(orientedFaces[faceIdx].begin(), orientedFaces[faceIdx].end());
				finalOrientation[faceIdx] *= -1;
			}
		}
	}

	// Log results
	int flippedFaces = 0;
	for (int sign : finalOrientation) {
		if (sign == -1) {
			++flippedFaces;
		}
	}
	std::cout << "    Winding fix: " << windingResult.components.size() << " components, "
			  << flippedFaces << " faces flipped, "
			  << componentFlips << " component flips, "
			  << windingResult.conflicts << " edge conflicts\n";

	return FlattenFaces(orientedFaces);
}

openvdb::FloatGrid::Ptr CreateSDFFromMesh(const MeshData& input) {
	return MeshToGrid(input);
}

std::tuple<pxr::VtArray<int>, pxr::VtArray<int>, int> FixOrientationWithVDBOracle(
	const pxr::VtArray<pxr::GfVec3f>& points,
	const pxr::VtArray<int>& counts,
	const pxr::VtArray<int>& indices,
	openvdb::FloatGrid::Ptr sdfGrid) {

	if (!sdfGrid) {
		return {counts, indices, 0};
	}

	std::vector<Face> faces = UnpackFaces(counts, indices);
	std::vector<pxr::GfVec3f> pointsVec(points.begin(), points.end());

	// Get the grid's accessor for fast sampling
	auto accessor = sdfGrid->getAccessor();
	auto& transform = sdfGrid->transform();

	int flippedCount = 0;
	const double epsilon = 5;  // Small offset along normal for sampling

	for (std::size_t faceIdx = 0; faceIdx < faces.size(); ++faceIdx) {
		const Face& face = faces[faceIdx];
		if (face.size() < 3) {
			continue;
		}

		// Compute face center
		pxr::GfVec3d faceCenter(0.0);
		for (int v : face) {
			faceCenter += pxr::GfVec3d(pointsVec[v]);
		}
		faceCenter /= static_cast<double>(face.size());

		// Compute face normal (using first 3 vertices)
		pxr::GfVec3d p0(pointsVec[face[0]]);
		pxr::GfVec3d p1(pointsVec[face[1]]);
		pxr::GfVec3d p2(pointsVec[face[2]]);
		pxr::GfVec3d normal = pxr::GfCross(p1 - p0, p2 - p0);
		double length = normal.GetLength();
		if (length < 1e-12) {
			continue;  // Degenerate face
		}
		normal /= length;

		// Sample SDF at point slightly offset along the normal
		pxr::GfVec3d samplePoint = faceCenter + epsilon * normal;
		openvdb::Vec3d worldPos(samplePoint[0], samplePoint[1], samplePoint[2]);
		openvdb::Vec3d indexPos = transform.worldToIndex(worldPos);

		// Sample the distance field
		float distance = accessor.getValue(openvdb::Coord(
			static_cast<int>(std::round(indexPos[0])),
			static_cast<int>(std::round(indexPos[1])),
			static_cast<int>(std::round(indexPos[2]))
		));

		// If the point along the normal is INSIDE (negative distance),
		// the normal points inward, so flip the face
		if (distance < 0.0f) {
			std::reverse(faces[faceIdx].begin(), faces[faceIdx].end());
			++flippedCount;
		}
	}

	std::cout << "    VDB Oracle: flipped " << flippedCount << " faces based on SDF\n";

	return std::make_tuple(
		std::get<0>(FlattenFaces(faces)),
		std::get<1>(FlattenFaces(faces)),
		flippedCount
	);
}

ProcessedMesh VoxelizeMesh(const MeshData& input) {
	ProcessedMesh result;

	if (input.IsEmpty()) {
		return result;
	}

	openvdb::FloatGrid::Ptr grid = MeshToGrid(input);
	if (!grid) {
		return result;
	}

	if (!GridToMesh(*grid, result.points, result.faceVertexCounts, result.faceVertexIndices)) {
		result = ProcessedMesh(); // Clear partial data
		return result;
	}

	result.normals = ComputeFaceNormals(result.points, result.faceVertexCounts, result.faceVertexIndices);
	result.extent = ComputeExtentFromPoints(result.points);

	return result;
}

pxr::VtArray<pxr::GfVec3f> ComputeFaceNormals(const pxr::VtArray<pxr::GfVec3f>& points,
											  const pxr::VtArray<int>& counts,
											  const pxr::VtArray<int>& indices) {
	pxr::VtArray<pxr::GfVec3f> normals;
	if (points.empty() || counts.empty() || indices.empty()) {
		return normals;
	}

	normals.reserve(counts.size());
	std::size_t offset = 0;
	for (int vertsInFace : counts) {
		pxr::GfVec3f n(0.0f);
		if (vertsInFace >= 3 && offset + 2 < indices.size()) {
			const int aIdx = indices[offset];
			const int bIdx = indices[offset + 1];
			const int cIdx = indices[offset + 2];
			if (aIdx >= 0 && bIdx >= 0 && cIdx >= 0 &&
				static_cast<std::size_t>(std::max({aIdx, bIdx, cIdx})) < points.size()) {
				const pxr::GfVec3f ab = points[bIdx] - points[aIdx];
				const pxr::GfVec3f ac = points[cIdx] - points[aIdx];
				n = pxr::GfCross(ab, ac);
				const double length = n.GetLength();
				if (length > 1e-6) {
					n /= static_cast<float>(length);
				} else {
					n = pxr::GfVec3f(0.0f);
				}
			}
		}

		normals.push_back(n);
		offset += vertsInFace;
	}

	return normals;
}

pxr::VtArray<pxr::GfVec3f> ComputeExtentFromPoints(const pxr::VtArray<pxr::GfVec3f>& points) {
	pxr::VtArray<pxr::GfVec3f> extent;
	if (points.empty()) {
		return extent;
	}

	pxr::GfRange3d bounds;
	for (const auto& p : points) {
		bounds.UnionWith(pxr::GfVec3d(p));
	}

	if (!bounds.IsEmpty()) {
		extent.resize(2);
		extent[0] = pxr::GfVec3f(bounds.GetMin());
		extent[1] = pxr::GfVec3f(bounds.GetMax());
	}

	return extent;
}

} // namespace FixMeshNormals
