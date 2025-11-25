#include "usd_io.h"
#include "topology.h"

#include <pxr/base/gf/vec3d.h>
#include <pxr/base/tf/token.h>
#include <pxr/usd/sdf/path.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usdGeom/tokens.h>

#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace FixMeshNormals {

namespace {

template <typename VecType>
void AppendPoints(const pxr::VtArray<VecType>& src, MeshData& dst) {
	dst.points.reserve(dst.points.size() + src.size());
	for (const auto& p : src) {
		dst.points.emplace_back(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
		dst.bounds.UnionWith(pxr::GfVec3d(p[0], p[1], p[2]));
	}
}

bool ShouldSkipPrim(const pxr::UsdPrim& prim) {
	const std::string primName = prim.GetPath().GetName();

	// Skip meshes ending with "_voxelised"
	if (primName.size() >= std::string("_voxelised").size() &&
		primName.rfind("_voxelised") == primName.size() - std::string("_voxelised").size()) {
		return true;
	}

	// Skip meshes ending with "_fixed"
	if (primName.size() >= std::string("_fixed").size() &&
		primName.rfind("_fixed") == primName.size() - std::string("_fixed").size()) {
		return true;
	}

	return false;
}

std::vector<std::size_t> ComputeFaceOffsets(const pxr::VtArray<int>& counts) {
	std::vector<std::size_t> offsets;
	offsets.reserve(counts.size());
	std::size_t cursor = 0;
	for (int c : counts) {
		offsets.push_back(cursor);
		cursor += static_cast<std::size_t>(c);
	}
	return offsets;
}

std::vector<Face> BuildFaces(const pxr::VtArray<int>& counts, const pxr::VtArray<int>& indices) {
	std::vector<Face> faces;
	faces.reserve(counts.size());
	std::size_t cursor = 0;
	for (int c : counts) {
		Face f;
		for (int i = 0; i < c && cursor + i < indices.size(); ++i) {
			f.push_back(indices[cursor + i]);
		}
		faces.push_back(std::move(f));
		cursor += static_cast<std::size_t>(c);
	}
	return faces;
}

std::vector<std::vector<int>> CollectComponents(const pxr::VtArray<int>& counts,
												const pxr::VtArray<int>& indices,
												bool processElements) {
	std::vector<std::vector<int>> components;
	if (!processElements) {
		std::vector<int> single;
		single.reserve(counts.size());
		for (std::size_t i = 0; i < counts.size(); ++i) {
			single.push_back(static_cast<int>(i));
		}
		components.push_back(std::move(single));
		return components;
	}

	std::vector<Face> faces = BuildFaces(counts, indices);
	WindingResult wr = UnifyWinding(faces);
	return wr.components;
}

struct ComponentData {
	pxr::VtArray<pxr::GfVec3f> points;
	pxr::VtArray<int> counts;
	pxr::VtArray<int> indices;
};

ComponentData ExtractComponent(const pxr::VtArray<pxr::GfVec3f>& allPoints,
							   const pxr::VtArray<int>& counts,
							   const pxr::VtArray<int>& indices,
							   const std::vector<int>& faceIndices,
							   const std::vector<std::size_t>& faceOffsets) {
	ComponentData comp;
	comp.counts.reserve(faceIndices.size());

	std::vector<int> remap(allPoints.size(), -1);

	for (int faceIdx : faceIndices) {
		if (faceIdx < 0 || static_cast<std::size_t>(faceIdx) >= counts.size()) {
			continue;
		}
		int vertsInFace = counts[faceIdx];
		const std::size_t offset = faceOffsets[faceIdx];
		if (offset + static_cast<std::size_t>(vertsInFace) > indices.size()) {
			continue;
		}

		comp.counts.push_back(vertsInFace);
		for (int i = 0; i < vertsInFace; ++i) {
			const int origIndex = indices[offset + i];
			if (origIndex < 0 || static_cast<std::size_t>(origIndex) >= allPoints.size()) {
				continue;
			}
			int mapped = remap[origIndex];
			if (mapped == -1) {
				mapped = static_cast<int>(comp.points.size());
				comp.points.push_back(allPoints[origIndex]);
				remap[origIndex] = mapped;
			}
			comp.indices.push_back(mapped);
		}
	}

	return comp;
}

pxr::GfVec3f ComputePivot(const pxr::VtArray<pxr::GfVec3f>& points) {
	pxr::GfRange3d bounds;
	for (const auto& p : points) {
		bounds.UnionWith(pxr::GfVec3d(p));
	}
	if (bounds.IsEmpty()) {
		return pxr::GfVec3f(0.0f);
	}
	pxr::GfVec3d center = (bounds.GetMin() + bounds.GetMax()) * 0.5;
	return pxr::GfVec3f(center);
}

void TranslatePoints(pxr::VtArray<pxr::GfVec3f>& points, const pxr::GfVec3f& delta) {
	for (auto& p : points) {
		p += delta;
	}
}

struct WeldResult {
	pxr::VtArray<pxr::GfVec3f> points;
	pxr::VtArray<int> counts;
	pxr::VtArray<int> indices;
};

struct QuantizedKey {
	long long x;
	long long y;
	long long z;
	bool operator==(const QuantizedKey& other) const {
		return x == other.x && y == other.y && z == other.z;
	}
};

struct QuantizedKeyHash {
	std::size_t operator()(const QuantizedKey& k) const noexcept {
		std::size_t hx = std::hash<long long>{}(k.x);
		std::size_t hy = std::hash<long long>{}(k.y);
		std::size_t hz = std::hash<long long>{}(k.z);
		return hx ^ (hy << 1) ^ (hz << 2);
	}
};

WeldResult WeldMesh(const pxr::VtArray<pxr::GfVec3f>& points,
					const pxr::VtArray<int>& counts,
					const pxr::VtArray<int>& indices,
					double tolerance) {
	WeldResult res;
	if (points.empty() || counts.empty() || indices.empty()) {
		return res;
	}

	std::unordered_map<QuantizedKey, std::vector<int>, QuantizedKeyHash> bins;
	std::vector<int> mapIdx(points.size(), -1);

	for (std::size_t i = 0; i < points.size(); ++i) {
		const pxr::GfVec3f& p = points[i];
		QuantizedKey key{
			static_cast<long long>(std::llround(p[0] / tolerance)),
			static_cast<long long>(std::llround(p[1] / tolerance)),
			static_cast<long long>(std::llround(p[2] / tolerance))
		};

		int targetIdx = -1;
		auto it = bins.find(key);
		if (it != bins.end()) {
			for (int candidate : it->second) {
				pxr::GfVec3f diff = p - res.points[candidate];
				if (diff.GetLength() <= tolerance) {
					targetIdx = candidate;
					break;
				}
			}
		}

		if (targetIdx == -1) {
			targetIdx = static_cast<int>(res.points.size());
			res.points.push_back(p);
			bins[key].push_back(targetIdx);
		}

		mapIdx[i] = targetIdx;
	}

	std::size_t cursor = 0;
	for (int c : counts) {
		if (cursor + static_cast<std::size_t>(c) > indices.size()) {
			break;
		}
		pxr::VtArray<int> mappedFace;
		mappedFace.reserve(static_cast<std::size_t>(c));
		std::unordered_set<int> uniqueVerts;
		for (int i = 0; i < c; ++i) {
			int originalIdx = indices[cursor + i];
			if (originalIdx < 0 || static_cast<std::size_t>(originalIdx) >= mapIdx.size()) {
				mappedFace.clear();
				break;
			}
			int mapped = mapIdx[originalIdx];
			mappedFace.push_back(mapped);
			uniqueVerts.insert(mapped);
		}

		if (!mappedFace.empty() && static_cast<int>(uniqueVerts.size()) >= 3) {
			res.counts.push_back(c);
			for (int v : mappedFace) {
				res.indices.push_back(v);
			}
		}
		cursor += static_cast<std::size_t>(c);
	}

	return res;
}

MeshData BuildMeshData(const pxr::VtArray<pxr::GfVec3f>& points,
					   const pxr::VtArray<int>& counts,
					   const pxr::VtArray<int>& indices) {
	MeshData meshData;
	meshData.points.reserve(points.size());
	for (const auto& p : points) {
		meshData.points.emplace_back(p[0], p[1], p[2]);
		meshData.bounds.UnionWith(pxr::GfVec3d(p));
	}

	std::size_t offset = 0;
	for (std::size_t face = 0; face < counts.size(); ++face) {
		int vertsInFace = counts[face];
		if (vertsInFace < 3) {
			offset += vertsInFace;
			continue;
		}

		if (offset + static_cast<std::size_t>(vertsInFace) > indices.size()) {
			break;
		}

		const int v0 = indices[offset];
		if (vertsInFace == 3) {
			meshData.triangles.emplace_back(v0, indices[offset + 1], indices[offset + 2]);
		} else if (vertsInFace == 4) {
			meshData.quads.emplace_back(v0, indices[offset + 1], indices[offset + 2], indices[offset + 3]);
		} else {
			for (int i = 1; i + 1 < vertsInFace; ++i) {
				meshData.triangles.emplace_back(v0, indices[offset + i], indices[offset + i + 1]);
			}
		}

		offset += vertsInFace;
	}

	return meshData;
}

} // anonymous namespace

bool ReadMeshFromUsd(const pxr::UsdGeomMesh& mesh, MeshData& outData) {
	pxr::VtArray<pxr::GfVec3f> pointsf;
	pxr::UsdAttribute pointsAttr = mesh.GetPointsAttr();
	if (pointsAttr.Get(&pointsf)) {
		AppendPoints(pointsf, outData);
	} else {
		pxr::VtArray<pxr::GfVec3d> pointsd;
		if (!pointsAttr.Get(&pointsd)) {
			std::cerr << "  Mesh " << mesh.GetPath().GetString() << " has no points; skipping.\n";
			return false;
		}
		AppendPoints(pointsd, outData);
	}

	pxr::VtArray<int> counts;
	pxr::VtArray<int> indices;
	if (!mesh.GetFaceVertexCountsAttr().Get(&counts) || !mesh.GetFaceVertexIndicesAttr().Get(&indices)) {
		std::cerr << "  Mesh " << mesh.GetPath().GetString() << " missing topology; skipping.\n";
		return false;
	}

	std::size_t offset = 0;
	for (std::size_t face = 0; face < counts.size(); ++face) {
		int vertsInFace = counts[face];
		if (vertsInFace < 3) {
			offset += vertsInFace;
			continue;
		}

		if (offset + static_cast<std::size_t>(vertsInFace) > indices.size()) {
			std::cerr << "  Mesh " << mesh.GetPath().GetString() << " has inconsistent topology; skipping.\n";
			return false;
		}

		const int v0 = indices[offset];
		if (vertsInFace == 3) {
			outData.triangles.emplace_back(v0, indices[offset + 1], indices[offset + 2]);
		} else if (vertsInFace == 4) {
			outData.quads.emplace_back(v0, indices[offset + 1], indices[offset + 2], indices[offset + 3]);
		} else {
			// Fan triangulation for n-gons to keep meshToSignedDistanceField happy.
			for (int i = 1; i + 1 < vertsInFace; ++i) {
				outData.triangles.emplace_back(v0, indices[offset + i], indices[offset + i + 1]);
			}
		}

		offset += vertsInFace;
	}

	return !outData.IsEmpty();
}

bool WriteMeshToUsd(pxr::UsdStageRefPtr stage, const pxr::SdfPath& originalPath, const ProcessedMesh& mesh, bool overwrite) {
	if (mesh.points.empty()) {
		return false;
	}

	const std::string voxelisedName = originalPath.GetName() + "_voxelised";
	const pxr::SdfPath voxelisedPath = originalPath.GetParentPath().AppendChild(pxr::TfToken(voxelisedName));

	pxr::UsdPrim existingPrim = stage->GetPrimAtPath(voxelisedPath);
	if (existingPrim) {
		if (!overwrite) {
			std::cout << "    Skipping existing voxelised prim: " << voxelisedPath.GetString() << "\n";
			return false;
		}
		// Remove existing prim to overwrite
		stage->RemovePrim(voxelisedPath);
		std::cout << "    Overwriting existing voxelised prim: " << voxelisedPath.GetString() << "\n";
	}

	pxr::UsdGeomMesh voxelMesh = pxr::UsdGeomMesh::Define(stage, voxelisedPath);
	voxelMesh.CreateSubdivisionSchemeAttr().Set(pxr::UsdGeomTokens->none);
	voxelMesh.CreatePointsAttr().Set(mesh.points);
	voxelMesh.CreateFaceVertexCountsAttr().Set(mesh.faceVertexCounts);
	voxelMesh.CreateFaceVertexIndicesAttr().Set(mesh.faceVertexIndices);

	if (!mesh.normals.empty()) {
		voxelMesh.CreateNormalsAttr().Set(mesh.normals);
		voxelMesh.SetNormalsInterpolation(pxr::TfToken("uniform"));
	}

	if (!mesh.extent.empty()) {
		voxelMesh.CreateExtentAttr().Set(mesh.extent);
	}

	std::cout << "    Created voxelised mesh: " << voxelisedPath.GetString()
			  << " (" << mesh.points.size() << " verts)\n";
	return true;
}

std::size_t ProcessStageFile(const std::filesystem::path& inputPath,
							 const std::filesystem::path& outputPath,
							 bool overwrite,
							 bool processElements,
							 bool writeVoxelised,
							 double weldTolerance) {
	// Use generic_string to keep forward slashes, which USD handles more reliably on Windows.
	std::string inputPathAsString = inputPath.generic_string();
	std::cout << "Opening USD: " << inputPathAsString << "\n";

	pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(inputPathAsString);
	if (!stage) {
		std::cerr << "Failed to open USD stage: " << inputPath << "\n";
		return 0;
	}

	std::size_t createdFixed = 0;
	std::size_t createdVoxelised = 0;
	std::size_t totalComponentsProcessed = 0;
	std::cout << "Processing stage...\n";

	for (const auto& prim : stage->Traverse()) {
		if (!prim.IsA<pxr::UsdGeomMesh>()) {
			continue;
		}

		if (ShouldSkipPrim(prim)) {
			continue;
		}

		pxr::UsdGeomMesh mesh(prim);
		std::cout << "  Found mesh: " << prim.GetPath().GetString() << "\n";

		// Read mesh topology
		pxr::VtArray<pxr::GfVec3f> points;
		pxr::VtArray<int> counts;
		pxr::VtArray<int> indices;

		pxr::UsdAttribute pointsAttr = mesh.GetPointsAttr();
		pxr::VtArray<pxr::GfVec3d> pointsd;
		if (pointsAttr.Get(&points)) {
			// Already in GfVec3f format
		} else if (pointsAttr.Get(&pointsd)) {
			// Convert from GfVec3d
			points.reserve(pointsd.size());
			for (const auto& p : pointsd) {
				points.push_back(pxr::GfVec3f(p));
			}
		} else {
			std::cerr << "  Mesh " << mesh.GetPath().GetString() << " has no points; skipping.\n";
			continue;
		}

		if (!mesh.GetFaceVertexCountsAttr().Get(&counts) || !mesh.GetFaceVertexIndicesAttr().Get(&indices)) {
			std::cerr << "  Mesh " << mesh.GetPath().GetString() << " missing topology; skipping.\n";
			continue;
		}

		const std::vector<std::size_t> faceOffsets = ComputeFaceOffsets(counts);
		std::vector<std::vector<int>> components = CollectComponents(counts, indices, processElements);
		std::cout << "    Processing " << components.size() << " component(s)\n";
		std::size_t processedComponents = 0;

		// Pre-allocate combined arrays to avoid repeated reallocations
		// Estimate: original size is a reasonable upper bound for fixed mesh
		pxr::VtArray<pxr::GfVec3f> fixedPointsCombined;
		pxr::VtArray<int> fixedCountsCombined;
		pxr::VtArray<int> fixedIndicesCombined;
		fixedPointsCombined.reserve(points.size());
		fixedCountsCombined.reserve(counts.size());
		fixedIndicesCombined.reserve(indices.size());

		ProcessedMesh voxelCombined;
		// Voxelization typically produces ~100x more vertices, but we can't predict exactly
		// Reserve a reasonable initial capacity to reduce early reallocations
		if (writeVoxelised) {
			voxelCombined.points.reserve(points.size() * 10);
			voxelCombined.faceVertexCounts.reserve(counts.size() * 10);
			voxelCombined.faceVertexIndices.reserve(indices.size() * 10);
		}

		int totalOracleFlips = 0;

		// // Only show per-component progress for small numbers of components
		// const bool verboseProgress = components.size() <= 20;

		for (const auto& compFaces : components) {
			ComponentData comp = ExtractComponent(points, counts, indices, compFaces, faceOffsets);
			if (comp.points.empty() || comp.counts.empty() || comp.indices.empty()) {
				continue;
			}
			const std::size_t compIndex = processedComponents + 1;

			std::cout << "    Mesh " << prim.GetPath().GetString()
					  << " component " << compIndex << "/" << components.size()
					  << " (" << comp.counts.size() << " faces, " << comp.points.size() << " verts)\n";


			pxr::GfVec3f pivot = ComputePivot(comp.points);
			TranslatePoints(comp.points, -pivot);

			WeldResult welded = WeldMesh(comp.points, comp.counts, comp.indices, weldTolerance);
			if (welded.points.empty() || welded.counts.empty() || welded.indices.empty()) {
				continue;
			}

			// Check mesh topology after welding and fill holes if needed
			std::vector<Face> weldedFaces = BuildFaces(welded.counts, welded.indices);
			EdgeMap edgeMap = BuildEdgeMap(weldedFaces);
			ManifoldDiagnostic diag = CheckManifold(weldedFaces);

			// Fill holes if mesh is not watertight
			// Skip hole filling for flat planes (1-2 faces) - they're intentionally open geometry
			bool isFlatPlane = weldedFaces.size() <= 2;
			if (!diag.isWatertight && !isFlatPlane) {
				std::vector<pxr::GfVec3f> pointsVec(welded.points.begin(), welded.points.end());
				HoleFillResult fillResult = FillHoles(pointsVec, weldedFaces, edgeMap);

				if (fillResult.holesFilled > 0) {

					std::cout << "    Filled " << fillResult.holesFilled << " hole(s) with "
							  << fillResult.filledFaces.size() << " triangle(s)\n";


					// Add fill faces to welded mesh
					for (const auto& face : fillResult.filledFaces) {
						welded.counts.push_back(static_cast<int>(face.size()));
						for (int idx : face) {
							welded.indices.push_back(idx);
						}
					}

					// Re-check topology after hole filling
					weldedFaces = BuildFaces(welded.counts, welded.indices);
					edgeMap = BuildEdgeMap(weldedFaces);
					diag = CheckManifold(weldedFaces);
				}
			}

			if (!diag.isWatertight || !diag.isManifold) {
				if (!isFlatPlane) {
					PrintManifoldDiagnostic(diag, prim.GetPath().GetString() + " (after hole fill)");
				}
			}

			auto [windingCounts, windingIndices] = FixWindingAndOrientation(welded.points, welded.counts, welded.indices);

			// Skip SDF/voxelization for flat planes - they're intentionally 2D geometry
			pxr::VtArray<int> oracleCounts = windingCounts;
			pxr::VtArray<int> oracleIndices = windingIndices;
			int oracleFlips = 0;

			if (!isFlatPlane) {
				MeshData sdfMesh = BuildMeshData(welded.points, windingCounts, windingIndices);
				if (!sdfMesh.IsEmpty()) {
					openvdb::FloatGrid::Ptr sdfGrid = CreateSDFFromMesh(sdfMesh);
					if (sdfGrid) {
						auto [oCounts, oIndices, oFlips] = FixOrientationWithVDBOracle(
							welded.points, windingCounts, windingIndices, sdfGrid);
						oracleCounts = oCounts;
						oracleIndices = oIndices;
						oracleFlips = oFlips;
					}
				}
			}
			totalOracleFlips += oracleFlips;

			// Build fixed geometry (restore pivot)
			pxr::VtArray<pxr::GfVec3f> fixedPoints = welded.points;
			TranslatePoints(fixedPoints, pivot);

			std::size_t baseFixed = fixedPointsCombined.size();
			// Bulk append points
			std::size_t oldPointsSize = fixedPointsCombined.size();
			fixedPointsCombined.resize(oldPointsSize + fixedPoints.size());
			std::copy(fixedPoints.begin(), fixedPoints.end(), fixedPointsCombined.begin() + oldPointsSize);
			// Bulk append counts
			std::size_t oldCountsSize = fixedCountsCombined.size();
			fixedCountsCombined.resize(oldCountsSize + oracleCounts.size());
			std::copy(oracleCounts.begin(), oracleCounts.end(), fixedCountsCombined.begin() + oldCountsSize);
			// Indices need offset adjustment - bulk resize then fill
			std::size_t oldIndicesSize = fixedIndicesCombined.size();
			fixedIndicesCombined.resize(oldIndicesSize + oracleIndices.size());
			for (std::size_t i = 0; i < oracleIndices.size(); ++i) {
				fixedIndicesCombined[oldIndicesSize + i] = static_cast<int>(baseFixed) + oracleIndices[i];
			}

			// Voxelisation - skip for flat planes (intentionally 2D geometry) or if not writing voxelised output
			if (!isFlatPlane && writeVoxelised) {
				MeshData voxelMeshInput = BuildMeshData(welded.points, oracleCounts, oracleIndices);
				ProcessedMesh voxelComp = VoxelizeMesh(voxelMeshInput);
				if (!voxelComp.points.empty()) {
					TranslatePoints(voxelComp.points, pivot);
					// Recompute extent after translation
					voxelComp.extent = ComputeExtentFromPoints(voxelComp.points);

					std::size_t baseVoxel = voxelCombined.points.size();
					// Bulk append points
					std::size_t oldVoxelPointsSize = voxelCombined.points.size();
					voxelCombined.points.resize(oldVoxelPointsSize + voxelComp.points.size());
					std::copy(voxelComp.points.begin(), voxelComp.points.end(),
						voxelCombined.points.begin() + oldVoxelPointsSize);
					// Bulk append counts
					std::size_t oldVoxelCountsSize = voxelCombined.faceVertexCounts.size();
					voxelCombined.faceVertexCounts.resize(oldVoxelCountsSize + voxelComp.faceVertexCounts.size());
					std::copy(voxelComp.faceVertexCounts.begin(), voxelComp.faceVertexCounts.end(),
						voxelCombined.faceVertexCounts.begin() + oldVoxelCountsSize);
					// Indices need offset adjustment - bulk resize then fill
					std::size_t oldVoxelIndicesSize = voxelCombined.faceVertexIndices.size();
					voxelCombined.faceVertexIndices.resize(oldVoxelIndicesSize + voxelComp.faceVertexIndices.size());
					for (std::size_t i = 0; i < voxelComp.faceVertexIndices.size(); ++i) {
						voxelCombined.faceVertexIndices[oldVoxelIndicesSize + i] =
							static_cast<int>(baseVoxel) + voxelComp.faceVertexIndices[i];
					}
				}
			}
			++processedComponents;
		}

		
		std::cout << "                                                              \r";
	

		if (fixedPointsCombined.empty() || fixedCountsCombined.empty() || fixedIndicesCombined.empty()) {
			std::cerr << "  Mesh " << mesh.GetPath().GetString() << " has no valid processed geometry; skipping.\n";
			continue;
		}
		std::cout << "    Processed " << processedComponents << " component(s)\n";
		totalComponentsProcessed += processedComponents;

		// Compute face-varying normals for fixed mesh
		pxr::VtArray<pxr::GfVec3f> faceVaryingNormals;
		std::size_t offset = 0;
		for (int faceCount : fixedCountsCombined) {
			if (offset + static_cast<std::size_t>(faceCount) > fixedIndicesCombined.size()) {
				break;
			}
			pxr::GfVec3f faceNormal(0.0f);
			if (faceCount >= 3 && offset + 2 < fixedIndicesCombined.size()) {
				const pxr::GfVec3f& p0 = fixedPointsCombined[fixedIndicesCombined[offset]];
				const pxr::GfVec3f& p1 = fixedPointsCombined[fixedIndicesCombined[offset + 1]];
				const pxr::GfVec3f& p2 = fixedPointsCombined[fixedIndicesCombined[offset + 2]];
				pxr::GfVec3f n = pxr::GfCross(p1 - p0, p2 - p0);
				float len = n.GetLength();
				if (len > 1e-6f) {
					faceNormal = n / len;
				}
			}
			for (int i = 0; i < faceCount; ++i) {
				faceVaryingNormals.push_back(faceNormal);
			}
			offset += faceCount;
		}

		const std::string fixedName = prim.GetPath().GetName() + "_fixed";
		const pxr::SdfPath fixedPath = prim.GetPath().GetParentPath().AppendChild(pxr::TfToken(fixedName));
		pxr::UsdPrim existingFixed = stage->GetPrimAtPath(fixedPath);
		if (existingFixed) {
			if (!overwrite) {
				std::cout << "    Skipping existing fixed prim: " << fixedPath.GetString() << "\n";
			} else {
				stage->RemovePrim(fixedPath);
				std::cout << "    Overwriting existing fixed prim: " << fixedPath.GetString() << "\n";
			}
		}

		if (!existingFixed || overwrite) {
			pxr::UsdGeomMesh fixedMesh = pxr::UsdGeomMesh::Define(stage, fixedPath);
			fixedMesh.CreateSubdivisionSchemeAttr().Set(pxr::UsdGeomTokens->none);
			fixedMesh.CreatePointsAttr().Set(fixedPointsCombined);
			fixedMesh.CreateFaceVertexCountsAttr().Set(fixedCountsCombined);
			fixedMesh.CreateFaceVertexIndicesAttr().Set(fixedIndicesCombined);
			fixedMesh.GetNormalsAttr().Set(faceVaryingNormals);
			fixedMesh.SetNormalsInterpolation(pxr::TfToken("faceVarying"));

			pxr::VtArray<pxr::GfVec3f> extent = ComputeExtentFromPoints(fixedPointsCombined);
			if (!extent.empty()) {
				fixedMesh.CreateExtentAttr().Set(extent);
			}

			std::cout << "    Created fixed mesh: " << fixedPath.GetString()
					  << " (" << totalOracleFlips << " faces flipped by oracle, "
					  << fixedPointsCombined.size() << " verts)\n";
			++createdFixed;
		}

		if (writeVoxelised && !voxelCombined.points.empty()) {
			voxelCombined.extent = ComputeExtentFromPoints(voxelCombined.points);
			voxelCombined.normals = ComputeFaceNormals(voxelCombined.points, voxelCombined.faceVertexCounts, voxelCombined.faceVertexIndices);
			if (WriteMeshToUsd(stage, prim.GetPath(), voxelCombined, overwrite)) {
				++createdVoxelised;
			}
		}
	}

	std::cout << "Processing complete. Created " << createdFixed << " fixed mesh(es)"
			  << " and " << createdVoxelised << " voxelised mesh(es). "
			  << "Processed " << totalComponentsProcessed << " component(s) total.\n";

	const std::string outputPathAsString = outputPath.generic_string();
	if (!stage->Export(outputPathAsString)) {
		std::cerr << "Failed to export stage to: " << outputPathAsString << "\n";
		return 0;
	}

	std::cout << "Wrote voxelised stage to: " << outputPathAsString << "\n";
	return createdFixed;
}

} // namespace FixMeshNormals
