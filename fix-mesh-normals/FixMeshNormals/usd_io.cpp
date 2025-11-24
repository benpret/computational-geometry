#include "usd_io.h"

#include <pxr/base/gf/vec3d.h>
#include <pxr/base/tf/token.h>
#include <pxr/usd/sdf/path.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usdGeom/tokens.h>

#include <iostream>

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

std::size_t ProcessStageFile(const std::filesystem::path& inputPath, const std::filesystem::path& outputPath, bool overwrite) {
	// Use generic_string to keep forward slashes, which USD handles more reliably on Windows.
	std::string inputPathAsString = inputPath.generic_string();
	std::cout << "Opening USD: " << inputPathAsString << "\n";

	pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(inputPathAsString);
	if (!stage) {
		std::cerr << "Failed to open USD stage: " << inputPath << "\n";
		return 0;
	}

	std::size_t createdMeshes = 0;
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

		// Step 1: Fix winding consistency
		auto [windingFixedCounts, windingFixedIndices] = FixWindingAndOrientation(points, counts, indices);

		// Step 2: Convert to MeshData for SDF creation
		MeshData meshData;
		meshData.points.reserve(points.size());
		for (const auto& p : points) {
			meshData.points.emplace_back(p[0], p[1], p[2]);
			meshData.bounds.UnionWith(pxr::GfVec3d(p));
		}

		// Convert topology to triangles/quads for VDB
		std::size_t offset = 0;
		for (std::size_t face = 0; face < windingFixedCounts.size(); ++face) {
			int vertsInFace = windingFixedCounts[face];
			if (vertsInFace < 3) {
				offset += vertsInFace;
				continue;
			}

			if (offset + static_cast<std::size_t>(vertsInFace) > windingFixedIndices.size()) {
				std::cerr << "  Mesh " << mesh.GetPath().GetString() << " has inconsistent topology; skipping.\n";
				continue;
			}

			const int v0 = windingFixedIndices[offset];
			if (vertsInFace == 3) {
				meshData.triangles.emplace_back(v0, windingFixedIndices[offset + 1], windingFixedIndices[offset + 2]);
			} else if (vertsInFace == 4) {
				meshData.quads.emplace_back(v0, windingFixedIndices[offset + 1], windingFixedIndices[offset + 2], windingFixedIndices[offset + 3]);
			} else {
				// Fan triangulation for n-gons
				for (int i = 1; i + 1 < vertsInFace; ++i) {
					meshData.triangles.emplace_back(v0, windingFixedIndices[offset + i], windingFixedIndices[offset + i + 1]);
				}
			}

			offset += vertsInFace;
		}

		if (meshData.IsEmpty()) {
			std::cerr << "  Mesh " << mesh.GetPath().GetString() << " has no valid geometry; skipping.\n";
			continue;
		}

		// Step 3: Create SDF from winding-fixed mesh (as oracle)
		openvdb::FloatGrid::Ptr sdfGrid = CreateSDFFromMesh(meshData);
		if (!sdfGrid) {
			std::cerr << "    Failed to create SDF for " << prim.GetPath().GetString() << "\n";
			continue;
		}

		// Step 4: Use SDF oracle to fix original mesh face orientations
		auto [oracleCounts, oracleIndices, oracleFlips] = FixOrientationWithVDBOracle(
			points, windingFixedCounts, windingFixedIndices, sdfGrid);

		// Step 5: Create "_fixed" mesh with corrected topology
		const std::string fixedName = prim.GetPath().GetName() + "_fixed";
		const pxr::SdfPath fixedPath = prim.GetPath().GetParentPath().AppendChild(pxr::TfToken(fixedName));

		// Check if _fixed mesh already exists
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
			// Create new _fixed mesh
			pxr::UsdGeomMesh fixedMesh = pxr::UsdGeomMesh::Define(stage, fixedPath);
			fixedMesh.CreateSubdivisionSchemeAttr().Set(pxr::UsdGeomTokens->none);
			fixedMesh.CreatePointsAttr().Set(points);
			fixedMesh.CreateFaceVertexCountsAttr().Set(oracleCounts);
			fixedMesh.CreateFaceVertexIndicesAttr().Set(oracleIndices);

			// Compute and set face-varying normals for the corrected mesh
			pxr::VtArray<pxr::GfVec3f> faceVaryingNormals;
			offset = 0;
			for (int faceCount : oracleCounts) {
				pxr::GfVec3f faceNormal(0.0f);
				if (faceCount >= 3 && offset + 2 < oracleIndices.size()) {
					const pxr::GfVec3f& p0 = points[oracleIndices[offset]];
					const pxr::GfVec3f& p1 = points[oracleIndices[offset + 1]];
					const pxr::GfVec3f& p2 = points[oracleIndices[offset + 2]];
					pxr::GfVec3f n = pxr::GfCross(p1 - p0, p2 - p0);
					float len = n.GetLength();
					if (len > 1e-6f) {
						faceNormal = n / len;
					}
				}
				// Repeat normal for each vertex in the face (face-varying)
				for (int i = 0; i < faceCount; ++i) {
					faceVaryingNormals.push_back(faceNormal);
				}
				offset += faceCount;
			}

			fixedMesh.GetNormalsAttr().Set(faceVaryingNormals);
			fixedMesh.SetNormalsInterpolation(pxr::TfToken("faceVarying"));

			// Compute extent
			pxr::GfRange3d bounds;
			for (const auto& p : points) {
				bounds.UnionWith(pxr::GfVec3d(p));
			}
			if (!bounds.IsEmpty()) {
				pxr::VtArray<pxr::GfVec3f> extent;
				extent.resize(2);
				extent[0] = pxr::GfVec3f(bounds.GetMin());
				extent[1] = pxr::GfVec3f(bounds.GetMax());
				fixedMesh.CreateExtentAttr().Set(extent);
			}

			std::cout << "    Created fixed mesh: " << fixedPath.GetString()
					  << " (" << oracleFlips << " faces flipped by oracle, " << points.size() << " verts)\n";
			++createdMeshes;
		}

		// Step 6: Optionally create voxelized mesh as well
		ProcessedMesh voxelized = VoxelizeMesh(meshData);
		if (!voxelized.points.empty()) {
			WriteMeshToUsd(stage, prim.GetPath(), voxelized, overwrite);
		}
	}

	std::cout << "Processing complete. Created " << createdMeshes << " voxelised mesh(es).\n";

	const std::string outputPathAsString = outputPath.generic_string();
	if (!stage->Export(outputPathAsString)) {
		std::cerr << "Failed to export stage to: " << outputPathAsString << "\n";
		return 0;
	}

	std::cout << "Wrote voxelised stage to: " << outputPathAsString << "\n";
	return createdMeshes;
}

} // namespace FixMeshNormals
