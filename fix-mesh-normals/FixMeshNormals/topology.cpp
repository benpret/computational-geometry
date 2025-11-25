#include "topology.h"

#include <pxr/base/gf/vec3d.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <stack>

namespace FixMeshNormals {

EdgeMap BuildEdgeMap(const std::vector<Face>& faces) {
	EdgeMap edgeMap;
	for (std::size_t faceIdx = 0; faceIdx < faces.size(); ++faceIdx) {
		const Face& face = faces[faceIdx];
		if (face.size() < 2) {
			continue;
		}
		for (std::size_t i = 0; i < face.size(); ++i) {
			int v0 = face[i];
			int v1 = face[(i + 1) % face.size()];
			EdgeKey key = (v0 < v1) ? EdgeKey(v0, v1) : EdgeKey(v1, v0);
			edgeMap[key].emplace_back(static_cast<int>(faceIdx), v0, v1);
		}
	}
	return edgeMap;
}

WindingResult UnifyWinding(const std::vector<Face>& faces) {
	EdgeMap edgeMap = BuildEdgeMap(faces);
	std::vector<int> orientation(faces.size(), 0);  // 0 = unset, +1/-1 = orientation
	std::vector<std::vector<int>> components;
	int conflicts = 0;

	for (std::size_t seed = 0; seed < faces.size(); ++seed) {
		if (orientation[seed] != 0) {
			continue;
		}

		orientation[seed] = 1;
		std::stack<int> stack;
		stack.push(static_cast<int>(seed));
		std::vector<int> component;

		while (!stack.empty()) {
			int faceIdx = stack.top();
			stack.pop();
			component.push_back(faceIdx);

			const Face& face = faces[faceIdx];
			int sign = orientation[faceIdx];
			if (face.size() < 2) {
				continue;
			}

			for (std::size_t i = 0; i < face.size(); ++i) {
				int v0 = face[i];
				int v1 = face[(i + 1) % face.size()];

				// Actual edge direction considering current orientation
				std::pair<int, int> actualEdge = (sign == 1) ? std::make_pair(v0, v1) : std::make_pair(v1, v0);
				EdgeKey key = (v0 < v1) ? EdgeKey(v0, v1) : EdgeKey(v1, v0);

				auto it = edgeMap.find(key);
				if (it == edgeMap.end()) {
					continue;
				}

				for (const auto& [neighborIdx, nbV0, nbV1] : it->second) {
					if (neighborIdx == faceIdx) {
						continue;
					}

					// Neighbor should traverse this edge in opposite direction
					std::pair<int, int> desired(actualEdge.second, actualEdge.first);
					std::pair<int, int> baseEdge(nbV0, nbV1);

					int neededOrientation = (baseEdge == desired) ? 1 : -1;

					if (orientation[neighborIdx] == 0) {
						orientation[neighborIdx] = neededOrientation;
						stack.push(neighborIdx);
					} else if (orientation[neighborIdx] != neededOrientation) {
						++conflicts;
					}
				}
			}
		}

		components.push_back(component);
	}

	// Fill any still-unset faces (isolated edges/points)
	for (std::size_t idx = 0; idx < orientation.size(); ++idx) {
		if (orientation[idx] == 0) {
			orientation[idx] = 1;
		}
	}

	return WindingResult{orientation, components, conflicts};
}

std::vector<Face> ApplyOrientation(const std::vector<Face>& faces, const std::vector<int>& orientation) {
	std::vector<Face> adjusted;
	adjusted.reserve(faces.size());
	for (std::size_t i = 0; i < faces.size(); ++i) {
		if (orientation[i] == -1) {
			Face reversed(faces[i].rbegin(), faces[i].rend());
			adjusted.push_back(reversed);
		} else {
			adjusted.push_back(faces[i]);
		}
	}
	return adjusted;
}

pxr::GfVec3d PolygonAreaVector(const std::vector<pxr::GfVec3f>& points, const Face& face) {
	if (face.size() < 3) {
		return pxr::GfVec3d(0.0);
	}

	pxr::GfVec3d p0(points[face[0]]);
	pxr::GfVec3d vec(0.0);

	for (std::size_t i = 1; i + 1 < face.size(); ++i) {
		pxr::GfVec3d p1(points[face[i]]);
		pxr::GfVec3d p2(points[face[i + 1]]);
		vec += pxr::GfCross(p1 - p0, p2 - p0);
	}

	return vec;
}

double SignedVolume(const std::vector<pxr::GfVec3f>& points,
					const std::vector<Face>& faces,
					const std::vector<int>& faceIndices) {
	double vol = 0.0;
	for (int faceIdx : faceIndices) {
		const Face& face = faces[faceIdx];
		if (face.size() < 3) {
			continue;
		}

		pxr::GfVec3d p0(points[face[0]]);
		for (std::size_t i = 1; i + 1 < face.size(); ++i) {
			pxr::GfVec3d p1(points[face[i]]);
			pxr::GfVec3d p2(points[face[i + 1]]);
			vol += pxr::GfDot(p0, pxr::GfCross(p1, p2)) / 6.0;
		}
	}
	return vol;
}

double CentroidFlux(const std::vector<pxr::GfVec3f>& points,
					const std::vector<Face>& faces,
					const std::vector<int>& faceIndices) {
	// Collect unique vertices in this component
	std::set<int> usedVertices;
	for (int faceIdx : faceIndices) {
		for (int v : faces[faceIdx]) {
			usedVertices.insert(v);
		}
	}

	if (usedVertices.empty()) {
		return 0.0;
	}

	// Compute component centroid
	pxr::GfVec3d compCentroid(0.0);
	for (int v : usedVertices) {
		compCentroid += pxr::GfVec3d(points[v]);
	}
	compCentroid /= static_cast<double>(usedVertices.size());

	double flux = 0.0;
	for (int faceIdx : faceIndices) {
		const Face& face = faces[faceIdx];
		if (face.size() < 3) {
			continue;
		}

		// Face centroid
		pxr::GfVec3d faceCentroid(0.0);
		for (int v : face) {
			faceCentroid += pxr::GfVec3d(points[v]);
		}
		faceCentroid /= static_cast<double>(face.size());

		pxr::GfVec3d areaVec = PolygonAreaVector(points, face);
		flux += pxr::GfDot(areaVec, faceCentroid - compCentroid);
	}

	return flux;
}

std::pair<int, int> FaceRadialSigns(const std::vector<pxr::GfVec3f>& points,
									const std::vector<Face>& faces,
									const std::vector<int>& faceIndices) {
	pxr::GfVec3d compCentroid = ComponentCentroid(points, faces, faceIndices);

	int pos = 0;
	int neg = 0;

	for (int faceIdx : faceIndices) {
		const Face& face = faces[faceIdx];
		if (face.size() < 3) {
			continue;
		}

		pxr::GfVec3d areaVec = PolygonAreaVector(points, face);

		pxr::GfVec3d faceCentroid(0.0);
		for (int v : face) {
			faceCentroid += pxr::GfVec3d(points[v]);
		}
		faceCentroid /= static_cast<double>(face.size());

		double dot = pxr::GfDot(areaVec, faceCentroid - compCentroid);
		if (dot > 0) {
			++pos;
		} else if (dot < 0) {
			++neg;
		}
	}

	return {pos, neg};
}

pxr::GfVec3d ComponentCentroid(const std::vector<pxr::GfVec3f>& points,
							   const std::vector<Face>& faces,
							   const std::vector<int>& faceIndices) {
	std::set<int> usedVertices;
	for (int faceIdx : faceIndices) {
		for (int v : faces[faceIdx]) {
			usedVertices.insert(v);
		}
	}

	if (usedVertices.empty()) {
		// Fallback to all points
		pxr::GfVec3d centroid(0.0);
		for (const auto& p : points) {
			centroid += pxr::GfVec3d(p);
		}
		return centroid / static_cast<double>(points.size());
	}

	pxr::GfVec3d centroid(0.0);
	for (int v : usedVertices) {
		centroid += pxr::GfVec3d(points[v]);
	}
	return centroid / static_cast<double>(usedVertices.size());
}

ManifoldDiagnostic CheckManifold(const std::vector<Face>& faces) {
	EdgeMap edgeMap = BuildEdgeMap(faces);

	ManifoldDiagnostic diag{};
	diag.totalEdges = static_cast<int>(edgeMap.size());

	for (const auto& [edge, faceRefs] : edgeMap) {
		int faceCount = static_cast<int>(faceRefs.size());
		if (faceCount == 1) {
			++diag.boundaryEdges;
			diag.boundaryEdgeList.push_back(edge);
		} else if (faceCount == 2) {
			++diag.manifoldEdges;
		} else {
			++diag.nonManifoldEdges;
			diag.nonManifoldEdgeList.push_back(edge);
		}
	}

	diag.isWatertight = (diag.boundaryEdges == 0);
	diag.isManifold = (diag.nonManifoldEdges == 0);

	return diag;
}

void PrintManifoldDiagnostic(const ManifoldDiagnostic& diag, const std::string& meshName) {
	std::cout << "    Manifold check for " << meshName << ":\n";
	std::cout << "      Total edges: " << diag.totalEdges << "\n";
	std::cout << "      Manifold edges (2 faces): " << diag.manifoldEdges << "\n";
	std::cout << "      Boundary edges (1 face): " << diag.boundaryEdges << "\n";
	std::cout << "      Non-manifold edges (3+ faces): " << diag.nonManifoldEdges << "\n";
	std::cout << "      Watertight: " << (diag.isWatertight ? "YES" : "NO") << "\n";
	std::cout << "      Manifold: " << (diag.isManifold ? "YES" : "NO") << "\n";

	if (!diag.isWatertight) {
		std::cout << "      WARNING: Mesh has " << diag.boundaryEdges
				  << " boundary edge(s) - voxelization may fail!\n";
	}
}

} // namespace FixMeshNormals
