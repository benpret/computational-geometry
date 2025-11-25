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

namespace {

// Build a map from vertex to connected boundary vertices (undirected)
// Returns both the adjacency map and a set of all boundary edges
std::pair<std::map<int, std::set<int>>, std::set<EdgeKey>> BuildBoundaryGraph(
	const std::vector<Face>& faces,
	const EdgeMap& edgeMap) {

	std::map<int, std::set<int>> adjacency;
	std::set<EdgeKey> boundaryEdges;

	for (const auto& [edgeKey, faceRefs] : edgeMap) {
		if (faceRefs.size() != 1) {
			continue;  // Not a boundary edge
		}

		int v0 = edgeKey.first;
		int v1 = edgeKey.second;

		// Store undirected adjacency
		adjacency[v0].insert(v1);
		adjacency[v1].insert(v0);
		boundaryEdges.insert(edgeKey);
	}

	return {adjacency, boundaryEdges};
}

// Trace a boundary loop starting from a vertex using undirected graph
// Returns the loop as a list of vertex indices, or empty if no valid loop
std::vector<int> TraceBoundaryLoopUndirected(int startVertex,
											  std::map<int, std::set<int>>& adjacency,
											  std::set<int>& visitedVertices) {
	std::vector<int> loop;
	int current = startVertex;
	int prev = -1;

	while (true) {
		loop.push_back(current);
		visitedVertices.insert(current);

		auto it = adjacency.find(current);
		if (it == adjacency.end() || it->second.empty()) {
			// Dead end - not a valid loop
			return {};
		}

		// Find next vertex (not the one we came from)
		int next = -1;
		for (int neighbor : it->second) {
			if (neighbor != prev) {
				next = neighbor;
				break;
			}
		}

		if (next == -1) {
			// No unvisited neighbor
			return {};
		}

		if (next == startVertex && loop.size() >= 3) {
			// Completed the loop
			break;
		}

		if (visitedVertices.count(next) && next != startVertex) {
			// Hit a visited vertex that's not the start - complex topology
			return {};
		}

		prev = current;
		current = next;

		// Safety check for infinite loops
		if (loop.size() > 10000) {
			return {};
		}
	}

	return loop;
}

// Find all boundary loops in the mesh
std::vector<std::vector<int>> FindBoundaryLoops(const std::vector<Face>& faces,
												const EdgeMap& edgeMap) {
	auto [adjacency, boundaryEdges] = BuildBoundaryGraph(faces, edgeMap);
	std::set<int> visitedVertices;
	std::vector<std::vector<int>> loops;

	// Find all starting vertices (any vertex with boundary edge connections)
	for (const auto& [v, neighbors] : adjacency) {
		if (visitedVertices.count(v)) {
			continue;
		}

		// Each boundary vertex should have exactly 2 neighbors for a simple loop
		if (neighbors.size() != 2) {
			// Complex topology - skip for now
			continue;
		}

		std::vector<int> loop = TraceBoundaryLoopUndirected(v, adjacency, visitedVertices);
		if (loop.size() >= 3) {
			loops.push_back(loop);
		}
	}

	return loops;
}

// Compute the normal of a boundary loop based on its vertices
pxr::GfVec3d ComputeLoopNormal(const std::vector<pxr::GfVec3f>& points,
							   const std::vector<int>& loop) {
	if (loop.size() < 3) {
		return pxr::GfVec3d(0, 0, 1);
	}

	// Use Newell's method for robust normal computation
	pxr::GfVec3d normal(0.0);
	for (std::size_t i = 0; i < loop.size(); ++i) {
		const pxr::GfVec3f& curr = points[loop[i]];
		const pxr::GfVec3f& next = points[loop[(i + 1) % loop.size()]];

		normal[0] += (curr[1] - next[1]) * (curr[2] + next[2]);
		normal[1] += (curr[2] - next[2]) * (curr[0] + next[0]);
		normal[2] += (curr[0] - next[0]) * (curr[1] + next[1]);
	}

	double len = normal.GetLength();
	if (len > 1e-10) {
		normal /= len;
	}

	return normal;
}

// Determine the correct winding for a fill face based on adjacent face normals
// Returns true if the loop should be reversed
bool ShouldReverseLoop(const std::vector<pxr::GfVec3f>& points,
					   const std::vector<Face>& faces,
					   const std::vector<int>& loop,
					   const EdgeMap& edgeMap) {
	if (loop.size() < 3) {
		return false;
	}

	// Find an adjacent face to determine correct winding
	// Look at the first edge of the loop
	int v0 = loop[0];
	int v1 = loop[1];
	EdgeKey key = (v0 < v1) ? EdgeKey(v0, v1) : EdgeKey(v1, v0);

	auto it = edgeMap.find(key);
	if (it == edgeMap.end() || it->second.empty()) {
		return false;
	}

	// Get the adjacent face's edge direction
	const auto& [faceIdx, adjV0, adjV1] = it->second[0];

	// The boundary loop should have opposite winding to the adjacent face
	// If the adjacent face has edge (v0->v1), the fill should have (v1->v0)
	// The loop is traced following boundary edges, so it already has the face's winding
	// We need to reverse it to get the opposite winding (facing outward to fill the hole)

	// Actually, boundary edges are traced in the direction they appear in faces
	// So the loop has the same winding as adjacent faces
	// To create a fill face that closes the hole, we need opposite winding
	return true;
}

// Triangulate a polygon loop using ear clipping
// Returns triangles as faces with 3 vertices each
std::vector<Face> TriangulateLoop(const std::vector<pxr::GfVec3f>& points,
								  const std::vector<int>& loop) {
	std::vector<Face> triangles;

	if (loop.size() < 3) {
		return triangles;
	}

	if (loop.size() == 3) {
		triangles.push_back({loop[0], loop[1], loop[2]});
		return triangles;
	}

	// Simple fan triangulation for convex or nearly-convex polygons
	// For more complex holes, ear clipping would be better
	// Fan from first vertex
	for (std::size_t i = 1; i + 1 < loop.size(); ++i) {
		triangles.push_back({loop[0], loop[i], loop[i + 1]});
	}

	return triangles;
}

// More robust ear clipping triangulation
std::vector<Face> EarClipTriangulate(const std::vector<pxr::GfVec3f>& points,
									 std::vector<int> loop,
									 const pxr::GfVec3d& normal) {
	std::vector<Face> triangles;

	if (loop.size() < 3) {
		return triangles;
	}

	if (loop.size() == 3) {
		triangles.push_back({loop[0], loop[1], loop[2]});
		return triangles;
	}

	// Project points onto 2D plane for ear detection
	// Find two axes perpendicular to the normal
	pxr::GfVec3d axisU, axisV;
	if (std::abs(normal[0]) < 0.9) {
		axisU = pxr::GfCross(normal, pxr::GfVec3d(1, 0, 0));
	} else {
		axisU = pxr::GfCross(normal, pxr::GfVec3d(0, 1, 0));
	}
	axisU.Normalize();
	axisV = pxr::GfCross(normal, axisU);
	axisV.Normalize();

	// Project loop vertices to 2D
	auto project2D = [&](int vertIdx) -> std::pair<double, double> {
		pxr::GfVec3d p(points[vertIdx]);
		return {pxr::GfDot(p, axisU), pxr::GfDot(p, axisV)};
	};

	// Check if a triangle is counter-clockwise in 2D
	auto isCounterClockwise = [&](int i0, int i1, int i2) -> bool {
		auto [x0, y0] = project2D(loop[i0]);
		auto [x1, y1] = project2D(loop[i1]);
		auto [x2, y2] = project2D(loop[i2]);
		double cross = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
		return cross > 0;
	};

	// Check if point is inside triangle
	auto pointInTriangle = [&](int pi, int i0, int i1, int i2) -> bool {
		auto [px, py] = project2D(loop[pi]);
		auto [x0, y0] = project2D(loop[i0]);
		auto [x1, y1] = project2D(loop[i1]);
		auto [x2, y2] = project2D(loop[i2]);

		auto sign = [](double x1, double y1, double x2, double y2, double x3, double y3) {
			return (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
		};

		double d1 = sign(px, py, x0, y0, x1, y1);
		double d2 = sign(px, py, x1, y1, x2, y2);
		double d3 = sign(px, py, x2, y2, x0, y0);

		bool hasNeg = (d1 < 0) || (d2 < 0) || (d3 < 0);
		bool hasPos = (d1 > 0) || (d2 > 0) || (d3 > 0);

		return !(hasNeg && hasPos);
	};

	// Ear clipping
	std::vector<int> remaining;
	for (std::size_t i = 0; i < loop.size(); ++i) {
		remaining.push_back(static_cast<int>(i));
	}

	int maxIterations = static_cast<int>(loop.size() * loop.size());
	int iterations = 0;

	while (remaining.size() > 3 && iterations < maxIterations) {
		++iterations;
		bool earFound = false;

		for (std::size_t i = 0; i < remaining.size(); ++i) {
			std::size_t prevIdx = (i + remaining.size() - 1) % remaining.size();
			std::size_t nextIdx = (i + 1) % remaining.size();

			int i0 = remaining[prevIdx];
			int i1 = remaining[i];
			int i2 = remaining[nextIdx];

			// Check if this is a convex vertex (ear candidate)
			if (!isCounterClockwise(i0, i1, i2)) {
				continue;
			}

			// Check if any other vertex is inside this triangle
			bool isEar = true;
			for (std::size_t j = 0; j < remaining.size(); ++j) {
				if (j == prevIdx || j == i || j == nextIdx) {
					continue;
				}
				if (pointInTriangle(remaining[j], i0, i1, i2)) {
					isEar = false;
					break;
				}
			}

			if (isEar) {
				triangles.push_back({loop[i0], loop[i1], loop[i2]});
				remaining.erase(remaining.begin() + static_cast<std::ptrdiff_t>(i));
				earFound = true;
				break;
			}
		}

		if (!earFound) {
			// No ear found - fall back to fan triangulation for remaining vertices
			break;
		}
	}

	// Handle remaining vertices (3 or couldn't find ears)
	if (remaining.size() == 3) {
		triangles.push_back({loop[remaining[0]], loop[remaining[1]], loop[remaining[2]]});
	} else if (remaining.size() > 3) {
		// Fallback: fan triangulation
		for (std::size_t i = 1; i + 1 < remaining.size(); ++i) {
			triangles.push_back({loop[remaining[0]], loop[remaining[i]], loop[remaining[i + 1]]});
		}
	}

	return triangles;
}

} // anonymous namespace

HoleFillResult FillHoles(const std::vector<pxr::GfVec3f>& points,
						 const std::vector<Face>& faces,
						 const EdgeMap& edgeMap) {
	HoleFillResult result;
	result.holesFound = 0;
	result.holesFilled = 0;

	// Find all boundary loops
	std::vector<std::vector<int>> loops = FindBoundaryLoops(faces, edgeMap);
	result.holesFound = static_cast<int>(loops.size());

	for (const auto& loop : loops) {
		if (loop.size() < 3) {
			continue;
		}

		// Compute loop normal for triangulation
		pxr::GfVec3d loopNormal = ComputeLoopNormal(points, loop);

		// Determine if we need to reverse the loop for correct winding
		std::vector<int> fillLoop = loop;
		if (ShouldReverseLoop(points, faces, loop, edgeMap)) {
			std::reverse(fillLoop.begin(), fillLoop.end());
			loopNormal = -loopNormal;
		}

		// Triangulate the loop
		std::vector<Face> fillTriangles = EarClipTriangulate(points, fillLoop, loopNormal);

		if (!fillTriangles.empty()) {
			for (const auto& tri : fillTriangles) {
				result.filledFaces.push_back(tri);
			}
			++result.holesFilled;
		}
	}

	return result;
}

HoleFillResult FillHoles(const std::vector<pxr::GfVec3f>& points,
						 const std::vector<Face>& faces) {
	EdgeMap edgeMap = BuildEdgeMap(faces);
	return FillHoles(points, faces, edgeMap);
}

} // namespace FixMeshNormals
