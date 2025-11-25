#pragma once

#include <filesystem>
#include <optional>
#include <vector>

namespace FixMeshNormals {

struct Options {
	std::vector<std::filesystem::path> inputs;
	std::optional<std::filesystem::path> out;
	bool overwrite = true;  // Overwrite existing _voxelised meshes by default
	bool processElements = false; // Process each connected element separately
	bool writeVoxelised = true;   // Write _voxelised outputs by default
	double weldTolerance = 0.001; // Default weld tolerance in stage units (1 mm)
};

void PrintUsage();
Options ParseArgs(int argc, char** argv);

} // namespace FixMeshNormals
