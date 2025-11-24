#pragma once

#include <filesystem>
#include <optional>
#include <vector>

namespace FixMeshNormals {

struct Options {
	std::vector<std::filesystem::path> inputs;
	std::optional<std::filesystem::path> out;
	bool overwrite = true;  // Overwrite existing _voxelised meshes by default
};

void PrintUsage();
Options ParseArgs(int argc, char** argv);

} // namespace FixMeshNormals
