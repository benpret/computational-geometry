#include "cli.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace FixMeshNormals {

void PrintUsage() {
	std::cout << "Usage: FixMeshNormals [options] <input.usd...>\n"
			  << "Options:\n"
			  << "  --out <file>            Output USD (only valid with a single input)\n"
			  << "  --skip-existing         Skip existing _voxelised meshes (default: overwrite)\n"
			  << "  --process-elements      Process each connected element within a mesh\n"
			  << "  --no-process-elements   Disable element processing (process whole mesh)\n"
			  << "  --write-voxelised       Write _voxelised meshes (default: on)\n"
			  << "  --no-write-voxelised    Skip writing _voxelised meshes\n"
			  << "  --weld-tolerance <m>    Vertex weld tolerance in stage units (default: 0.001 = 1mm)\n"
			  << "  -h, --help              Show this help text\n";
}

Options ParseArgs(int argc, char** argv) {
	Options opts;
	for (int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if (arg == "-h" || arg == "--help") {
			PrintUsage();
			std::exit(0);
		} else if (arg == "--out") {
			if (i + 1 >= argc) {
				throw std::runtime_error("--out expects a file path argument");
			}
			opts.out = std::filesystem::path(argv[++i]);
		} else if (arg == "--skip-existing") {
			opts.overwrite = false;
		} else if (arg == "--process-elements") {
			opts.processElements = true;
		} else if (arg == "--no-process-elements") {
			opts.processElements = false;
		} else if (arg == "--write-voxelised") {
			opts.writeVoxelised = true;
		} else if (arg == "--no-write-voxelised") {
			opts.writeVoxelised = false;
		} else if (arg == "--weld-tolerance") {
			if (i + 1 >= argc) {
				throw std::runtime_error("--weld-tolerance expects a numeric argument");
			}
			try {
				opts.weldTolerance = std::stod(argv[++i]);
			} catch (const std::exception&) {
				throw std::runtime_error("Invalid value for --weld-tolerance");
			}
			if (opts.weldTolerance <= 0.0) {
				throw std::runtime_error("--weld-tolerance must be positive");
			}
		} else if (!arg.empty() && arg[0] == '-') {
			throw std::runtime_error("Unknown argument: " + arg);
		} else {
			opts.inputs.push_back(std::filesystem::path(arg));
		}
	}

	return opts;
}

} // namespace FixMeshNormals
