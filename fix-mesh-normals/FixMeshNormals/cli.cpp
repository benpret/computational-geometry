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
		} else if (!arg.empty() && arg[0] == '-') {
			throw std::runtime_error("Unknown argument: " + arg);
		} else {
			opts.inputs.push_back(std::filesystem::path(arg));
		}
	}

	return opts;
}

} // namespace FixMeshNormals
