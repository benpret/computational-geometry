#include "cli.h"
#include "usd_io.h"

#include <openvdb/openvdb.h>

#include <filesystem>
#include <iostream>

int main(int argc, char** argv) {
	try {
		openvdb::initialize();

		FixMeshNormals::Options opts = FixMeshNormals::ParseArgs(argc, argv);
		if (opts.inputs.empty()) {
			FixMeshNormals::PrintUsage();
			return 1;
		}
		if (opts.out && opts.inputs.size() != 1) {
			std::cerr << "--out can only be used with a single input file.\n";
			return 1;
		}

		int exitCode = 0;
		for (const auto& rawInputPath : opts.inputs) {
			std::filesystem::path inputPath = rawInputPath;
			if (inputPath.is_relative()) {
				inputPath = std::filesystem::absolute(inputPath);
			}

			if (!std::filesystem::exists(inputPath)) {
				std::cerr << "Input file does not exist: " << inputPath << "\n";
				exitCode = 1;
				continue;
			}

			std::filesystem::path outputPath = opts.out ? *opts.out : inputPath;
			if (outputPath.is_relative()) {
				outputPath = std::filesystem::absolute(outputPath);
			}

			std::size_t created = FixMeshNormals::ProcessStageFile(
				inputPath,
				outputPath,
				opts.overwrite,
				opts.processElements,
				opts.writeVoxelised,
				opts.weldTolerance);
			if (created == 0) {
				std::cout << "  No meshes created for " << inputPath << "\n";
			}
		}

		std::cout << "Parsed " << opts.inputs.size() << " input file(s).\n";
		if (opts.out) {
			std::cout << "Writing output to: " << opts.out->string() << "\n";
		}
		return exitCode;
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	} catch (...) {
		std::cerr << "Unknown error encountered.\n";
		return 1;
	}
}
