/**
 * Copyright (c) 2026, Diego Celio, diego.celio@me.com
 *
 * This file is part of the Seamline project.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License version 3 for more details.
 *
 * You should have received a copy of the GNU General Public License
 * version 3 along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

 #include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <cstdlib>

#include "seamline.h"
#include "io_tiff.h"


using std::vector;
namespace fs = std::filesystem;

int main(int argc, char **argv)
{
    vector<std::string> img_names;

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_folder> <output_folder> [visibleSeam=true|false]\n";
        return EXIT_FAILURE;
    }

    fs::path input_folder = argv[1];
    fs::path output_folder = argv[2];

    bool visibleSeam = false;
    if (argc > 3) {
        std::string arg = argv[3];
        std::transform(arg.begin(), arg.end(), arg.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        visibleSeam = (arg == "true" || arg == "1" || arg == "yes");
    }

    std::cout << "Input folder:  " << input_folder << "\n";
    std::cout << "Output folder: " << output_folder << "\n";
    std::cout << "Visible seam:  " << (visibleSeam ? "true" : "false") << "\n";

    if (!fs::exists(input_folder)) {
        std::cerr << "Error: input file does not exist: " << input_folder << std::endl;
        return EXIT_FAILURE;
    }

    // Create output folder if it doesn't exist
    if (!fs::exists(output_folder)) {
        fs::create_directories(output_folder);
    }

    // All allowed image extensions  -- for now only tif
    vector<std::string> exts = {".tif", ".tiff"};
    for (const auto& entry : fs::directory_iterator(input_folder)) {
        if (entry.is_regular_file()) {
            std::string ext = entry.path().extension().string();
            transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
                
            // Check if extension is handled
            if (find(exts.begin(), exts.end(), ext) != exts.end()) {   
                // Construct new filename with _Seamline
                fs::path filename = entry.path().stem();
                filename += "_Seamline";
                filename += entry.path().extension();

                fs::path dst = output_folder / filename;

                try {
                    fs::copy_file(entry.path(), dst, fs::copy_options::overwrite_existing);
                } catch (fs::filesystem_error& e){
                    std::cerr << "Failed to copy " << entry.path() << " to " 
                    << dst << ": " << e.what() << std::endl;
                }

                img_names.push_back(dst.string());
            } 
        }
    }

    for (size_t i = 0; i < img_names.size(); ++i) {
        for (size_t j = i + 1; j < img_names.size(); ++j) {
            std::cout << "Start semaline computation for : " << img_names[i] << " and " << img_names[j] << "\n";
            if (Seamline(img_names[i], img_names[j], i, j, visibleSeam) != 0) {
                std::cerr << "Seamlines not computed ! \n";
            }
        }
    }

    return 0;
}
