/**
 * Copyright (c) 2026, Diego Celio <diego.celio@me.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file MSR_original_flex.cpp
 * @brief Multiscale Retinex with color restoration (flexible 8-bit/16-bit PNG/TIFF)
 *
 * Works with both PNG and TIFF images, automatically preserving the input bit depth
 * and ignoring nodata value pixels (typical in GeoTiff files)
 * Uses auxiliary_lib for input/output conversions.
 *
 * Author: Adapted from original IPOL code
 */

#include <tiffio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <filesystem>

#include "library/io_png.h"
#include "library/io_tiff.h"
#include "library/auxiliary_lib.h"
#include "library/auxiliary_lib_rgb.h"
#include "library/parser.h"

namespace fs = std::filesystem;
using namespace std;

int main(int argc, char **argv)
{
    size_t nx = 0, ny = 0, nc = 0;
    size_t image_size = 0;

    unsigned char *data_in_u8 = nullptr, *data_out_u8 = nullptr;
    uint16_t *data_in_u16 = nullptr, *data_out_u16 = nullptr;

    bool isTIF = false;
    uint16_t bitsPerSample = 0, outbitsPerSample = 0;
    float s1 = 0.0f, s2 = 0.0f;

    vector<double> nodata_values;

    // Command-line options
    vector<OptStruct *> options;
    OptStruct oL  = {"L:", 0, "0.01", NULL, "percentage saturation left"}; options.push_back(&oL);
    OptStruct oR  = {"R:", 0, "0.01", NULL, "percentage saturation right"}; options.push_back(&oR);
    OptStruct oB  = {"B:", 0, "0", NULL, "output image bit depth"}; options.push_back(&oB);

    vector<ParStruct *> pparameters;
    ParStruct pinput  = {"input",  NULL, "input image file (PNG or TIF)"}; pparameters.push_back(&pinput);
    ParStruct poutput = {"output_folder", NULL, "output folder path"}; pparameters.push_back(&poutput);

    if (!parsecmdline("Global_stretch",
                      "Simple global stretch",
                      argc, argv, options, pparameters))
        return EXIT_FAILURE;

    s1 = atof(oL.value);
    s2 = atof(oR.value);
    outbitsPerSample = static_cast<uint16_t>(atoi(oB.value));
    if (outbitsPerSample != 0 && outbitsPerSample != 8 && outbitsPerSample != 16) {
        fprintf(stderr, "Error: output bit depth must be 8 or 16\n");
        return EXIT_FAILURE;
    }

    fs::path input_file = pinput.value;
    fs::path output_folder = poutput.value;

    // Check input/output paths
    if (!fs::exists(input_file)) {
        cerr << "Error: input file does not exist: " << input_file << endl;
        return EXIT_FAILURE;
    }

    if (!fs::exists(output_folder)) {
        fs::create_directories(output_folder);
    }

    // --- Build output filename automatically ---
    string base_name = input_file.stem().string();
    string ext = input_file.extension().string();
    transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    fs::path output_file = output_folder / (base_name + "_SCB" + ext);

    // Detect file type and read input
    string ext_lc = ext;
    if (ext_lc == ".tif" || ext_lc == ".tiff") {
        isTIF = true;

        TIFF *tif = TIFFOpen(input_file.string().c_str(), "r");
        if (!tif) {
            cerr << "Error: cannot open TIFF file: " << input_file << endl;
            return EXIT_FAILURE;
        }

        if (!TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsPerSample)) {
            cerr << "Error: could not read BitsPerSample in " << input_file << endl;
            TIFFClose(tif);
            return EXIT_FAILURE;
        }
        TIFFClose(tif);

        if (bitsPerSample == 8) {
            data_in_u8 = io_tiff_read<unsigned char>(input_file.string().c_str(), &nx, &ny, &nc, &nodata_values);
        } else if (bitsPerSample == 16) {
            data_in_u16 = io_tiff_read<uint16_t>(input_file.string().c_str(), &nx, &ny, &nc, &nodata_values);
        } else {
            cerr << "Error: unsupported TIFF bit depth (" << bitsPerSample << ") in " << input_file << endl;
            return EXIT_FAILURE;
        }
    } else {
        bitsPerSample = 8;
        data_in_u8 = io_png_read_u8(input_file.string().c_str(), &nx, &ny, &nc);
        if (!data_in_u8) {
            cerr << "Error: failed to read PNG: " << input_file << endl;
            return EXIT_FAILURE;
        }
        nodata_values.assign(nc, numeric_limits<double>::quiet_NaN());
    }

    if (outbitsPerSample == 0)
        outbitsPerSample = bitsPerSample;

    image_size = nx * ny;

    // Allocate RGB and nodata buffers
    vector<uint16_t> R(image_size), G(image_size), B(image_size);
    std::vector<uint8_t> nodataR(image_size, 0), nodataG(image_size, 0), nodataB(image_size, 0);

    if (outbitsPerSample == 16)
        data_out_u16 = (uint16_t *)malloc(3 * image_size * sizeof(uint16_t));
    else
        data_out_u8 = (unsigned char *)malloc(3 * image_size * sizeof(unsigned char));

    // Convert input to RGB arrays
    if (bitsPerSample == 16)
        input_rgb(data_in_u16, R.data(), G.data(), B.data(), image_size);
    else
        input_rgb(data_in_u8, R.data(), G.data(), B.data(), image_size);

    // Compute nodata mask
    for (size_t n = 0; n < image_size; ++n) {
        uint8_t isNoData = (R[n] == 0 && G[n] == 0 && B[n] == 0) ? 1 : 0;
        nodataR[n] = (std::isnan(nodata_values[0])) ? isNoData : (R[n] == nodata_values[0]) ? 1 : 0;
        nodataG[n] = (std::isnan(nodata_values[1])) ? isNoData : (G[n] == nodata_values[1]) ? 1 : 0;
        nodataB[n] = (std::isnan(nodata_values[2])) ? isNoData : (B[n] == nodata_values[2]) ? 1 : 0;
    }

    // Apply color balance stretch
    simplest_color_balance(R.data(), R.data(), nodataR.data(), image_size, s1, s2, outbitsPerSample);
    simplest_color_balance(G.data(), G.data(), nodataG.data(), image_size, s1, s2, outbitsPerSample);
    simplest_color_balance(B.data(), B.data(), nodataB.data(), image_size, s1, s2, outbitsPerSample);

    // Write output
    if (outbitsPerSample == 16)
        rgb_output(R.data(), G.data(), B.data(), data_out_u16, image_size);
    else
        rgb_output(R.data(), G.data(), B.data(), data_out_u8, image_size);

    if (isTIF) {
        if (outbitsPerSample == 16)
            io_tiff_write<uint16_t>(output_file.string().c_str(), data_out_u16, nx, ny, nc, input_file.string().c_str());
        else
            io_tiff_write<unsigned char>(output_file.string().c_str(), data_out_u8, nx, ny, nc, input_file.string().c_str());
    } else {
        io_png_write_u8(output_file.string().c_str(), data_out_u8, nx, ny, nc);
    }

    cout << "Output written to: " << output_file << endl;
    
    // Cleanup
    if (outbitsPerSample == 16) free(data_out_u16);
    else free(data_out_u8);

    if (bitsPerSample == 16) {
        if (isTIF) delete[] data_in_u16;
        else free(data_in_u16);
    } else {
        if (isTIF) delete[] data_in_u8;
        else free(data_in_u8);
    }

    return EXIT_SUCCESS;
}