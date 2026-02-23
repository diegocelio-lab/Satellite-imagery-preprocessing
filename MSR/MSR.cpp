/*
 *Original code:
 * Copyright 2013 IPOL Image Processing On Line http://www.ipol.im/
 *
 * This file implements an algorithm possibly linked to the patents:
 *
 *  - US 5991456, "Method of improving a digital image," Issued Nov 23, 1999
 *  - US 6834125, "Method of improving a digital image as a function of its
 *  dynamic range," Issued Dec 21, 2004
 *  - US 6842543 B2, "Method of improving a digital image having white
 *  zones," Issued Jan 11, 2005
 *  - US 8111943, "Smart Image Enhancement Process," Issued Feb 7, 2012
 *  - EP 0901671, "Method of improving a digital image,"
 *  Issued September 3, 2003
 *  - AUS 713076, "Method of improving a digital image,"
 *  Issued February 26, 1998
 *  - WO 1997045809 A1, "Method of improving a digital image," July 4, 2006
 *  - JPO 4036391 B2, "Method of improving a digital image"
 *
 * This file is made available for the exclusive aim of serving as
 * scientific tool to verify the soundness and completeness of the
 * algorithm description. Compilation, execution and redistribution of
 * this file may violate patents rights in certain countries. The
 * situation being different for every country and changing
 * over time, it is your responsibility to determine which patent rights
 * restrictions apply to you before you compile, use, modify, or
 * redistribute this file. A patent lawyer is qualified to make this
 * determination. If and only if they don't conflict with any patent
 * terms, you can benefit from the following license terms attached to this
 * file.
 *
 * Modifications by Diego Celio, 2026:
 *  - Added support for TIF 8bit and 16bit inputs
 *  - Added handling for nodata pixels
 *
 * Your modifications are licensed under the GNU General Public License v3.
 * See LICENSE file for full details.
 */

#include <tiffio.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <memory>
#include "library/io_png.h"
#include "library/io_tiff.h"
#include "library/MSR_lib.h"
#include "library/auxiliary_lib.h"
#include "library/auxiliary_lib_rgb.h"
#include "library/parser.h"

namespace fs = std::filesystem;
using namespace std;

int main(int argc, char **argv)
{
    int nscales;
    double w, scale[3];
    size_t nx = 0, ny = 0, nc = 0, image_size = 0;

    // Using smart pointers for input data
    std::unique_ptr<unsigned char[]> data_in_u8;
    std::unique_ptr<uint16_t[]> data_in_u16;
    std::unique_ptr<unsigned char[]> data_outG_u8;
    std::unique_ptr<uint16_t[]> data_outG_u16;

    // Work arrays as std::vector for automatic cleanup
    std::vector<double> R, G, B, gray, grayout;
    std::vector<uint16_t> Rout, Gout, Bout;
    std::vector<uint8_t> nodata;
    std::vector<double> nodata_values;

    float s1, s2;
    bool isTIF = false;
    uint16_t bitsPerSample = 0;
    uint16_t outbitsPerSample = 0;

    fs::path input_file, output_folder;

    // --- Command line setup ---
    std::vector<OptStruct *> options;
    OptStruct oS = {"S:", 0, "3", NULL, "number of scales (1-3)"}; options.push_back(&oS);
    OptStruct oL = {"L:", 0, "15", NULL, "Low scale"}; options.push_back(&oL);
    OptStruct oM = {"M:", 0, "80", NULL, "Medium scale"}; options.push_back(&oM);
    OptStruct oH = {"H:", 0, "250", NULL, "High scale"}; options.push_back(&oH);
    OptStruct ol = {"l:", 0, "0.1", NULL, "percentage saturation left"}; options.push_back(&ol);
    OptStruct oR = {"R:", 0, "0.1", NULL, "percentage saturation right"}; options.push_back(&oR);
    OptStruct oB = {"B:", 0, "0", NULL, "output image bit depth"}; options.push_back(&oB);

    std::vector<ParStruct *> pparameters;
    ParStruct pinput = {"input", NULL, "input image file (PNG or TIFF)"}; pparameters.push_back(&pinput);
    ParStruct poutput = {"output", NULL, "output image folder"}; pparameters.push_back(&poutput);

    if (!parsecmdline("MSR",
                      "Multiscale Retinex with color restoration (flexible 8/16-bit PNG/TIFF)",
                      argc, argv, options, pparameters))
        return EXIT_FAILURE;

    // --- Parse options ---
    nscales = atoi(oS.value);
    if (nscales < 1 || nscales > 3) {
        fprintf(stderr, "The number of scales must be 1, 2, or 3\n");
        return EXIT_FAILURE;
    }

    if (nscales == 1) scale[0] = atof(oL.value);
    else if (nscales == 2) { scale[0] = atof(oL.value); scale[1] = atof(oM.value); }
    else { scale[0] = atof(oL.value); scale[1] = atof(oM.value); scale[2] = atof(oH.value); }

    s1 = atof(ol.value);
    s2 = atof(oR.value);

    outbitsPerSample = atoi(oB.value);
    if (outbitsPerSample != 0 && outbitsPerSample != 8 && outbitsPerSample != 16) {
        fprintf(stderr, "The output image bit depth must be 0, 8, or 16\n");
        return EXIT_FAILURE;
    }

    input_file = pinput.value;
    output_folder = poutput.value;

    // Create output folder if it doesn't exist
    if (!fs::exists(output_folder)) {
        fs::create_directories(output_folder);
    }

    // --- Build output filename ---
    string base_name = input_file.stem().string();
    string ext = input_file.extension().string();
    transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    fs::path output_file = output_folder / (base_name + "_MSR" + ext);

    // --- Read input file (PNG/TIFF) ---
    if (ext == ".tif" || ext == ".tiff") {
        isTIF = true;

        TIFF* tif = TIFFOpen(input_file.string().c_str(), "r");
        if (!tif) {
            fprintf(stderr, "Failed to open TIFF file %s\n", input_file.string().c_str());
            return EXIT_FAILURE;
        }

        if (!TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsPerSample)) {
            fprintf(stderr, "Could not read BitsPerSample tag in %s\n", input_file.string().c_str());
            TIFFClose(tif);
            return EXIT_FAILURE;
        }
        TIFFClose(tif);

        if (bitsPerSample == 8) {
            data_in_u8.reset(io_tiff_read<unsigned char>(input_file.string().c_str(), &nx, &ny, &nc, &nodata_values));
            if (!data_in_u8) return EXIT_FAILURE;
        } else if (bitsPerSample == 16) {
            data_in_u16.reset(io_tiff_read<uint16_t>(input_file.string().c_str(), &nx, &ny, &nc, &nodata_values));
            if (!data_in_u16) return EXIT_FAILURE;
        } else {
            fprintf(stderr, "Unsupported TIFF bit depth (%d)\n", bitsPerSample);
            return EXIT_FAILURE;
        }
    } else {
        bitsPerSample = 8;
        data_in_u8.reset(io_png_read_u8(input_file.string().c_str(), &nx, &ny, &nc));
        nodata_values.assign(nc, std::numeric_limits<double>::quiet_NaN());
        if (!data_in_u8) {
            fprintf(stderr, "Failed to read PNG file %s\n", input_file.string().c_str());
            return EXIT_FAILURE;
        }
    }

    if (outbitsPerSample == 0) outbitsPerSample = bitsPerSample;

    image_size = nx * ny;

    // --- Allocate working arrays ---
    R.resize(image_size);
    G.resize(image_size);
    B.resize(image_size);
    gray.resize(image_size);
    grayout.resize(image_size);
    Rout.resize(image_size);
    Gout.resize(image_size);
    Bout.resize(image_size);
    nodata.resize(image_size);

    if (outbitsPerSample == 16)
        data_outG_u16 = std::make_unique<uint16_t[]>(3 * image_size);
    else
        data_outG_u8 = std::make_unique<unsigned char[]>(3 * image_size);

    // --- Convert input to RGB arrays ---
    if (bitsPerSample == 16)
        input_rgb(data_in_u16.get(), R.data(), G.data(), B.data(), image_size);
    else
        input_rgb(data_in_u8.get(), R.data(), G.data(), B.data(), image_size);

    // --- Compute nodata mask ---
    for (size_t n = 0; n < image_size; n++) {
        if (std::isnan(nodata_values[0]) &&
            std::isnan(nodata_values[1]) &&
            std::isnan(nodata_values[2])) {
            nodata[n] = (R[n] == 0 && G[n] == 0 && B[n] == 0) ? 1 : 0;
        } else {
            nodata[n] = (R[n] == nodata_values[0] ||
                         G[n] == nodata_values[1] ||
                         B[n] == nodata_values[2]) ? 1 : 0;
        }
    }

    // --- Compute gray intensity ---
    gray_intensity(gray.data(), R.data(), G.data(), B.data(), image_size);
    w = 1.0 / nscales;

    // --- Apply Multiscale Retinex ---
    MSRetinex(grayout.data(), gray.data(), nodata.data(), scale, nscales, w, nx, ny);

    // --- Stretch gray output ---
    simplest_color_balance_double(grayout.data(), grayout.data(), nodata.data(), image_size, s1, s2, bitsPerSample);  // 0,1 to remove outliers

    // --- Restore color from gray ---
    compute_color_from_grayscale(Rout.data(), Gout.data(), Bout.data(),
                                 R.data(), G.data(), B.data(),
                                 gray.data(), grayout.data(), image_size, bitsPerSample);

    // --- Flatten RGB for output ---
    if (outbitsPerSample == 16)
        rgb_output(Rout.data(), Gout.data(), Bout.data(), nodata.data(), data_outG_u16.get(), image_size);
    else
        rgb_output(Rout.data(), Gout.data(), Bout.data(), nodata.data(), data_outG_u8.get(), image_size);

    // --- Write output file ---
    if (isTIF) {
        if (outbitsPerSample == 16)
            io_tiff_write<uint16_t>(output_file.string().c_str(), data_outG_u16.get(), nx, ny, nc, input_file.string().c_str());
        else
            io_tiff_write<unsigned char>(output_file.string().c_str(), data_outG_u8.get(), nx, ny, nc, input_file.string().c_str());
    } else {
        io_png_write_u8(output_file.string().c_str(), data_outG_u8.get(), nx, ny, nc);
    }

    return EXIT_SUCCESS;
}
