/**
 *
 * Copyright (c) 2020, Jose-Luis Lisani, joseluis.lisani@uib.es
 * Copyright (c) 2026, Diego celio diego.celio@me.com
 *
 * This file is a modified version of the original work.
 * Modifications made by Diego Celio on 2026:
 *  - Adapted the code to support 8-bit and 16-bit TIFF input images with nodata pixels.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <tiffio.h>
#include "library/io_tiff.h"
#include "library/io_RGB.h"

#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <limits>

//Use Marco Mondelli' implementation of MCM filter
//Available at: http://www.ipol.im/pub/art/2013/53/
#include "library/MCMMondelli/FDS_MCM.h"

//Use Pascal Getreuer' implementation of Gaussian filter
//Available at: http://www.ipol.im/pub/art/2013/87/
#include "library/gaussianGetreuer/gaussian_conv.h"

//Use Sylvain Paris and Frédo Durand fast implementation of bilateral filter
//Available at: http://people.csail.mit.edu/sparis/bf/
#include <cstdint>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include "library/TRUNCATED_KERNEL_BF/include/geom.h"
#include "library/TRUNCATED_KERNEL_BF/include/fast_lbf.h"

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

namespace fs = std::filesystem;
typedef Array_2D<float> image_type;

// -----------------------------
// Linear stretch
// -----------------------------
void stretch_range(vector<float>& A, const vector<uint8_t>& nodata, int w, int h)
{
    float minVal = 255.0f;
    float maxVal = 0.0f;

    size_t size = w * h;
    for (size_t n = 0; n < size; n++) {
        if (!nodata[n]) {
            minVal = std::min(minVal, A[n]);
            maxVal = std::max(maxVal, A[n]);
        }
    }

    if (maxVal == minVal) return;

    for (size_t n = 0; n < size; n++) {
        if (!nodata[n]) A[n] = 255.0f * (A[n] - minVal) / (maxVal - minVal);
    }
}

// -----------------------------
// Read input image
// -----------------------------
void read_input_image(const char* name,
                      vector<float>& R,
                      vector<float>& G,
                      vector<float>& B,
                      vector<float>& I,
                      vector<uint8_t>& nodata,
                      int& w, int& h, int& c,
                      float bitdepth)
{
    size_t ww, hh, cc;
    std::vector<double> nodata_values;
    uint16_t* img = nullptr;

    img = io_tiff_read<uint16_t>(name, &ww, &hh, &cc, &nodata_values);
    c = static_cast<int>(cc);


    if (!img) {
        throw std::runtime_error(string("Error reading image ") + name);
    }

    w = static_cast<int>(ww);
    h = static_cast<int>(hh);
    size_t size = w * h;

    R.resize(size);
    G.resize(size);
    B.resize(size);
    I.resize(size);
    nodata.resize(size);

    for (size_t n = 0; n < size; n++) {
        if(std::isnan(nodata_values[0]) && 
           std::isnan(nodata_values[1]) && 
           std::isnan(nodata_values[2])) {
            nodata[n] = (img[n] == 0 && 
                         img[size + n] == 0 && 
                         img[2 * size + n] == 0) ? 1 : 0;
        } else {
            nodata[n] = (img[n] == nodata_values[0] || 
                         img[size + n] == nodata_values[1] || 
                         img[2 * size + n] == nodata_values[2]) ? 1 : 0;
        }
    }

    input2RGB(img, R.data(), G.data(), B.data(), size, bitdepth);
    RGBtoI(R.data(), G.data(), B.data(), I.data(), size);

    delete[] img;
}



// -----------------------------
// Local contrast correction
// -----------------------------
//option == 1 -> Gaussian filter
//option == 2 -> MCM filter
//option == 3 -> Bilateral filter
void loglocal_correction(vector<float>& A,
                         const vector<uint8_t>& nodata,
                         int w, int h,
                         float sigmaI, float sigmaS, float gradth,
                         float gammalog, int option)
{
    stretch_range(A, nodata, w, h);

    size_t size = w * h;
    vector<float> temp(size);
    vector<float> mask(size);

    // Normalize intensity
    for (size_t n = 0; n < size; n++) temp[n] = A[n] / 255.0f;

    // Filter selection
    if (option == 1) {
        Gaussian2D(temp.data(), mask.data(), nodata.data(), w, h, sigmaS);
    } else if (option == 2) {
        mcm_main(temp.data(), mask.data(), w, h, sigmaS, gradth / 255.0f);
    } else { // option == 3, bilateral
        image_type image(w, h);
        for (size_t n = 0; n < size; n++) image(n % w, n / w) = temp[n];

        image_type filtered_image(w, h);
        Image_filter::fast_LBF(image, image, sigmaS, sigmaI / 255.0f, false, &filtered_image, &filtered_image);

        for (size_t n = 0; n < size; n++) mask[n] = filtered_image(n % w, n / w);
    }

    // Compute alpha map
    for (size_t n = 0; n < size; n++) {
        if (mask[n] <= 0.5f)
            temp[n] = 0.5f - 0.5f * std::pow(mask[n] / 0.5f, gammalog);
        else
            temp[n] = -(0.5f - 0.5f * std::pow((1.0f - mask[n]) / 0.5f, gammalog));
    }

    // Apply tone mapping
    for (size_t n = 0; n < size; n++) {
        if (nodata[n]) continue;

        float alpha = temp[n];
        float aagmax = 255.0f / std::log(255.0f * std::fabs(alpha) + 1.0f);
        float aag;

        if (alpha > 0) aag = aagmax * std::log(alpha * A[n] + 1.0f);
        else if (alpha < 0) aag = 255.0f - aagmax * std::log(-alpha * (255.0f - A[n]) + 1.0f);
        else aag = A[n];

        A[n] = std::clamp(aag, 0.0f, 255.0f);
    }
}


// -----------------------------
// Color processing
// -----------------------------
void color_processing(vector<float>& R, vector<float>& G, vector<float>& B,
                      const vector<float>& I, const vector<float>& Iout,
                      const vector<uint8_t>& nodata, int w, int h)
{
    size_t total = w * h;
    float min_r = 255.0f, min_g = 255.0f, min_b = 255.0f;

    for (size_t n = 0; n < total; n++) {
        if (!nodata[n] && I[n] != 0.0f) {
            float factorI = Iout[n] / I[n];

            float rr = R[n] * factorI;
            float gg = G[n] * factorI;
            float bb = B[n] * factorI;

            float outmax = std::max({rr, gg, bb});
            if (outmax > 255.0f) {
                rr *= 255.0f / outmax;
                gg *= 255.0f / outmax;
                bb *= 255.0f / outmax;
            }

            R[n] = rr; G[n] = gg; B[n] = bb;

            min_r = std::min(min_r, R[n]);
            min_g = std::min(min_g, G[n]);
            min_b = std::min(min_b, B[n]);
        }
    }

    for (size_t n = 0; n < total; n++) {
        if (nodata[n] || I[n] == 0.0f) {
            R[n] = min_r; G[n] = min_g; B[n] = min_b;
        }
    }
}


//Main function: implements Algorithm 1
int main(int argc, const char **argv)
{

    if (argc < 3) {
        printf("Usage: llcc input.tif output_file_path [output bit depth=16] [option=1] [sigmaS=5] [sigmaI=70/gradth=10]\n");
        printf("       option == 1 --> Use Gaussian filter\n");
        printf("       option == 2 --> Use MCM filter\n");
        printf("       option == 3 --> Use Bilinear filter\n");
        printf("       Parameter of Gaussian filter: sigmaS (spatial scale)\n");
        printf("       Parameter of MCM filter: sigmaS (spatial scale)\n");
        printf("                              : gradth (gradient threshold)\n");
        printf("       Parameter of Bilinear filter: sigmaI (intensity scale)\n");
        printf("                                     sigmaS (spatial scale)\n");
        return EXIT_FAILURE;
    }

    // Input and output paths
    fs::path input_file = argv[1];
    fs::path output_folder = argv[2];

    if (!fs::exists(input_file)) {
        cerr << "Error: input file does not exist: " << input_file << endl;
        return EXIT_FAILURE;
    }

    // Create output folder if it doesn't exist
    if (!fs::exists(output_folder)) {
        fs::create_directories(output_folder);
    }

    // File stem and extension
    string base_name = input_file.stem().string();
    string ext = input_file.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    fs::path output_file = output_folder / (base_name + "_LLCC" + ext);

    // Parameters
    int bitdepth_arg = (argc > 3) ? std::atoi(argv[3]) : 16;
    float outbitdepth = (bitdepth_arg == 8) ? 255.0f : 65535.0f;

    int option = (argc > 4) ? std::atoi(argv[4]) : 1;
    if (option < 1 || option > 3) option = 1;

    float sigmaS = (argc > 5) ? std::atof(argv[5]) : 5.0f;
    float sigmaI = 70.0f;
    float gradth = 10.0f;

    if (argc > 6) {
        if (option == 2) gradth = std::atof(argv[6]);   // MCM
        if (option == 3) sigmaI = std::atof(argv[6]);   // Bilateral
    }

    float gammalog = 0.05f;

    // Determine input image type
    if(!(ext == ".tif" || ext == ".tiff")){
        cerr << "Error: input file is not .tif " << input_file << endl;
        return EXIT_FAILURE;
    }

    float input_bitdepth = 255.0f;

    
    TIFF* tif = TIFFOpen(input_file.string().c_str(), "r");
    if (!tif) {
        cerr << "Error: cannot open TIFF file " << input_file << endl;
        return EXIT_FAILURE;
    }
    uint16_t bitsPerSample = 0;
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    input_bitdepth = (bitsPerSample == 8) ? 255.0f : 65535.0f;
    TIFFClose(tif);
    

    // Read input image
    vector<float> R, G, B, I;
    vector<uint8_t> nodata;
    int w, h, c;

    try {
        read_input_image(input_file.string().c_str(), R, G, B, I, nodata, w, h, c, input_bitdepth);
    } catch (const std::exception& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }

    // Apply LLCC to intensity channel
    vector<float> Iout = I;
    loglocal_correction(Iout, nodata, w, h, sigmaI, sigmaS, gradth, gammalog, option);

    // Recover chrominance
    color_processing(R, G, B, I, Iout, nodata, w, h);

    // Save result
    if (outbitdepth == 255.0f) {
        vector<unsigned char> out(w * h * c);
        RGB2output<unsigned char>(R.data(), G.data(), B.data(), nodata.data(), out.data(), w * h, outbitdepth);

        io_tiff_write<unsigned char>(output_file.string().c_str(), out.data(), w, h, c, input_file.string().c_str());
    } else {
        vector<uint16_t> out(w * h * c);
        RGB2output<uint16_t>(R.data(), G.data(), B.data(), nodata.data(), out.data(), w * h, outbitdepth);

        io_tiff_write<uint16_t>(output_file.string().c_str(), out.data(), w, h, c, input_file.string().c_str());
    }

    cout << "Processed: " << input_file.filename() << " -> " << output_file << endl;

    return EXIT_SUCCESS;
}
               
               
               
