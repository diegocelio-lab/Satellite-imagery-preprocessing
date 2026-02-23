/*
 * Copyright (c) 2015, Thierry Guillemot <thierry.guillemot.work@gmail.com>
 * All rights reserved.
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

#include <iostream>
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdint>
#include <tiffio.h>
#include <cmath>

#include "io_tiff.h"
#include "Midway/midway.h"
#include "LibImages.h"

using namespace std;
namespace fs = filesystem;

/**
 * @brief print usage function
 */
void PrintUsage()
{
  cerr << "usage : input_img_folder output_img_folder [options]" << endl
       << "options :" << endl
       << "-dithering sigma : add a gaussian noise of variance sigma before computing midway" << endl
       << "-verbose : activate verbose mode" << endl;
}

/**
 * @brief main function call
 */
int main(int argc, char *const *argv)
{
  vector<vector<float>>   im;           
  vector<string>          img_names;
  vector<vector<bool>>    img_nodata;
  std::vector<double>     nodata_values;
  vector<ImageSize>       imSize; 
  string                  in_folder_path;
  string                  out_folder_path;   
  float                   sigma_dithering = .0f;
  size_t                  image_dynamic   = 8;
  bool                    verbose         = false; 
  std::vector<bool>       isTIF;
  uint16_t                bitsPerSample;
  std::vector<float>      checkbitdepth;
  float                   bitdepth        = 255.0f;
  
  //! Parsing and analysing the function arguments
  for (int pos = 3; pos < argc; ++pos)
  {
    if (strcmp(argv[pos], "-dithering") == 0)
    {
      if (++pos == argc)
      {
        cerr << "Wrong dithering value" << endl;
        PrintUsage();
        return EXIT_FAILURE;
      }
      sigma_dithering = atof(argv[pos]);
    }
    else if (strcmp(argv[pos], "-verbose") == 0)
    {
      cout << "Verbose mode enabled\n";
      verbose = true;
    }
    else
    {
      PrintUsage();
      return EXIT_FAILURE;
    }
  }

  //! Read and load all input images

  // All allowed image extensions
  vector<string> exts = {".png", ".tif", ".tiff"};

  // Go through provided image folder
  in_folder_path = argv[1];
  out_folder_path = argv[2];
  try {
    for (const auto& entry : fs::directory_iterator(in_folder_path)) {
      if (entry.is_regular_file()) {
        string ext = entry.path().extension().string();
        transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        // Check if extension is handled
        if (find(exts.begin(), exts.end(), ext) != exts.end()) {
          // Save image name
          string img_name = entry.path().string();
          img_names.push_back(img_name);

          // Initialise img variable
          vector<float> im_n;
          ImageSize imSize_n;

          // Handle TIF images
          if (ext == ".tif" || ext == ".tiff") {
            isTIF.push_back(true);

            // Get bitdepth
            TIFF* tif = TIFFOpen(img_name.c_str(), "r");
            if (!tif) {
                fprintf(stderr, "Failed to open TIFF file %s\n", img_name.c_str());
                return EXIT_FAILURE;
            }
            if (!TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsPerSample)) {
                fprintf(stderr, "Could not read BitsPerSample tag in %s\n", img_name.c_str());
                TIFFClose(tif);
                return EXIT_FAILURE;
            }
            TIFFClose(tif);
            checkbitdepth.push_back((bitsPerSample == 8) ? 255.0f: 65535.0f);

            // Load
            if (loadImage(img_name.c_str(), im_n, imSize_n, verbose, true, &nodata_values) != EXIT_SUCCESS) {
              cerr << "Failed reading image " << img_name << endl;
              return EXIT_FAILURE;
            }
          } 
          // Handle png images
          else {
            isTIF.push_back(false);
            checkbitdepth.push_back(255.0f);

            // Load
            if (loadImage(img_name.c_str(), im_n, imSize_n, verbose, false, &nodata_values) != EXIT_SUCCESS) {
              cerr << "Failed reading image " << img_name << endl;
              return EXIT_FAILURE;
            }
          }

          // Save in image vector
          im.push_back(im_n);
          imSize.push_back(imSize_n);

        }
      }
    }
  } catch (const fs::filesystem_error& e) {
    cerr << "Error while loading images :" << e.what() << endl;
  }

  // I HAVE TO HANDLE MULTIPLE DIFFERENT INPUT IMAGE BITDEPTH
  bitdepth = checkbitdepth[0];
  image_dynamic = (bitdepth > 255.0f) ? 16 : 8;

  // Compute nodata area of both images
  img_nodata.resize(im.size());
  
  for (size_t i = 0; i < img_nodata.size(); i++) {
    img_nodata[i].resize(imSize[i].whc);
    size_t img_size =  imSize[i].wh;

    for (size_t n = 0; n < img_size; n++) {
      bool isNoData = (im[i][n] == 0.0f && im[i][img_size + n] == 0.0f && im[i][2*img_size + n] == 0.0f);
      img_nodata[i][n] = (isnan(nodata_values[0])) ? isNoData : im[i][n] == nodata_values[0];
      img_nodata[i][img_size + n] = (isnan(nodata_values[1])) ? isNoData : im[i][n] == nodata_values[1];
      img_nodata[i][2*img_size + n] = (isnan(nodata_values[2])) ? isNoData : im[i][n] == nodata_values[2];
    }
  }

  //! Apply dithering
  for (size_t i = 0; i < im.size(); i++) {
    if (sigma_dithering)
    {
      addNoise(im[i], im[i], img_nodata[i], sigma_dithering, bitdepth, verbose);
    }
  }
  
  //! Launch Midway equalization
  if (Launch(im, imSize, 1<<image_dynamic, bitdepth, verbose, &im, img_nodata) != EXIT_SUCCESS)
  {
    cerr << "Failed applying midway equalization" << endl;
    return EXIT_FAILURE;
  }

  //! Stretch image range(?) and save output
  for (size_t i = 0; i < im.size(); i++) {
    //Stretch(im[i], img_nodata[i], imSize[i].wh, 0.01, 0.01, bitdepth);

    if (saveImage(img_names[i], out_folder_path, img_names[i].c_str(), im[i], imSize[i], 1.f, bitdepth, isTIF[i], img_nodata[i]) != EXIT_SUCCESS)
    {
    cerr << "Failed writing output image of " << img_names[i] << endl;
    return EXIT_FAILURE;
    }
  }
  
  return EXIT_SUCCESS;
}
