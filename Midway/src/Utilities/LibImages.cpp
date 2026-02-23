/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Copyright (c) 2015, Thierry Guillemot <thierry.guillemot.work@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file LibImages.cpp
 * @brief Usefull functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibImages.h"

#include <iostream>
#include <filesystem>
#include <sstream>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm> 

#include "io_png.h"
#include "io_tiff.h"
#include "mt19937ar.h"
#include "utils.h"


using namespace std;
namespace fs = std::filesystem;

int loadImage(
    const char*          p_name,
    vector<float>       &o_im,
    ImageSize           &o_imSize,
    const bool           p_verbose,
    const bool           isTIF,
    std::vector<double>* nodata_values
    )
{
  //! read input image
  if (p_verbose)
  {
    cout << endl << "Read input image '" << p_name << "'...";
  }

  float *imTmp = NULL;
  size_t w, h, c;

  if (isTIF)
  {
    imTmp = io_tiff_read<float>(p_name, &w, &h, &c, nodata_values);
  } else {
    imTmp = read_png_f32(p_name, &w, &h, &c);
    nodata_values->assign(c, std::numeric_limits<double>::quiet_NaN());

  }
  
  if (!imTmp)
  {
    cout << "error :: " << p_name << " not found or not a correct png image" << endl;
    return EXIT_FAILURE;
  }
  
  if (p_verbose)
  {
    cout << "done." << endl;
  }

  //! test if image is really a color image and exclude the alpha channel
  if (c > 2)
  {
    unsigned k = 0;
    while (k < w * h && imTmp[k] == imTmp[w * h + k] && imTmp[k] == imTmp[2 * w * h + k])
    {
      k++;
    }
    c = (k == w * h ? 1 : 3);
  }

  //! Some image informations
  if (p_verbose)
  {
    cout << "image size :" << endl;
    cout << " - width          = " << w << endl;
    cout << " - height         = " << h << endl;
    cout << " - nb of channels = " << c << endl;
  }

  //! Initializations
  o_imSize.width      = w;
  o_imSize.height     = h;
  o_imSize.nChannels  = c;
  o_imSize.wh         = w * h;
  o_imSize.whc        = w * h * c;
  o_im.resize(w * h * c);
  for (unsigned k = 0; k < w * h * c; k++)
    o_im[k] = imTmp[k];

  //! Free Memory
  delete[] imTmp;

  return EXIT_SUCCESS;
}

int saveImage(
    const std::string        &filename,
    const std::string        &outname,
    const char*               fileref,
    const vector<float>      &i_im,
    const ImageSize          &p_imSize,
    const float               p_min,
    const float               p_max,
    const bool                isTIF,
    const vector<bool>       &nodata
    )
{
  // Create the output filename
  fs::path out_folder = outname;
  fs::path p = filename;

  if (!fs::exists(out_folder)) {
    fs::create_directories(out_folder);
  }

  std::string base = p.stem().string();  // "image"
  std::string ext = p.extension().string(); // ".png"

  std::string new_name = base + "_Midway" + ext;
  fs::path out_path = out_folder / new_name;
    
  //! Allocate Memory                         ---------- THIS IS FOR 16BIT
  uint16_t* imTmp = new uint16_t[p_imSize.whc];

  //! Check for boundary problems
  for (unsigned k = 0; k < p_imSize.whc; k++)
  {
    if (nodata[k]) {
      imTmp[k] = i_im[k];
      continue;
    }
    imTmp[k] = static_cast<uint16_t>(i_im[k] < p_min ? p_min : (i_im[k] > p_max ? p_max : i_im[k]));
  }

  if (isTIF) {
    if (io_tiff_write<uint16_t>(out_path.string().c_str(), imTmp, p_imSize.width, p_imSize.height, p_imSize.nChannels, fileref) != 0)
    {
      cout << "... failed to save tif image :'" << new_name << "'" << endl;
      return EXIT_FAILURE;
    }
  } else {
    /*
    if (write_png_f32(p_name.c_str(), imTmp, p_imSize.width, p_imSize.height, p_imSize.nChannels) != 0)
    {
      cout << "... failed to save png image :'" << p_name << "'" << endl;
      return EXIT_FAILURE;
    }
      */
  }

  //! Free Memory
  delete[] imTmp;

  return EXIT_SUCCESS;
}

void addNoise(
    const vector<float> &i_im,
    vector<float>       &o_imNoisy,
    const vector<bool>  &nodata,
    const float          p_sigma,
    const float          bitdepth,
    const bool           p_verbose
    )
{
  if (p_verbose)
  {
    cout << "Add noise [sigma = " << p_sigma << "] ...";
  }

  //! Initialization
  o_imNoisy = i_im;
  mt_init_genrand((unsigned long int) time (NULL) + (unsigned long int) getpid());

  //! Add noise
  for (unsigned k = 0; k < i_im.size(); ++k)
  {
    if (nodata[k]) continue;

    const double a = mt_genrand_res53();
    const double b = mt_genrand_res53();

    o_imNoisy[k] += p_sigma * (float) (sqrtl(-2.0l * log(a)) * cos(2.0l * M_PI * b));
    o_imNoisy[k] = Crop(o_imNoisy[k], 0.f, bitdepth); 
  }

  if (p_verbose)
  {
    cout << "done." << endl;
  }
}

int computePsnr(
    const vector<float> &i_im1,
    const vector<float> &i_im2,
    float               &o_psnr,
    float               &o_rmse,
    const float          bitdepth,
    const char*          p_imageName,
    const bool           p_verbose
    )
{
  if (i_im1.size() != i_im2.size())
  {
    cout << "Can't compute PSNR & RMSE: images have different sizes: " << endl;
    cout << "i_im1 : " << i_im1.size() << endl;
    cout << "i_im2 : " << i_im2.size() << endl;
    return EXIT_FAILURE;
  }

  float sum = 0.f;
  for (unsigned k = 0; k < i_im1.size(); k++)
    sum += (i_im1[k] - i_im2[k]) * (i_im1[k] - i_im2[k]);

  o_rmse = sqrtf(sum / (float) i_im1.size());
  o_psnr = 20.f * log10f(bitdepth / o_rmse);

  if (p_verbose)
  {
    cout << p_imageName << endl;
    cout << "PSNR = " << o_psnr << endl;
    cout << "RMSE = " << o_rmse << endl;
  }

  return EXIT_SUCCESS;
}

vector<float> ComputeNormalizedCumulativeHistogram(
    const float         *data_ptr,
    const size_t        size,
    const float         min,
    const float         max,
    const size_t        histogram_size,
    const vector<bool>  &nodata
    )
{
  //! Initialization of the histogram with 0 
  vector<float> histogram(histogram_size, 0);
  size_t size_max = histogram_size - 1; 
  
  //! Computation of the histogram
  for(size_t pos=0; pos<size; ++pos)
  {
    if (nodata[pos]) continue;
    ++histogram[(size_t) ( size_max*(data_ptr[pos]- min)/(max - min))];
  }
  
  //! Normalization of the cumulative histogram
  int num_valid = count(nodata.begin(), nodata.begin()+size, false);
  for(size_t pos=0; pos<histogram_size-1; ++pos)
  {
    histogram[pos+1] += histogram[pos];
    histogram[pos]   /= num_valid;
  }
  histogram[histogram_size-1] /= num_valid;

  return histogram;
}

vector<float> ComputeNormalizedHistogram(
    const float  *data_ptr,
    const size_t  size,
    const float   min,
    const float   max,
    const size_t  histogram_size
    )
{
  //! Initialization of the histogram with 0 
  vector<float> histogram(histogram_size, 0);
  size_t size_max = histogram_size - 1; 
  
  //! Computation of the histogram
  for(size_t pos=0; pos<size; ++pos)
  {
    ++histogram[(size_t) ( size_max*(data_ptr[pos]- min)/(max - min))];
  }
  
  //! Normalization of the cumulative histogram
  for(size_t pos=0; pos<histogram_size; ++pos)
  {
    histogram[pos]   /= size;
  }

  return histogram;

  
}

int myComparisonFunction(const void *x, const void *y)
{
    float dx, dy;

    dx = *(float *)x;
    dy = *(float *)y;

    if (dx < dy) return -1;
    else if (dx > dy) return +1;

    return 0;
}

void Stretch(
    vector<float>   &data, 
    vector<bool>    &nodata,
    size_t          image_size, 
    float           s1,
    float           s2,  
    float           bitdepth
  )
{
    // Count valid pixels
    size_t num_valid = count(nodata.begin(), nodata.begin()+image_size, false);

    auto get_thresholds = [data, num_valid, image_size, nodata, s1, s2](int channel, float &min_val, float &max_val)
    {
        // Copy only valid pixels to temporary array for sorting
        float *sortdata = (float*) malloc(num_valid * sizeof(float));
        size_t j = 0;
        for(size_t index = channel; index < channel + image_size; index++) {
            if(!nodata[index]) sortdata[j++] = data[index];
        }

        // Sort valid pixels
        qsort(sortdata, num_valid, sizeof sortdata[0], &myComparisonFunction);

        int per1 = (int)(s1 * num_valid / 100.0f);
        int per2 = (int)(s2 * num_valid / 100.0f);

        min_val = sortdata[per1];
        max_val = sortdata[num_valid - 1 - per2];

        free(sortdata);
    };
    

    float r_min, r_max, g_min, g_max, b_min, b_max;
    get_thresholds(0, r_min, r_max);
    get_thresholds(image_size, g_min, g_max);
    get_thresholds(2*image_size, b_min, b_max);

    float r_scale = (r_max != r_min) ? (bitdepth / (r_max - r_min)) : 0.0f;
    float g_scale = (g_max != g_min) ? (bitdepth / (g_max - g_min)) : 0.0f;
    float b_scale = (b_max != b_min) ? (bitdepth / (b_max - b_min)) : 0.0f;

    for(size_t n = 0; n < image_size; n++)
    {
      float r = (data[n] - r_min) * r_scale;
      float g = (data[image_size + n] - g_min) * g_scale;
      float b = (data[2 * image_size + n] - b_min) * b_scale;

      // Clamp to safe range
      data[n] = std::max(0.0f, std::min(r, bitdepth));
      data[image_size + n] = std::max(0.0f, std::min(g, bitdepth));
      data[2 * image_size + n] = std::max(0.0f, std::min(b, bitdepth));

    }

}
