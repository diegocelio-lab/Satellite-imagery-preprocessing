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

#include "midway.h"

#include <math.h>
#include <limits.h>
#include <iostream>
#include "utils.h"

using namespace std;

int Launch(
    const vector<vector<float>>  &imIn, 
    const vector<ImageSize>      &imSize,
    const size_t                  histogram_size,
    const float                   bitdepth,
    const bool                    verbose,
    vector<vector<float>>         *imOut,
    const vector<vector<bool>>   &img_nodata
    )
{
  //! Check if the images have the same number of channels
  size_t nchannels = imSize[0].nChannels;
  for (size_t i = 0; i < imIn.size(); i++) {
    if (nchannels != imSize[i].nChannels)
    {
      cerr << "Different number of channel between the images" << endl;
      return EXIT_FAILURE;
    }
  }

  //! Apply midway equalization over each channel independantly
  if (verbose)
  {
    cout << "Apply Midway equalization...";
  }
  
  vector<const float*> dataIn_ptr(imIn.size());
  vector<float*> dataOut_ptr(imIn.size());

  for (size_t channel = 0; channel < nchannels; channel++)
  {
    for (size_t i = 0; i < imIn.size(); i++) {
      dataIn_ptr[i]  = imIn[i].data() + channel*imSize[i].wh;
      dataOut_ptr[i] = (*imOut)[i].data() + channel*imSize[i].wh;
    }

    if (ComputeMidwayEqualization(dataIn_ptr, imSize,
                                  0.f, bitdepth, histogram_size,  // here set min to 1.0f?
                                  dataOut_ptr, img_nodata) == EXIT_FAILURE)
    {
      return EXIT_FAILURE;
    }
  }
  
  if (verbose)
  {
    cout << "done." << endl;
  }
  
  return EXIT_SUCCESS;
}

int ComputeMidwayEqualization(
    vector<const float*>       &dataIn_ptr,
    const vector<ImageSize>    &imSize,
    const float                 minimum_value,
    const float                 maximum_value,
    const size_t                histogram_size,
    vector<float*>             &dataOut_ptr,
    const vector<vector<bool>> &img_nodata
    )
{
  //! Check the number of bins of the histogram is not 0
  if(!histogram_size)
  {
    cerr << "Histogram size must be at least equal to 1" << endl;
    return EXIT_FAILURE;
  }

  //! Compute the cumulative histograms
  vector<vector<float>> cumulative_histograms(dataIn_ptr.size());
  for (size_t i = 0; i < dataIn_ptr.size(); i++) {
    cumulative_histograms[i] = ComputeNormalizedCumulativeHistogram(dataIn_ptr[i], imSize[i].wh,
       minimum_value, maximum_value, histogram_size, img_nodata[i]);
  }

  //! Compute the contrast change function using a lookup table
  vector<vector<size_t>> f_nn = ComputeContrastChangeFunction(cumulative_histograms);

  //! Apply contrast changes to the data
  if(maximum_value == minimum_value) {
    cerr << "Invalid value range: maximum_value == minimum_value" << endl;
    return EXIT_FAILURE;
  }
  float scale = (histogram_size-1) / (maximum_value - minimum_value);


  for (size_t i = 0; i < dataIn_ptr.size(); i++) {
    for (size_t pos = 0; pos < imSize[i].wh; pos++)
    {
      if (img_nodata[i][pos]) continue;
      size_t histogram_bin = static_cast<size_t> (Crop(scale*(dataIn_ptr[i][pos] - minimum_value), 0, histogram_size-1)); 
      float  x             = static_cast<float> (f_nn[i][histogram_bin]) / (histogram_size);
      dataOut_ptr[i][pos]  = Crop(x*(maximum_value-minimum_value) + minimum_value, 1.0f, maximum_value);  //crop to 1 so that 0 is only for nodata pixels
    }
  }

  return EXIT_SUCCESS;
}

vector<vector<size_t>> ComputeContrastChangeFunction(
    const vector<vector<float>> &histograms
    )
{
  size_t num_histograms = histograms.size();
  if(num_histograms == 0) return {};

  size_t N = histograms[0].size(); // number of bins per histogram (same for all)

  // Result: a vector of mapping tables, one per histogram
  vector<vector<size_t>> contrast_functions(num_histograms, vector<size_t>(N));

  // Iterate over each histogram to compute its mapping
  for(size_t h = 0; h < num_histograms; ++h)
  {
    // Iterate over each bin in the current histogram
    for(size_t pos = 0; pos < N; ++pos)
    {
      float sum_positions = 0.0f;
      
      // Compare the current bin value against all histograms
      for(const auto &hist : histograms)
      {
        size_t pos2;

        // Find the first bin in hist where the value >= current bin value
        for(pos2 = 0; pos2 < hist.size() && hist[pos2] < histograms[h][pos]; ++pos2);

        sum_positions += static_cast<float>(pos2);
      }

      // Store the average position across all histograms
      contrast_functions[h][pos] = static_cast<size_t>(sum_positions / num_histograms);

    }
  }

  return contrast_functions;
}
