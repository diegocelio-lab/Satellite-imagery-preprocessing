/*
 *  auxiliary_lib.c
 *
 *  Created by Catalina Sbert Juan on 13/11/12.
 *  Copyright 2012 Universitat de les Illes Balears. All rights reserved.
 *
 *  Modifications by Diego Celio, 2026:
 *    - Added handling for nodata pixels
 *    - Added handling for 16bit inputs
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdint.h>
#ifndef _AUXILIARY_LIB_H
#define _AUXILIARY_LIB_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Compute RGB channels from modified grayscale */
void compute_color_from_grayscale(uint16_t* Rout,
                                  uint16_t* Gout,
                                  uint16_t* Bout,
                                  const double* R,
                                  const double* G,
                                  const double* B,
                                  double* gray,
                                  double* gray1,
                                  size_t image_size,
                                  uint16_t bitsPerSample);

/* Simplest color balance */
void simplest_color_balance_double(double* data_out,
                            double* data_in,
                            uint8_t* nodata,
                            size_t image_size,
                            float s1,
                            float s2,
                            int bitdepth);

void simplest_color_balance_uint16(uint16_t* data_out,
                            uint16_t* data_in,
                            uint8_t* nodata,
                            size_t image_size,
                            float s1,
                            float s2,
                            int bitdepth);

/* Compute gray intensity from RGB channels */
double* gray_intensity(double* gray,
                       const double* R,
                       const double* G,
                       const double* B,
                       size_t image_size);

/* Comparison function for qsort */
int myComparisonFunction(const void* x, const void* y);

int myComparisonFunction_uint16(const void *x, const void *y);

#ifdef __cplusplus
}
#endif

#endif /* _AUXILIARY_LIB_H */
