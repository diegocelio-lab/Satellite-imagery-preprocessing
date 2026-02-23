 /**
 * 
 * Copyright (c) 2013, Catalina Sbert Juan, IPOL Image Processing On Line
 * Copyright (c) 2026, Diego Celio, diego.celio@me.com
 *
 * This file is a modified version of the original work.
 * Modifications made by Diego Celio on 2026:
 *  - Adapted the code to support 8-bit and 16-bit TIFF input images with nodata pixels.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of:
 *   - the GNU General Public License as published by the Free Software Foundation,
 *     either version 3 of the License, or (at your option) any later version, or
 *   - the simplified BSD license.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received copies of these licenses along with this program.
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */


#include "auxiliary_lib.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>


#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) < (b) ? (b):(a))



/**
 * @brief  comparison function
 * given x and y  pointers to doubles.
 * Returns  -1 if x < y
 * 0 if x == y
 * +1 if x > y

*/
int myComparisonFunction(const void *x, const void *y)
{
    uint16_t dx = *(const uint16_t *)x;
    uint16_t dy = *(const uint16_t *)y;

    if (dx < dy) return -1;
    else if (dx > dy) return +1;
    return 0;
}



// sane as above but for uint16_t input
void simplest_color_balance(uint16_t *data_out, uint16_t *data_in, uint8_t *nodata,
                               size_t image_size, float s1,float s2, int bitdepth)
{
    double min, max, scale;
    double depth;

    // Determine depth
    depth = (bitdepth == 16) ? 65535.0 : 255.0;

    // Count valid pixels
    size_t valid_count = 0;
    for(size_t index = 0; index < image_size; index++)
        if(!nodata[index]) valid_count++;

    // Copy only valid pixels to temporary array for sorting
    uint16_t *sortdata = (uint16_t*) malloc(valid_count * sizeof(uint16_t));
    size_t j = 0;
    for(size_t index = 0; index < image_size; index++) {
        if(!nodata[index]) sortdata[j++] = data_in[index];
    }

    // Sort valid pixels
    qsort(sortdata, valid_count, sizeof sortdata[0], &myComparisonFunction);

    int per1 = (int)(s1 * valid_count / 100.0);
    int per2 = (int)(s2 * valid_count / 100.0);

    min = sortdata[per1];
    max = sortdata[valid_count - 1 - per2];

    free(sortdata);

    if(max <= min)
        for(size_t index=0; index < image_size; index++)
            data_out[index]=max;
    else
    {
        scale=depth/(max-min);
        for(size_t index=0; index < image_size; index++)
        {
            if(nodata[index]) continue;
            if(data_in[index] <= min) data_out[index]=0. + 1.; //Add +1 to make sure that only nodata pixel are 0
            else if(data_in[index]> max) data_out[index]=depth;
            else data_out[index]=scale*(data_in[index]-min);
        }
    }

}



