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


#include "auxiliary_lib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>


#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) < (b) ? (b):(a))


/**
 *  @brief compute the R G B components of the output image from its gray level
 *
 *   Given a color image C=(R, G, B), given its gray level
 *
 * @f$ gray= (R+ G+ B)/3 \f$
 *
 * Given a modified gray image gray1
 *
 * This function computes an output color image C1=(R1,G1,B1) which each channel is proportinal
 * to the input channel and whose gray level is gray1,
 *
 * @f$ R1=\frac{gray1}{gray} R    G1=\frac{gray1}{gray} G    B1= \frac{gray1}{gray} B \f$
 *
 * Note that we make a restriction and does not permit a factor  greater than 3
 *
 *
 * @param R  red channel of the input color image
 * @param G  green channel of the input color image
 * @param B  blue channel of the input color image
 * @param gray gray level of the input color image
 * @param gray1 modified gray image
 * @param dim size of the image
 *
 * @return R_out new red channel
 * @return G_out new green channel
 * @return B_out new blue channel
 */

void compute_color_from_grayscale(uint16_t *Rout, uint16_t *Gout, uint16_t *Bout,
                                    const double *R, const double *G, const double *B, 
                                    double *gray, double *gray1, size_t image_size, uint16_t bitsPerSample)
{
    double  factor, maxv;
    double depth = (bitsPerSample == 16) ? 65535.0f : 255.0f;

    for(size_t index=0; index < image_size; index++)
    {
        if(gray[index] <= 1.) gray[index]=1.;
        factor = gray1[index] / gray[index];
        if( factor > 3.) factor=3.;
        double Rtmp = factor * R[index];
        double Gtmp = factor * G[index];
        double Btmp = factor * B[index];

        maxv = fmax(Rtmp, fmax(Gtmp, Btmp));

        if (maxv > depth && maxv > 0.0) {
            double max_in = fmax(R[index], fmax(G[index], B[index]));
            factor = depth / max_in;
            Rtmp = factor * R[index];
            Gtmp = factor * G[index];
            Btmp = factor * B[index];
        }
        
        Rout[index] = (uint16_t)round(fmin(depth, Rtmp));
        Gout[index] = (uint16_t)round(fmin(depth, Gtmp));
        Bout[index] = (uint16_t)round(fmin(depth, Btmp));
    }

}

/**
 * @brief  comparison function
 * given x and y  pointers to doubles.
 * Returns  -1 if x < y
 * 0 if x == y
 * +1 if x > y

*/
int myComparisonFunction(const void *x, const void *y)
{
    double dx, dy;

    dx = *(double *)x;
    dy = *(double *)y;

    if (dx < dy) return -1;
    else if (dx > dy) return +1;

    return 0;
}

int myComparisonFunction_uint16(const void *x, const void *y)
{
    uint16_t dx = *(const uint16_t *)x;
    uint16_t dy = *(const uint16_t *)y;

    if (dx < dy) return -1;
    else if (dx > dy) return +1;
    return 0;
}


/**
 * @brief Simplest color balance
 *
 * Sort the pixels values.
 * Compute the minimium as the sorted array at position dim x s1/100
 * Compute the maximum as the sorted array at position N(1-s2/100)-1.
 * Saturate the pixels according to the computed minimum and maximum.
 * Affine transformation between [minimum, maximum] and [0,255]
 *
 * @param data input array
 * @param s1 the percentage of saturated pixels on the left
 * @param s2 the percentage of saturated pixels on the right
 * @param dim size of the array
 *
 * @return data_out the scaled array.
 */
void simplest_color_balance_double(double *data_out, double *data_in, uint8_t *nodata,
                               size_t image_size, float s1, float s2, int bitdepth)
{
    double min, max, scale;
    double depth;

    // Determine depth
    depth = (bitdepth == 16) ? 65535.0f : 255.0f;

    // Count valid pixels
    size_t valid_count = 0;
    for(size_t index = 0; index < image_size; index++)
        if(!nodata[index]) valid_count++;

    // Copy only valid pixels to temporary array for sorting
    double *sortdata = (double*) malloc(valid_count * sizeof(double));
    size_t j = 0;
    for(size_t index = 0; index < image_size; index++) {
        if(!nodata[index]) sortdata[j++] = data_in[index];
    }

    // Sort valid pixels
    qsort(sortdata, valid_count, sizeof sortdata[0], &myComparisonFunction);

    int per1 = (int)(s1 * valid_count / 100.0f);
    int per2 = (int)(s2 * valid_count / 100.0f);

    min = sortdata[per1];
    max = sortdata[valid_count - 1 - per2];

    free(sortdata);


    //fprintf(stderr, "gray stretch min=%.3f, max=%.3f\n", min, max);

    if(max <= min)
        for(size_t index=0; index < image_size; index++)
            data_out[index]=max;
    else
    {
        scale=depth/(max-min);
        for(size_t index=0; index < image_size; index++)
        {
            if(data_in[index] < min) data_out[index]=0.;
            else if(data_in[index]> max) data_out[index]=depth;
            else data_out[index]=scale*(data_in[index]-min);
        }
    }

}

// sane as above but for uint16_t input
void simplest_color_balance_uint16(uint16_t *data_out, uint16_t *data_in, uint8_t *nodata,
                               size_t image_size, float s1,float s2, int bitdepth)
{
    double min, max, scale;
    double depth;

    // Determine depth
    depth = (bitdepth == 16) ? 65535.0f : 255.0f;

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
    qsort(sortdata, valid_count, sizeof sortdata[0], &myComparisonFunction_uint16);

    int per1 = (int)(s1 * valid_count / 100.0);
    int per2 = (int)(s2 * valid_count / 100.0);

    min = sortdata[per1];
    max = sortdata[valid_count - 1 - per2];

    free(sortdata);

    scale = (min != max) ? depth/(max-min) : 0.0f;


    for(size_t index=0; index < image_size; index++) {
        
        if(nodata[index]) {
            data_out[index] = 0;
            continue;
        }

        double v = (data_in[index]-min) * scale; // Add +1 to make sure only no data are 0
        if (v < 1.0)
            data_out[index] = 1;
        else if (v > depth)
            data_out[index] = (uint16_t)depth;
        else
            data_out[index] = (uint16_t)v;
    }

}

/**
 * @brief Computes the gray intensity value of a color image
*
* @f$ gray= (R+ G+ B)/3 \f$
*
* @param data_in input color image
* @param dim size of the image
*
* @return gray output gray image
*
*/

double *gray_intensity(double *gray, const double *R, const double *G, const double *B,
                       size_t image_size)
{
    for(size_t index=0; index < image_size; index++) {
        gray[index]=(double)(R[index]+G[index]+B[index])/3.;
    }

    return gray;

}


