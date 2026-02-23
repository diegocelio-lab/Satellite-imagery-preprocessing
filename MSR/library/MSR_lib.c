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

/**
 *  @file MSR_lib.c
 *
 *  @brief Libraries using in the MSR.cpp
 *
 *
 *  @author Catalina Sbert Juan
 */


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "MSR_lib.h"


#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif
#define PI2  6.283185307179586  /* 2*pi*/


/**
 * @brief Convolution with a Gaussian kernel using FFT.
 *
 *
 * @param input double array
 * @param scale the size  of the gaussian kernel
 * @param nx x-size of the array
 * @param ny y-size of the array
 *
 * @return output the convolved array
 */

double *convolution(double *input, double scale, double *output,
                    size_t nx, size_t ny)
{
    double *out;
    fftw_plan p;
    int image_size, image_size4;
    int i,j,index;
    double sigma,normx, normy;

    out = (double*) fftw_malloc(sizeof(double) * (nx*ny));

    /*compute the Fourier transform of the input data*/

    p= fftw_plan_r2r_2d((int)ny, (int)nx, input, out, FFTW_REDFT10,
                        FFTW_REDFT10,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    /*define the gaussian constants for the convolution*/

    sigma=scale*scale/2.;
    normx=M_PI/(double)nx;
    normy=M_PI/(double) ny;
    normx*=normx;
    normy*=normy;

    image_size=(int)nx * (int)ny;
    image_size4=4*image_size;

    for(j=0; j<(int)ny; j++)
    {
        index=j*(int)nx;
        for(i=0; i<(int)nx; i++)
            out[i+index]*=exp((double)(-sigma)*(normx*i*i+normy*j*j));
    }

    /*compute the Inverse Fourier transform */

    p=fftw_plan_r2r_2d((int)ny, (int)nx, out, output, FFTW_REDFT01,
                       FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(p);

    for(index=0; index<image_size; index++)
        output[index]/=image_size4;

    fftw_destroy_plan(p);
    fftw_free(out);

    return output;
}

/**
 * @brief Convolution with a Gaussian kernel using FFT, ignoring invalid pixels.
 *
 * This function performs a 2D convolution of the input image with a Gaussian kernel
 * of the given scale, taking into account a nodata mask. Pixels marked as invalid
 * in the mask are ignored during the convolution, and the result is properly
 * normalized near edges to avoid bias from missing pixels.
 *
 * @param input double array containing the input image
 * @param scale the standard deviation of the Gaussian kernel
 * @param output double array where the convolved result will be stored
 * @param nx x-size of the image
 * @param ny y-size of the image
 * @param nodata uint8_t array of size nx*ny; 1 for invalid pixels, 0 for valid
 *
 * @return output the convolved array, normalized for nodata pixels
 */
double *convolution_masked(double *input, uint8_t *nodata, double scale, double *output,
                            size_t nx, size_t ny)
{
    size_t image_size = nx * ny;
    double *tmp = (double*) malloc(image_size * sizeof(double));
    double *mask = (double*) malloc(image_size * sizeof(double));
    double *conv_mask = (double*) malloc(image_size * sizeof(double));

    // Prepare masked input and mask
    for(size_t i = 0; i < image_size; i++) {
        tmp[i]  = nodata[i] ? 0.0 : input[i]; // zero out nodata pixels
        mask[i] = nodata[i] ? 0.0 : 1.0;      // 1 for valid, 0 for invalid
    }

    // Convolve image
    convolution(tmp, scale, output, nx, ny);

    // Convolve mask
    convolution(mask, scale, conv_mask, nx, ny);

    // Normalize convolution near edges
    for(size_t i = 0; i < image_size; i++) {
        if(conv_mask[i] > 1e-12)  // avoid division by zero
            output[i] /= conv_mask[i];
        else
            output[i] = 0.0;      // if all neighbors invalid, set to 0
    }

    free(tmp);
    free(mask);
    free(conv_mask);

    return output;
}

/**
 * @brief The main part of the Multiscale Retinex
 *
 * @f$ MSRout= \sum w (\log(input)-\log(input* G_\sigma))\f$
 *
 * @param input input color channel
 * @param scale[nscales] the  scales for the convolution
 * @param nscales number of scales
 * @param w the weight for each scale
 * @param nx x-size of the image
 * @param ny y-size of the image
 *
 * @return out output of the multiscale retinex
 */

double *MSRetinex(double *out, double *input, uint8_t *nodata, double *scale, int nscales,
                  double w, size_t nx, size_t ny)
{
    int i, image_size, n;
    double *pas;

    image_size=(int)nx* (int)ny;

    pas=(double*) malloc(image_size*sizeof(double));

    /* initialization of the output*/
    
    for(i=0; i<image_size; i++)
        out[i]=0.;

    /* Compute Retinex output*/

    for(n=0; n<nscales; n++)
    {
        convolution_masked(input, nodata, scale[n], pas, nx, ny);
        for(i=0; i<(int)image_size; i++) {
            if(nodata[i])continue;

            if(input[i] <= 0.) input[i] = 1e-3;
            if(pas[i] <= 0.) pas[i] = 1e-3;
            out[i]+=w*(log(input[i])-log(pas[i]));
        }
    }

    free(pas);

    return out;
}