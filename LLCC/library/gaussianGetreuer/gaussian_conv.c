/**
 * \file gaussian_demo.c
 * \brief Gaussian convolution demo
 * \author Pascal Getreuer <getreuer@cmla.ens-cachan.fr>
 *
 * Copyright (c) 2012-2013, Pascal Getreuer
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under, at your option, the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version, or the terms of the
 * simplified BSD license.
 *
 * You should have received a copy of these licenses along this program.
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include <ctype.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gaussian_conv_sii.h"


void Gaussian1D(float *input, float *output, const uint8_t *nodata, int width, float sigma)
{
    /* Stacked integral images. */
    float *buffer = NULL;
    sii_coeffs c;
    
    int K=3; /* default value */
    sii_precomp(&c, sigma, K);
    
    buffer = (float *)malloc(sizeof(float) * sii_buffer_size(c,width));
    
    sii_gaussian_conv(c, output, buffer, input, nodata, width, 1);

    free(buffer);
    
    
}

void Gaussian2D(float *input, float *output, const uint8_t *nodata, int width, int height, float sigma)
{
    /* Stacked integral images. */
    float *buffer = NULL;
    sii_coeffs c;
    
    int K=3; /* default value */
    sii_precomp(&c, sigma, K);
    
    buffer = (float *)malloc(sizeof(float) * sii_buffer_size(c,((width >= height) ? width : height)));
    
    sii_gaussian_conv_image(c, output, buffer, input, nodata, width, height, 1);
    free(buffer);
    
}


