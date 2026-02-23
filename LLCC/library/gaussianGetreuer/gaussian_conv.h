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
#include <stdbool.h>

#ifndef GAUSSIAN_CONV_H
#define GAUSSIAN_CONV_H

#ifdef __cplusplus
extern "C" {
#endif
    

void Gaussian1D(float *input, float *output, const uint8_t *nodata, int w, float sigma);

void Gaussian2D(float *input, float *output, const uint8_t *nodata, int w, int h, float sigma);

#ifdef __cplusplus
}
#endif


#endif
