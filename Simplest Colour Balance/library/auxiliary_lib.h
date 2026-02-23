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

#include <stdbool.h>
#include <stdint.h>
#ifndef _AUXILIARY_LIB_H
#define _AUXILIARY_LIB_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Simplest color balance */
void simplest_color_balance(uint16_t* data_out,
                            uint16_t* data_in,
                            uint8_t* nodata,
                            size_t image_size,
                            float s1,
                            float s2,
                            int bitdepth);

int myComparisonFunction(const void *x, const void *y);

#ifdef __cplusplus
}
#endif

#endif /* _AUXILIARY_LIB_H */
