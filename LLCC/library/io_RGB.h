/**
 *
 * Copyright (c) 2020, Jose-Luis Lisani, joseluis.lisani@uib.es
 * Copyright (c) 2026, Diego celio diego.celio@me.com
 *
 * This file is a modified version of the original work.
 * Modifications made by Diego Celio on 2026:
 *  - Adapted the code to support 8-bit and 16-bit TIFF input images with nodata pixels.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef IO_RGB_HEADER
#define IO_RGB_HEADER


#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

// -------------------------------------------------------
// Helper: clamp float to a range
// -------------------------------------------------------
inline float clampf(float val, float min_val, float max_val) {
    return (val < min_val) ? min_val : (val > max_val) ? max_val : val;
}


// -------------------------------------------------------
// Convert RGB channels to intensity channel (simple average)
// -------------------------------------------------------
inline void RGBtoI(const float *R, const float *G, const float *B, float *I, const size_t size)
{
    for (size_t n = 0; n < size; ++n)
    {
        I[n] = (R[n] + G[n] + B[n]) / 3.0f;
    }
}

// -------------------------------------------------------
// Convert planar input buffer to R/G/B channels
// Assumes input layout: R-plane, G-plane, B-plane
// -------------------------------------------------------
inline void input2RGB(const uint16_t *input, float *R, float *G, float *B, const size_t size, const float bitdepth)
{
    for (size_t n = 0; n < size; ++n)
    {
        R[n] = static_cast<float>(input[n]) * 255.0f / bitdepth;
        G[n] = static_cast<float>(input[size + n]) * 255.0f / bitdepth;
        B[n] = static_cast<float>(input[2 * size + n]) * 255.0f / bitdepth;
    }
}

// -------------------------------------------------------
// Rescale RGB channels and write to output
// Template: U = output type (uint8_t, uint16_t, etc.)
// -------------------------------------------------------
template <typename U>
inline void RGB2output(const float *R, const float *G, const float *B, const uint8_t *nodata,
                       U *output, size_t size, float bitdepth)
{
    for (size_t n = 0; n < size; ++n)
    {
        if (nodata[n]) {
            output[n] = output[size + n] = output[2 * size + n] = 0;
            continue;
        }

        float rescale = bitdepth / 255.0f;

        output[n]       = static_cast<U>(clampf(R[n]*rescale, 1.0f, bitdepth));
        output[size+n]  = static_cast<U>(clampf(G[n]*rescale, 1.0f, bitdepth));
        output[2*size+n]= static_cast<U>(clampf(B[n]*rescale, 1.0f, bitdepth));

    }
}

#endif // IO_RGB_HEADER
