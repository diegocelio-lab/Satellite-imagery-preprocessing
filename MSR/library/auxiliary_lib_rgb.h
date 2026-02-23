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

 #ifndef _AUXILIARY_LIB_CPP_H
#define _AUXILIARY_LIB_CPP_H

#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <type_traits>

/* --- Convert interleaved input (8-bit or 16-bit) to separate RGB doubles --- */
template <typename T>
void input_rgb(const T* input, double* R, double* G, double* B, size_t size)
{
    static_assert(std::is_same<T, uint8_t>::value || std::is_same<T, uint16_t>::value,
                  "T must be uint8_t or uint16_t");
    
    for (size_t n = 0; n < size; n++) {
        R[n] = static_cast<double>(input[n]);
        G[n] = static_cast<double>(input[size + n]);
        B[n] = static_cast<double>(input[2*size + n]);
    }
}

/* --- Convert separate RGB doubles to interleaved output (8-bit or 16-bit) --- */
template <typename T>
void rgb_output(const uint16_t* Rout, const uint16_t* Gout, const uint16_t* Bout, const uint8_t* nodata, T* output, size_t size)
{
    static_assert(std::is_same<T, uint8_t>::value || std::is_same<T, uint16_t>::value,
                  "T must be uint8_t or uint16_t");

    for (size_t n = 0; n < size; n++) {
        if constexpr (std::is_same<T, uint16_t>::value) {
            if(nodata[n]){
                output[n] = 0;
                output[size + n] = 0;
                output[2*size + n] = 0;
                continue;
            }
            output[n] = std::clamp<uint16_t>(Rout[n], 1, 65535);
            output[size + n] = std::clamp<uint16_t>(Gout[n], 1, 65535);
            output[2*size + n] = std::clamp<uint16_t>(Bout[n], 1, 65535);
        } else {
            if(nodata[n]){
                output[n] = 0;
                output[size + n] = 0;
                output[2*size + n] = 0;
                continue;
            }
            output[n] = std::clamp<uint16_t>(Rout[n], 1, 255);
            output[size + n] = std::clamp<uint16_t>(Gout[n], 1, 255);
            output[2*size + n] = std::clamp<uint16_t>(Bout[n], 1, 255);
        }
    }
}

#endif /* _AUXILIARY_LIB_CPP_H */



