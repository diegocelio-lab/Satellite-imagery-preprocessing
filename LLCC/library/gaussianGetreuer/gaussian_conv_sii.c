/**
 * \file gaussian_conv_sii.c
 * \brief Gaussian convolution using stacked integral images
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
 * You should have received a copy of these licenses along with this program.
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include "gaussian_conv_sii.h"
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>


#ifndef M_PI
/** \brief The constant pi */
#define M_PI        3.14159265358979323846264338327950288
#endif

/**
 * \brief Half-sample symmetric boundary extension
 * \param N     signal length
 * \param n     requested sample, possibly outside {0,...,`N`-1}
 * \return reflected sample in {0,...,`N`-1}
 *
 * This function is used for boundary handling. Suppose that `src` is an array
 * of length `N`, then `src[extension(N, n)]` evaluates the symmetric
 * extension of `src` at location `n`.
 *
 * Half-sample symmetric extension is implemented by the pseudocode
 \verbatim
 repeat
 if n < 0, reflect n over -1/2
 if n >= N, reflect n over N - 1/2
 until 0 <= n < N
 \endverbatim
 * The loop is necessary as some `n` require multiple reflections to bring
 * them into the domain {0,...,`N`-1}.
 *
 * This function is used by all of the Gaussian convolution algorithms
 * included in this work except for DCT-based convolution (where symmetric
 * boundary handling is performed implicitly by the transform). For FIR, box,
 * extended box, SII, and Deriche filtering, this function could be replaced
 * to apply some other boundary extension (e.g., periodic or constant
 * extrapolation) without any further changes. However, Alvarez-Mazorra and
 * Vliet-Young-Verbeek are hard-coded for symmetric extension on the right
 * boundary, and would require specific modification to change the handling
 * on the right boundary.
 *
 * \par A note on efficiency
 * This function is a computational bottleneck, as it is used within core
 * filtering loops. As a small optimization, we encourage inlining by defining
 * the function as `static`. We refrain from further optimization since this
 * is a pedagogical implementation, and code readability is more important.
 * Ideally, filtering routines should take advantage of algorithm-specific
 * properties such as exploiting sequential sample locations (to update the
 * extension cheaply) and samples that are provably in the interior (where
 * boundary checks may omitted be entirely).
 */
static long extension(long N, long n)
{
    while (1)
        if (n < 0)
            n = -1 - n;         /* Reflect over n = -1/2.    */
        else if (n >= N)
            n = 2 * N - 1 - n;  /* Reflect over n = N - 1/2. */
        else
            break;
    
    return n;
}


/**
 * \brief Precompute filter coefficients for SII Gaussian convolution
 * \param c         sii_coeffs pointer to hold precomputed coefficients
 * \param sigma     Gaussian standard deviation
 * \param K         number of boxes = 3, 4, or 5
 * \return 1 on success, 0 on failure
 *
 * This routine reads Elboher and Werman's optimal SII radii and weights for
 * reference standard deviation \f$ \sigma_0 = 100/\pi \f$ from a table and
 * scales them to the specified value of sigma.
 */
void sii_precomp(sii_coeffs *c, double sigma, int K)
{
    /* Elboher and Werman's optimal radii and weights. */
    const double sigma0 = 100.0 / M_PI;
    static const short radii0[SII_MAX_K - SII_MIN_K + 1][SII_MAX_K] =
        {{76, 46, 23, 0, 0},
         {82, 56, 37, 19, 0},
         {85, 61, 44, 30, 16}};
    static const float weights0[SII_MAX_K - SII_MIN_K + 1][SII_MAX_K] =
        {{0.1618f, 0.5502f, 0.9495f, 0, 0},
         {0.0976f, 0.3376f, 0.6700f, 0.9649f, 0},
         {0.0739f, 0.2534f, 0.5031f, 0.7596f, 0.9738f}};

    const int i = K - SII_MIN_K;
    double sum;
    int k;
    
    assert(c && sigma > 0 && SII_VALID_K(K));
    c->K = K;
    
    for (k = 0, sum = 0; k < K; ++k)
    {
        c->radii[k] = (long)(radii0[i][k] * (sigma / sigma0) + 0.5);
        sum += weights0[i][k] * (2 * c->radii[k] + 1);
    }
    
    for (k = 0; k < K; ++k)
        c->weights[k] = (float)(weights0[i][k] / sum);
    
    return;
}

/**
 * \brief Determines the buffer size needed for SII Gaussian convolution
 * \param c     sii_coeffs created by sii_precomp()
 * \param N     number of samples
 * \return required buffer size in units of num samples
 *
 * This routine determines the minimum size of the buffer needed for use in
 * sii_gaussian_conv() or sii_gaussian_conv_image(). This size is the length
 * of the signal (or in 2D, max(width, height)) plus the twice largest box
 * radius, for padding.
 */
long sii_buffer_size(sii_coeffs c, long N)
{
    long pad = c.radii[0] + 1;
    return 2 * (N + 2 * pad);
}

/**
 * \brief Gaussian convolution (SII approximation) with nodata mask
 *
 * \param c         Filter coefficients from sii_precomp()
 * \param dest      Output array for convolved data
 * \param buffer    Temporary array (twice sii_buffer_size(c,N) floats)
 * \param src       Input data (may equal dest for in-place use)
 * \param nodata    Mask array; true marks invalid samples
 * \param N         Number of samples
 * \param stride    Stride between successive samples
 *
 * Performs a stacked integral image (SII) Gaussian convolution on a 1-D
 * signal, ignoring samples marked in \a nodata. Invalid samples do not
 * contribute to the cumulative sums; results are renormalized by the total
 * valid weight to preserve intensity near missing data. Boundary handling
 * is half-sample symmetric via extension().
 *
 * Output is set to 0.0f where no valid samples are present.
 */

void sii_gaussian_conv(const sii_coeffs c,
                       float *dest,
                       float *buffer,
                       const float *src,
                       const uint8_t *nodata,
                       long N,
                       long stride)
{
    assert(dest && buffer && src && nodata && N > 0 && stride > 0);

    long pad = c.radii[0] + 1;

    float *buf_val = buffer + pad;                  // size: N + 2*pad
    float *buf_wgt = buf_val + (N + 2*pad);    // size: N + 2*pad

    double accum_val = 0.0;
    double accum_wgt = 0.0;

    /* --- cumulative sums --- */
    for (long n = -pad; n < N + pad; ++n) {
        long e = extension(N, n);
        long idx = stride * e;
        bool valid = !nodata[idx];

        float value = valid ? src[idx] : 0.0f;
        float weight = valid ? 1.0f : 0.0f;

        accum_val += value;
        accum_wgt += weight;

        buf_val[n] = (float)accum_val;
        buf_wgt[n] = (float)accum_wgt;
    }

    /* --- Gaussian via stacked box filters --- */
    for (long n = 0; n < N; ++n, dest += stride) {
        long e = extension(N, n);
        long idx = stride * e;

        if (nodata[idx]) {
            *dest = 0.0f;  // nodata
            continue;
        }

        double val_sum = 0.0;
        double wgt_sum = 0.0;

        for (int k = 0; k < c.K; ++k) {
            long r = c.radii[k];
            double wk = c.weights[k];

            val_sum += wk * (buf_val[n + r] - buf_val[n - r - 1]);
            wgt_sum += wk * (buf_wgt[n + r] - buf_wgt[n - r - 1]);
        }

        *dest = (wgt_sum > 0.0) ? (float)(val_sum / wgt_sum) : NAN;
    }
}

/**
 * \brief 2D Gaussian convolution (SII approximation)
 *
 * Performs a 2D Gaussian-like convolution using stacked integral images.
 * Supports multiple channels and nodata masking.
 *
 * \param c             sii_coeffs created by sii_precomp()
 * \param dest          output array (can be same as src for in-place)
 * \param buffer        preallocated workspace (size >= 2*(max(width,height) + 2*(c.radii[0]+1)))
 * \param src           input image array
 * \param nodata        boolean mask of invalid pixels (true = nodata)
 * \param width         image width
 * \param height        image height
 * \param num_channels  number of image channels
 *
 * Notes:
 * - The buffer must be distinct from src/dest.
 * - Pixels marked as nodata in the mask are set to NAN in the output.
 * - Neighboring pixels are normalized by the sum of valid weights, so edges
 *   near nodata regions are computed correctly.
 */
void sii_gaussian_conv_image(const sii_coeffs c,
                             float *dest,
                             float *buffer,
                             const float *src,
                             const uint8_t *nodata,
                             int width,
                             int height,
                             int num_channels)
{
    assert(dest && buffer && src && width > 0 && height > 0 && num_channels > 0);

    long num_pixels = (long)width * height;

    for (int channel = 0; channel < num_channels; ++channel) {
        float *dest_y = dest;
        const float *src_y = src;
        const uint8_t *nodata_y = nodata;

        /* Filter rows */
        for (int y = 0; y < height; ++y) {
            sii_gaussian_conv(c, dest_y, buffer, src_y, nodata_y, width, 1);
            dest_y += width;
            src_y += width;
            nodata_y += width;
        }

        /* Filter columns */
        for (int x = 0; x < width; ++x)
            sii_gaussian_conv(c, dest + x, buffer, dest + x, nodata + x, height, width);

        dest += num_pixels;
        src += num_pixels;
        nodata += num_pixels;
    }
}


