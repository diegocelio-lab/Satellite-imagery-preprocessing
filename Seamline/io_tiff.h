/**
 * Copyright (c) 2026, Diego Celio, diego.celio@me.com
 *
 * This file is part of the Seamline project.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License version 3 for more details.
 *
 * You should have received a copy of the GNU General Public License
 * version 3 along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

 #ifndef _IO_TIFF_H
#define _IO_TIFF_H

#include <tiffio.h>
#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <queue>
#include <iostream>     // std::cerr, std::cout
#include <vector>       // std::vector
#include <array>        // std::array
#include <string>       // std::string
#include <cmath>        // std::abs, std::floor (optional)
#include <type_traits>  // std::is_same
#include <limits>   // for std::numeric_limits
#include <unordered_set>

// for debugging mask
#include <png.h>


// GDAL headers
#include "gdal_priv.h"  // GDALDataset, GDALRasterBand
#include "cpl_conv.h"   // CPLMalloc, CPLErr


struct Point { int r, c; };
struct PointWorld { double x; double y; }; // x=east, y=north

struct Node {
    Point p;
    double dist;
    bool operator>(const Node& other) const { return dist > other.dist; }
};

// ------------------------------- Utilities ----------------------------------
// 4-connected neighborhood (edge-touching only)
constexpr int dr4[4] = { -1, 1, 0, 0 };
constexpr int dc4[4] = { 0, 0, -1, 1 };

// 8-connected neighborhood (corner-touching too)
constexpr int dr8[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };
constexpr int dc8[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
const float edgeLen8[8] = {1, 1, 1, 1, std::sqrt(2.0f), std::sqrt(2.0f), std::sqrt(2.0f), std::sqrt(2.0f)};

// Flat memory-friendly templates
using Channel = std::vector<float>;
using ImageRGB = std::array<Channel, 3>;

template <typename T>
using Channel_T = std::vector<T>;

template <typename T>
using ImageRGB_T = std::array<Channel_T<T>, 3>;

// ----------------------------------------
// Build a flat binary mask of valid pixels
// mask[idx] = 1 if valid, 0 if nodata
// ----------------------------------------
template <typename T>
inline uint8_t* buildMask(
    const ImageRGB_T<T> &img, 
    const double (&nodata)[3], 
    size_t H, size_t W)
{
    uint8_t* mask = new uint8_t[H * W];
    for (size_t idx = 0; idx < H * W; ++ idx) {
        mask[idx] = 1;
        for (size_t ch = 0; ch < 3; ++ ch) {
            if (img[ch][idx] == static_cast<T>(nodata[ch])) {
                mask[idx] = 0;
                break;
            }
        }
    }
    return mask; // caller must delete[] mask
}

// ----------------------------------------
// Extract perimeter of a flat mask (1D vector)
// ----------------------------------------
inline Point* extractPerimeter(const uint8_t* mask, size_t H, size_t W, size_t & nPts)
{
    nPts = 0;
    Point* perimeter = new Point[H * W];

    for (size_t r = 0; r < H; ++ r) {
        for (size_t c = 0; c < W; ++ c) {
            size_t idx = r * W + c;
            if (mask[idx] == 1) {
                bool edge = false;
                for (size_t k = 0; k < 8; ++ k) {
                    int rr = static_cast<int>(r) + dr8[k];
                    int cc = static_cast<int>(c) + dc8[k];
                    if (rr < 0 || rr >= static_cast<int>(H) || cc < 0 || cc >= static_cast<int>(W) || mask[rr * W + cc] == 0) {
                        edge = true;
                        break;
                    }
                }
                if (edge) perimeter[nPts ++] = {static_cast<int>(r), static_cast<int>(c)};
            }
        }
    }
    return perimeter; // caller must delete[] perimeter
}


// ----------------------------------------
// Convert world coordinates → pixel coordinates
// using GDAL-like geotransform
// gt[0]: top-left x
// gt[1]: w-e pixel resolution
// gt[2]: row rotation (usually 0)
// gt[3]: top-left y
// gt[4]: column rotation (usually 0)
// gt[5]: n-s pixel resolution (usually negative)
// ----------------------------------------
inline Point worldToPixel(const PointWorld &pw, const double gt[6])
{
    double det = gt[1] * gt[5] - gt[2] * gt[4];
    if (std::abs(det) < 1e-12)
        throw std::runtime_error("Invalid geotransform: determinant is zero.");

    // Invert affine transform
    double inv_det = 1.0 / det;

    double dx = pw.x - gt[0];
    double dy = pw.y - gt[3];

    int c = static_cast<int>(std::floor(( gt[5] * dx - gt[2] * dy) * inv_det));
    int r = static_cast<int>(std::floor((-gt[4] * dx + gt[1] * dy) * inv_det));

    return { r, c };
}


// ----------------------------------------
// Convert pixel coordinates → world
// ----------------------------------------
inline PointWorld pixelToWorld(const Point& p, const double gt[6])
{
    double x = gt[0] + p.c * gt[1] + p.r * gt[2];
    double y = gt[3] + p.c * gt[4] + p.r * gt[5];
    return { x, y };
}

// ----------------------------------------
// Converts an array of pixel points to world coords (batch).
// ----------------------------------------
inline PointWorld* toWorldCoords(const Point* pixels, size_t nPts, const double gt[6])
{
    PointWorld* world = new PointWorld[nPts];
    for (size_t i = 0; i < nPts; ++i)
        world[i] = pixelToWorld(pixels[i], gt);

    return world; // caller must delete[] world
}

// ----------------------------------------
// Compute candidate intersection points of two perimeters
// ----------------------------------------
inline PointWorld* computePerimeterIntersections(
    const PointWorld* p1, size_t n1,
    const PointWorld* p2, size_t n2,
    float tol, size_t &nPts)
{
    nPts = 0;
    PointWorld* candidates = new PointWorld[std::min(n1, n2)]; // rough upper bound

    for (size_t i = 0; i < n1; ++ i) {
        for (size_t j = 0; j < n2; ++ j) {
            double dx = p1[i].y - p2[j].y;
            double dy = p1[i].x - p2[j].x;
            if (std::hypot(dx, dy) <= tol/100) {                    // relevant to keep tolerance here? I divided by 100 otherwise not working
                candidates[nPts ++] = { (p1[i].x + p2[j].x) / 2, (p1[i].y + p2[j].y) / 2 };
            }
        }
    }
    return candidates; // caller must delete[] candidates
}

// ----------------------------------------
// Compute union mask of two images
// 0 = pixel not corresponding to anything (in union mask because of rectangle nature of array)
// 1 = valid only in image 1
// 2 = valid only in image 2
// 3 = nodata pixels (for both images)
// 4 = overlapping valid pixels
// ----------------------------------------
template <typename T>
inline std::vector<uint8_t> computeXORMask(
    const ImageRGB_T<T> &I1, const double (&nodata1)[3], size_t H1, size_t W1,
    const double gt1[6],
    const ImageRGB_T<T> &I2, const double (&nodata2)[3], size_t H2, size_t W2,
    const double gt2[6],
    size_t &unionH, size_t &unionW,
    Point &shift1_union, Point &shift2_union)
{
    // --- Compute bounding boxes in world coordinates
    PointWorld tl1 = pixelToWorld({0, 0}, gt1);
    PointWorld br1 = pixelToWorld({static_cast<int>(H1 - 1), static_cast<int>(W1 - 1)}, gt1);
    PointWorld tl2 = pixelToWorld({0, 0}, gt2);
    PointWorld br2 = pixelToWorld({static_cast<int>(H2 - 1), static_cast<int>(W2 - 1)}, gt2);

    // --- Union bounds
    double x_min = std::min(tl1.x, tl2.x);
    double x_max = std::max(br1.x, br2.x);
    double y_min = std::min(br1.y, br2.y);
    double y_max = std::max(tl1.y, tl2.y);

    // --- Union grid dimensions
    double pxW = gt1[1];
    double pxH = std::abs(gt1[5]);
    unionW = static_cast<size_t>(std::ceil((x_max - x_min) / pxW));
    unionH = static_cast<size_t>(std::ceil((y_max - y_min) / pxH));
    
    // --- Union geotransform
    double gt_union[6] = { x_min, pxW, 0, y_max, 0, -pxH };

    // --- Compute shifts for each image (top-left pixel offset)
    Point img1TopLeftInUnion = worldToPixel(pixelToWorld({0, 0}, gt1), gt_union);
    Point img2TopLeftInUnion = worldToPixel(pixelToWorld({0, 0}, gt2), gt_union);

    shift1_union = img1TopLeftInUnion;
    shift2_union = img2TopLeftInUnion;

    // --- Allocate zero-initialized union mask
    std::vector<uint8_t> diffMask(unionH * unionW, 0); // zero-initialized

    // --- Lambda to fill and XOR mask from a given image
    auto xorMask = [&](const ImageRGB_T<T>& I, const double (&nodata)[3],
                       size_t H, size_t W, const double gt[6], uint8_t true_val)
    {
        for (size_t r = 0; r < H; ++r) {
            for (size_t c = 0; c < W; ++c) {
                size_t idx = r * W + c;
                bool valid = true;
                for (size_t ch = 0; ch < 3; ++ch)
                    if (I[ch][idx] == static_cast<T>(nodata[ch])) { valid = false; break; }

                PointWorld pw = pixelToWorld({static_cast<int>(r), static_cast<int>(c)}, gt);
                Point pUnion = worldToPixel(pw, gt_union);

                if (pUnion.r >= 0 && pUnion.r < static_cast<int>(unionH) &&
                    pUnion.c >= 0 && pUnion.c < static_cast<int>(unionW))
                {
                    size_t uidx = pUnion.r * unionW + pUnion.c;

                    if (!valid) {
                        if (diffMask[uidx] == 0 || diffMask[uidx] == 3) diffMask[uidx] = 3;
                        continue;
                    }

                    // --- XOR toggle logic (flip state)
                    diffMask[uidx] = (diffMask[uidx] == 0 || diffMask[uidx] == 3) ? true_val : 4;
                }
            }
        }
    };

    // --- Fill from both images with XOR semantics
    xorMask(I1, nodata1, H1, W1, gt1, 1);
    xorMask(I2, nodata2, H2, W2, gt2, 2);

    return diffMask;
}


// ----------------------------------------
// Compute and store polygons composed of allowed mask values
// ----------------------------------------
inline void computePolygons(
    std::vector<std::vector<Point>>& Polygons,
    const size_t H, const size_t W,
    const std::vector<uint8_t>& diffMask,
    const std::unordered_set<uint8_t>& allowedValues)
{
    std::vector<std::vector<uint8_t>> visited(H, std::vector<uint8_t>(W, 0)); // zero-initialized

    for (size_t r = 0; r < H; r++) {
        for (size_t c = 0; c < W; c++) {
            uint8_t maskVal = diffMask[r * W + c];

            // Skip if pixel not in allowed set or already visited
            if (!allowedValues.count(maskVal) || visited[r][c])
                continue;

            // Start BFS flood-fill
            std::vector<Point> queue;
            queue.push_back({static_cast<int>(r), static_cast<int>(c)});
            visited[r][c] = 1;

            std::vector<Point> poly;
            size_t qStart = 0;

            while (qStart < queue.size()) {
                Point p = queue[qStart++];
                poly.push_back(p);

                // 4-connected neighbors
                for (size_t k = 0; k < 4; k++) {
                    int rr = p.r + dr4[k];
                    int cc = p.c + dc4[k];

                    if (rr >= 0 && rr < static_cast<int>(H) &&
                        cc >= 0 && cc < static_cast<int>(W) &&
                        !visited[rr][cc] &&
                        allowedValues.count(diffMask[rr * W + cc]))
                    {
                        visited[rr][cc] = 1;
                        queue.push_back({rr, cc});
                    }
                }
            }

            // Store polygon if non-empty
            if (!poly.empty())
                Polygons.push_back(std::move(poly));
        }
    }
}


// ----------------------------------------
// Select two optimal endpoints
// Only considers polygons entirely from one image (mask == 1 or mask == 2)
// ----------------------------------------
inline void selectOptimalEndpoints(
    const std::vector<uint8_t>& diffMask, size_t H, size_t W,
    const std::vector<Point>& candidates, size_t nCandidates,
    std::vector<Point>& Vb, std::vector<Point>& Ve)
{
    std::vector<std::vector<Point>> onePolygons;
    std::vector<std::vector<Point>> overlapPolygons;

    const std::unordered_set<uint8_t> oneValue = {1};
    const std::unordered_set<uint8_t> overlapValue = {4};

    // Compute polygons
    computePolygons(overlapPolygons, H, W, diffMask, overlapValue);
    computePolygons(onePolygons, H, W, diffMask, oneValue);

    std::vector<std::vector<int>> touchedOverlap(nCandidates);
    std::vector<std::vector<int>> touchedOne(nCandidates);

    // Identify which polygons/overlaps each candidate touches
    for (size_t i = 0; i < nCandidates; ++i) {
        const Point& cand = candidates[i];

        for (size_t k = 0; k < 8; ++k) {
            int rr = cand.r + dr8[k];
            int cc = cand.c + dc8[k];
            if (rr < 0 || rr >= static_cast<int>(H) || cc < 0 || cc >= static_cast<int>(W))
                continue;

            // --- Touching 1-polygons ---
            for (int m = 0; m < onePolygons.size(); m++) {
                for (const auto& p : onePolygons[m]) {
                    if (p.r == rr && p.c == cc) {
                        if (std::find(touchedOne[i].begin(), touchedOne[i].end(), m) == touchedOne[i].end()) {
                            touchedOne[i].push_back(m);
                        }
                        break;
                    }
                }
            }
        }
        for (size_t k = 0; k < 4; ++k) {
            int rr = cand.r + dr4[k];
            int cc = cand.c + dc4[k];
            if (rr < 0 || rr >= static_cast<int>(H) || cc < 0 || cc >= static_cast<int>(W))
                continue;

            // --- Touching overlap polygons ---
            for (int n = 0; n < overlapPolygons.size(); n++) {
                if (overlapPolygons[n].size() < 9) continue;
                for (const auto& p : overlapPolygons[n]) {
                    if (p.r == rr && p.c == cc) {
                        if (std::find(touchedOverlap[i].begin(), touchedOverlap[i].end(), n) == touchedOverlap[i].end()) {
                            touchedOverlap[i].push_back(n);
                        }
                        break;
                    }
                }
            }
        }
    }

    // Helper: find the two candidates that are farthest apart
    auto pickTwoFarthest = [&](const std::vector<Point>& cands) {
        if (cands.size() < 2) return;

        double maxDist = -1.0;
        size_t idxA = 0, idxB = 0;

        for (size_t i = 0; i < cands.size(); ++i) {
            for (size_t j = i + 1; j < cands.size(); ++j) {
                double d = std::hypot(double(cands[i].r - cands[j].r),
                                    double(cands[i].c - cands[j].c));
                if (d > maxDist) {
                    maxDist = d;
                    idxA = i;
                    idxB = j;
                }
            }
        }

        Vb.push_back(cands[idxA]);
        Ve.push_back(cands[idxB]);
    };

    // For each overlap polygon and 1-polygon, look for candidates touching both and compute best two candidates
    for (int n = 0; n < overlapPolygons.size(); n++) {
        if (overlapPolygons[n].size() < 5) continue; // all overlap pixels are reachable with 4-connectivity
        for (int m = 0; m < onePolygons.size(); m++) {
            std::vector<Point> relevantCands;
            for (size_t i = 0; i < nCandidates; i++) {
                if (std::find(touchedOverlap[i].begin(), touchedOverlap[i].end(), n) != touchedOverlap[i].end() &&
                    !touchedOne[i].empty() && touchedOne[i][0] == m) {
                        relevantCands.push_back(candidates[i]);
                }
            }
            pickTwoFarthest(relevantCands);
        }
    }
}


// -----------------------------------------------------------------------------
// Helper Functions to read a single band (full image + overlap region)
// -----------------------------------------------------------------------------
template <typename T>
bool read_band(GDALDataset* ds,
               int band_idx,
               Channel_T<T>& full,
               Channel_T<T>& overlap,
               int col_start, int row_start,
               size_t nx, size_t ny,
               size_t ovw, size_t ovh,
               double& nodata_val,
               GDALDataType dtype)
{
    GDALRasterBand* band = ds->GetRasterBand(band_idx);
    if (!band) return false;

    // Read full image
    if (band->RasterIO(GF_Read, 0, 0, static_cast<int>(nx), static_cast<int>(ny),
                       full.data(), static_cast<int>(nx), static_cast<int>(ny),
                       dtype, 0, 0) != CE_None)
        return false;

    // Read overlap
    if (band->RasterIO(GF_Read, col_start, row_start,
                       static_cast<int>(ovw), static_cast<int>(ovh),
                       overlap.data(), static_cast<int>(ovw), static_cast<int>(ovh),
                       dtype, 0, 0) != CE_None)
        return false;

    int hasNoData = 0;
    double nd = band->GetNoDataValue(&hasNoData);
    nodata_val = hasNoData ? nd : std::numeric_limits<double>::quiet_NaN();
    return true;
}

// ----------------------------------------
// Helper function to compute raster bounds in world coordinates (rotation-safe)
// ----------------------------------------
inline void computeRasterBounds(size_t H, size_t W, const double gt[6],
                                double& x_min, double& x_max,
                                double& y_min, double& y_max)
{
    PointWorld tl = pixelToWorld({0, 0}, gt);
    PointWorld tr = pixelToWorld({0, static_cast<int>(W - 1)}, gt);
    PointWorld bl = pixelToWorld({static_cast<int>(H - 1), 0}, gt);
    PointWorld br = pixelToWorld({static_cast<int>(H - 1), static_cast<int>(W - 1)}, gt);

    x_min = std::min({tl.x, tr.x, bl.x, br.x});
    x_max = std::max({tl.x, tr.x, bl.x, br.x});
    y_min = std::min({tl.y, tr.y, bl.y, br.y});
    y_max = std::max({tl.y, tr.y, bl.y, br.y});
}


inline void saveMaskAsPNGG(const std::string& filename,
                          const std::vector<uint8_t>& mask,
                          size_t H, size_t W)
{
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp) { perror("fopen"); return; }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop info = png_create_info_struct(png);
    if (setjmp(png_jmpbuf(png))) return;

    png_init_io(png, fp);
    png_set_IHDR(png, info, W, H, 8,
                 PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png, info);

    const uint8_t* ptr = mask.data();
    for (size_t y = 0; y < H; ++y)
        png_write_row(png, (png_bytep)(ptr + y * W));

    png_write_end(png, nullptr);
    fclose(fp);
    png_destroy_write_struct(&png, &info);

    printf("Wrote debug PNG: %s (%lux%lu)\n", filename.c_str(), W, H);
}


// ----------------------------------------
// Read two GeoTIFFs, extract overlap region and compute seam endpoints
// ----------------------------------------
template <typename T>
bool io_tiff_read_overlap(
    const std::string& fname1, const std::string& fname2,
    size_t& nx1, size_t& ny1,
    size_t& nx2, size_t& ny2,
    ImageRGB_T<T>& I1, ImageRGB_T<T>& I2,
    ImageRGB_T<T>& I1_overlap, ImageRGB_T<T>& I2_overlap,
    std::vector<uint8_t>& diffMask, size_t& unionW, size_t& unionH,
    std::vector<Point>& Vb, std::vector<Point>& Ve,
    size_t& overlap_width, size_t& overlap_height,
    Point& shift1, Point& shift1_union,
    Point& shift2, Point& shift2_union,
    std::vector<double>& nodata_values1, std::vector<double>& nodata_values2)
{
    GDALAllRegister();

    GDALDataset* ds1 = static_cast<GDALDataset*>(GDALOpen(fname1.c_str(), GA_ReadOnly));
    GDALDataset* ds2 = static_cast<GDALDataset*>(GDALOpen(fname2.c_str(), GA_ReadOnly));
    if (!ds1 || !ds2) {
        std::cerr << "Error opening datasets\n";
        return false;
    }

    // GeoTransforms
    double gt1[6], gt2[6];
    if (ds1->GetGeoTransform(gt1) != CE_None || ds2->GetGeoTransform(gt2) != CE_None) {
        std::cerr << "Error: Missing or invalid GeoTransform\n";
        GDALClose(ds1); GDALClose(ds2);
        return false;
    }

    // Raster sizes
    nx1 = ds1->GetRasterXSize(); ny1 = ds1->GetRasterYSize();
    nx2 = ds2->GetRasterXSize(); ny2 = ds2->GetRasterYSize();

    if (ds1->GetRasterCount() != 3 || ds2->GetRasterCount() != 3) {
        std::cerr << "Error: One or both images have fewer or more than 3 bands.\n";
        GDALClose(ds1); GDALClose(ds2);
        return false;
    }

    // --- compute raster bounds (rotation-safe) from four corners ---
    double x1_min, x1_max, y1_min, y1_max;
    computeRasterBounds(ny1, nx1, gt1, x1_min, x1_max, y1_min, y1_max);

    double x2_min, x2_max, y2_min, y2_max;
    computeRasterBounds(ny2, nx2, gt2, x2_min, x2_max, y2_min, y2_max);

    // --- overlap extents in world coords ---
    double x_min_overlap = std::max(x1_min, x2_min);
    double x_max_overlap = std::min(x1_max, x2_max);
    double y_min_overlap = std::max(y1_min, y2_min);
    double y_max_overlap = std::min(y1_max, y2_max);

    if (x_min_overlap >= x_max_overlap || y_min_overlap >= y_max_overlap) {
        std::cerr << "No overlap found!\n";
        GDALClose(ds1);
        GDALClose(ds2);
        return false;
    }

    // --- convert overlap corners to pixel coordinates for each image ---
    // world coords passed as (x, y)
    Point start_1 = worldToPixel({x_min_overlap, y_max_overlap}, gt1);
    Point end_1   = worldToPixel({x_max_overlap, y_min_overlap}, gt1);
    // ensure start <= end
    int col1_s = std::min(start_1.c, end_1.c);
    int col1_e = std::max(start_1.c, end_1.c);
    int row1_s = std::min(start_1.r, end_1.r);
    int row1_e = std::max(start_1.r, end_1.r);

    Point start_2 = worldToPixel({x_min_overlap, y_max_overlap}, gt2);
    Point end_2   = worldToPixel({x_max_overlap, y_min_overlap}, gt2);
    int col2_s = std::min(start_2.c, end_2.c);
    int col2_e = std::max(start_2.c, end_2.c);
    int row2_s = std::min(start_2.r, end_2.r);
    int row2_e = std::max(start_2.r, end_2.r);

    // --- compute overlap width/height (inclusive pixel indices) ---
    if (col1_e < col1_s || row1_e < row1_s) {
        std::cerr << "Invalid overlap pixel range for image 1.\n";
        GDALClose(ds1); GDALClose(ds2);
        return false;
    }
    overlap_width  = static_cast<size_t>(col1_e - col1_s + 1);
    overlap_height = static_cast<size_t>(row1_e - row1_s + 1);

    if (overlap_width == 0 || overlap_height == 0) {
        std::cerr << "Invalid overlap region.\n";
        GDALClose(ds1); GDALClose(ds2);
        return false;
    }

    // top-left shift into image 1 and image 2
    shift1 = { row1_s, col1_s };
    shift2 = { row2_s, col2_s };

    // --- GDAL data type mapping ---
    GDALDataType gdal_dtype =
        std::is_same<T, unsigned char>::value ? GDT_Byte :
        std::is_same<T, uint16_t>::value      ? GDT_UInt16 :
        std::is_same<T, float>::value         ? GDT_Float32 :
        GDT_Unknown;

    if (gdal_dtype == GDT_Unknown) {
        std::cerr << "Unsupported pixel type.\n";
        GDALClose(ds1); GDALClose(ds2);
        return false;
    }

    // Allocate storage
    size_t full1 = nx1 * ny1;
    size_t full2 = nx2 * ny2;
    size_t ovsz = overlap_width * overlap_height;

    for (auto& ch : I1) ch.resize(full1);
    for (auto& ch : I2) ch.resize(full2);
    for (auto& ch : I1_overlap) ch.resize(ovsz);
    for (auto& ch : I2_overlap) ch.resize(ovsz);

    nodata_values1.assign(3, std::numeric_limits<double>::quiet_NaN());
    nodata_values2.assign(3, std::numeric_limits<double>::quiet_NaN());

    // --- read bands into full image buffers and overlap buffers ---
    for (int b = 0; b < 3; ++b) {
        if (!read_band<T>(ds1, b + 1, I1[b], I1_overlap[b], col1_s, row1_s,
                            nx1, ny1, overlap_width, overlap_height, nodata_values1[b], gdal_dtype) 
         || !read_band<T>(ds2, b + 1, I2[b], I2_overlap[b], col2_s, row2_s,
                            nx2, ny2, overlap_width, overlap_height, nodata_values2[b], gdal_dtype))
        {
            std::cerr << "Error reading band " << b + 1 << '\n';
            GDALClose(ds1); GDALClose(ds2);
            return false;
        }
    }

    GDALClose(ds1); GDALClose(ds2);

    // --- check that there is at least one pixel valid in both overlap images ---
    bool hasCommonValid = false;
    for (size_t i = 0; i < ovsz; ++i) {
        bool v1 = true, v2 = true;
        for (size_t ch = 0; ch < 3; ++ch) {
            if (I1_overlap[ch][i] == static_cast<T>(nodata_values1[ch])) v1 = false;
            if (I2_overlap[ch][i] == static_cast<T>(nodata_values2[ch])) v2 = false;
        }
        if (v1 && v2) { hasCommonValid = true; break; }
    }

    if (!hasCommonValid) {
        std::cerr << "Warning: no pixels are valid in both overlap images — disjoint coverage.\n";
        return false;
    }

    std::cout << "Overlap read: " << overlap_width << "x" << overlap_height << " pixels\n";


    // --- build full-image masks and extract perimeters ---
    uint8_t* mask1 = buildMask(I1, { nodata_values1[0], nodata_values1[1], nodata_values1[2] }, ny1, nx1);
    uint8_t* mask2 = buildMask(I2, { nodata_values2[0], nodata_values2[1], nodata_values2[2] }, ny2, nx2);

    size_t nPerim1 = 0, nPerim2 = 0;
    Point* perim1 = extractPerimeter(mask1, ny1, nx1, nPerim1);
    Point* perim2 = extractPerimeter(mask2, ny2, nx2, nPerim2);

    delete[] mask1; delete[] mask2;

    // --- convert perimeters to world coords ---
    PointWorld* perim1_w = toWorldCoords(perim1, nPerim1, gt1);
    PointWorld* perim2_w = toWorldCoords(perim2, nPerim2, gt2);

    delete[] perim1; delete[] perim2;

    // --- intersect perimeters in world space to get candidate points ---
    size_t nCandidates = 0;
    PointWorld* candidates_w = computePerimeterIntersections(perim1_w, nPerim1, perim2_w, nPerim2,
                                                             static_cast<float>(gt1[1]), nCandidates);

    delete[] perim1_w; delete[] perim2_w;

    std::vector<Point> candidates_px;
    size_t nValid = 0;

    for (size_t i = 0; i < nCandidates; i++) {
        Point px = worldToPixel(candidates_w[i], gt1); // px are in coordinates of image 1
        int r_overlap = px.r - shift1.r; // shift to coordinates of overlap
        int c_overlap = px.c - shift1.c;
        if (r_overlap >= 0 && r_overlap < overlap_height &&
            c_overlap >= 0 && c_overlap < overlap_width) 
        {
            bool notNull = true;
            for (size_t ch = 0; ch < 3; ch++){
                if(I1_overlap[ch][r_overlap * overlap_width + c_overlap] == nodata_values1[ch] ||
                   I2_overlap[ch][r_overlap * overlap_width + c_overlap] == nodata_values2[ch]) {
                    notNull = false;
                    break;
                }
            }
            if (notNull) { 
                nValid ++;
                candidates_px.push_back({ px.r, px.c }); 
            }
        }
    }

    delete[] candidates_w;

    // --- compute XOR diff mask covering the union of both images (pointer-returning function) ---
    diffMask = computeXORMask(I1, { nodata_values1[0], nodata_values1[1], nodata_values1[2] }, ny1, nx1, gt1,
                                    I2, { nodata_values2[0], nodata_values2[1], nodata_values2[2] }, ny2, nx2, gt2,
                                    unionH, unionW,
                                    shift1_union, shift2_union);

    // shift candidates into union-mask coordinates
    for (size_t i = 0; i < nValid; i++) {
        candidates_px[i].r += shift1_union.r;
        candidates_px[i].c += shift1_union.c;
    }

    // Debugging mask
    for (size_t i = 0; i < nValid; i++) {
        diffMask[candidates_px[i].r * unionW + candidates_px[i].c] = 5;
    }

    // --- select optimal endpoints on diff mask ---
    selectOptimalEndpoints(diffMask, unionH, unionW, candidates_px, nValid, Vb, Ve);

    auto out_of_bounds = [&](const Point& P){
        return (P.r < 0 || P.c < 0 ||
                P.r >= static_cast<int>(overlap_height) ||
                P.c >= static_cast<int>(overlap_width));
    };

    for (size_t i = 0; i < Vb.size(); i++) {
        // Debugging mask
        diffMask[Vb[i].r*unionW + Vb[i].c] = 6;
        diffMask[Ve[i].r*unionW + Ve[i].c] = 6;

        // convert endpoints back to overlap-image pixel coords
        Vb[i].r -= shift1_union.r + shift1.r;
        Vb[i].c -= shift1_union.c + shift1.c;
        Ve[i].r -= shift1_union.r + shift1.r;
        Ve[i].c -= shift1_union.c + shift1.c;

        if (out_of_bounds(Vb[i]) || out_of_bounds(Ve[i])) {
            std::cerr << "Warning: Vb or Ve out of overlap bounds\n"
                    << "Vb=(" << Vb[i].r << "," << Vb[i].c << ")  "
                    << "Ve=(" << Ve[i].r << "," << Ve[i].c << ")\n"
                    << "Overlap size=" << overlap_height << "x" << overlap_width << "\n";
        }
    }

    return true;
}

// --- TIFF write ---
template <typename T>
int io_tiff_modify(const char* fname, const ImageRGB& img, size_t nx, size_t ny) {
    GDALAllRegister();

    GDALDataset* ds = (GDALDataset*) GDALOpen(fname, GA_Update);
    if (!ds) {
        fprintf(stderr, "Failed to open GeoTIFF for update: %s\n", fname);
        return -1;
    }

    // Verify number of bands
    if (ds->GetRasterCount() < 3) {
        fprintf(stderr, "Expected at least 3 bands (RGB)\n");
        GDALClose(ds);
        return -1;
    }

    // Detect GDAL data type
    GDALDataType gdal_dtype =
        std::is_same<T, unsigned char>::value ? GDT_Byte :
        std::is_same<T, uint16_t>::value      ? GDT_UInt16 :
        std::is_same<T, float>::value         ? GDT_Float32 :
        GDT_Unknown;

    if (gdal_dtype == GDT_Unknown) {
        fprintf(stderr, "Unsupported data type\n");
        GDALClose(ds);
        return -1;
    }

    // Overwrite pixel data in each band
    for (size_t i = 0; i < 3; ++i) {
        GDALRasterBand* band = ds->GetRasterBand(static_cast<int>(i + 1));
        if (!band) {
            fprintf(stderr, "Missing band %zu\n", i + 1);
            GDALClose(ds);
            return -1;
        }

        std::vector<T> buf(nx * ny);
        for (size_t j = 0; j < nx * ny; ++j) {
            float val = img[i][j];
            if constexpr (std::is_integral<T>::value) {
                val = std::max(0.f, std::min(static_cast<float>(std::numeric_limits<T>::max()), val));
            }
            buf[j] = static_cast<T>(val);
        }

        if (band->RasterIO(GF_Write, 0, 0,
                           static_cast<int>(nx), static_cast<int>(ny),
                           buf.data(),
                           static_cast<int>(nx), static_cast<int>(ny),
                           gdal_dtype, 0, 0) != CE_None) {
            fprintf(stderr, "Failed to write band %zu\n", i + 1);
            GDALClose(ds);
            return -1;
        }
    }

    GDALClose(ds);
    return 0;
}



#endif // _IO_TIFF_H