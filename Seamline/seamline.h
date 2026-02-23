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


 #pragma once
// seamline_algorithm.h
// Header for hierarchical seamline optimization (Algorithm 2)

#include "io_tiff.h"
#include <vector>
#include <array>
#include <string>
#include <utility>
#include <cstddef>
#include <cmath>


// ------------------------------- Utilities ----------------------------------
inline std::size_t idx(std::size_t r, std::size_t c, std::size_t width) { return r * width + c; }
inline int clampi(int v, int a, int b) { return v < a ? a : (v > b ? b : v); }

// Directions for 8-connectivity
extern const int dr8[8];
extern const int dc8[8];
extern const float edgeLen8[8];

// ----------------------------------------------------------------------------
// Image processing functions
// ----------------------------------------------------------------------------
ImageRGB readImageOverview(
    const ImageRGB &in, 
    std::size_t width, std::size_t height, 
    std::size_t factor, 
    std::size_t &outW, std::size_t &outH);

ImageRGB moravec(
    const ImageRGB &img, 
    std::size_t width, std::size_t height);

// ----------------------------------------------------------------------------
// Energy computation
// ----------------------------------------------------------------------------
struct EnergyParams { float alpha = 1, beta = 1, gamma = 0; };
Channel computeEnergy(
    const ImageRGB &I1, const ImageRGB &I2,
    const std::vector<double> &nodata_values1,
    const std::vector<double> &nodata_values2,
    std::size_t width, std::size_t height,
    const EnergyParams &p,
    const ImageRGB *Wr = nullptr);

// ----------------------------------------------------------------------------
// Shortest path / seamline functions
// ----------------------------------------------------------------------------
std::vector<Point> shortestPath(
    const Channel &W, std::size_t width, std::size_t height, 
    Point Vb, Point Ve);

std::vector<std::vector<Point>> splitPathToBlocks(const std::vector<Point> &S, int blockLen = 20);

struct BBox { int r0, c0, r1, c1; };
BBox computeBoundingBox(
    const std::vector<Point> &seg, 
    int pad, std::size_t H, std::size_t W);

std::pair<Point, Point> computeEndPoints(const std::vector<Point> &seg);

std::pair<ImageRGB, ImageRGB> readFineImage(
    const ImageRGB &I1, const ImageRGB &I2, 
    std::size_t width, std::size_t height, const BBox &box, 
    std::size_t &outW, std::size_t &outH);

std::vector<Point> constructPathFromSegments(const std::vector<std::vector<Point>> &Sfine);

// ----------------------------------------------------------------------------
// Hierarchical seam line optimization
// ----------------------------------------------------------------------------
std::vector<Point> hierarchicalSeamLineOptimization(
    const ImageRGB &I1_overlap, const ImageRGB &I2_overlap,
    const std::vector<double> &nodata_values1, const std::vector<double> &nodata_values2,
    Point Vb, Point Ve,
    std::size_t width_overlap, std::size_t height_overlap,
    const EnergyParams &p);

// ----------------------------------------------------------------------------
// Union mask flood filling
// ----------------------------------------------------------------------------
void floodFillMask (
    std::vector<uint8_t>& diffMask,
    const std::vector<std::vector<Point>>& seamlines,
    const size_t unionW, const size_t unionH,
    const Point shift1, const Point shift1_union,
    const size_t nb1, const size_t nb2);

// ----------------------------------------------------------------------------
// Seam cutting
// ----------------------------------------------------------------------------
void cutAlongSeam(
    ImageRGB& I,
    const std::vector<uint8_t>& diffMask,
    const size_t width, const size_t height,
    const size_t unionW,
    const Point shift_union,
    const std::vector<double>& nodata_values,
    const size_t imgIdx,
    const bool keepBoarder);

// ----------------------------------------------------------------------------
// Main seamline interface
// ----------------------------------------------------------------------------
int Seamline(
    const std::string &img1, const std::string &img2, 
    size_t nb1, size_t nb2,
    const bool visibleSeam);