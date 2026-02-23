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

 // seamline_algorithm.cpp
// Complete implementation of Algorithm 2 (hierarchical seamline optimization).
// Assumptions:
// - Input images I1 and I2 are provided as std::vector<float> in row-major order
//   and have the same width and height.
// - No external image libraries; everything operates on float buffers.
// - OpenMP used for parallel refinement.

#include "seamline.h"
#include "io_tiff.h"

#include <iostream>
#include <queue>
#include <limits>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <filesystem>
#include <omp.h>

using std::vector;
using std::size_t;
namespace fs = std::filesystem;

// ----------------------------------------------------------------------------
// readImageOverview – downsample the input images for coarse level processing
// ----------------------------------------------------------------------------
ImageRGB readImageOverview(
    const ImageRGB &in, 
    size_t width, size_t height,
    size_t factor, 
    size_t &outW, size_t &outH) 
{
    // If factor <= 1, return original image
    if (factor <= 1) { outW = width; outH = height; return in; }

    // Compute output dimensions, rounding up
    outW = (width + factor - 1) / factor;
    outH = (height + factor - 1) / factor;

    // Allocate output image for 3 channels
    ImageRGB out;
    for (auto& ch : out) 
        ch.resize(outW * outH, 0.0f);

    // Loop over RGB channels
    for (size_t b = 0; b < 3; ++b) {
        // Loop over output pixels
        for (size_t r = 0; r < outH; ++r) {
            for (size_t c = 0; c < outW; ++c) {
                // Determine input block bounds corresponding to this output pixel
                size_t  r0 = r * factor, c0 = c * factor;
                size_t  r1 = std::min(r0 + factor, height);
                size_t  c1 = std::min(c0 + factor, width);

                // Compute average over the input block
                float sum = 0;
                int count = 0;
                for (size_t rr = r0; rr < r1; ++rr) {
                    for (size_t cc = c0; cc < c1; ++cc) {
                        sum += in[b][rr * width + cc];
                        ++count;
                    }
                }

                // Store average in output pixel
                out[b][r * outW + c] = sum / count;
            }
        }
    }

    return out;
}

// ----------------------------------------------------------------------------
// Moravec operator – computes an "informativeness" measure (corner strength)
// ----------------------------------------------------------------------------
ImageRGB moravec(
    const ImageRGB &img, 
    size_t width, size_t height)
{
    // Allocate output image for 3 channels
    ImageRGB out;
    for (auto& ch : out)
        ch.resize(width * height, 0.0f);

    // Define 4 neighbor shifts for the Moravec operator (left, right, up, down)
    const int shift[4][2] = { {0,-1}, {0,1}, {-1,0}, {1,0} };

    // Loop over channels
    for (size_t b = 0; b < 3; ++b) {
        // Loop over all pixels
        for (size_t r = 0; r < height; ++r) {
            for (size_t c = 0; c < width; ++c) {
                float best = 1e9; // Initialize minimum SSD for this pixel

                // Loop over shifts (directions)
                for (int s = 0; s < 4; ++s) {
                    float ssd = 0; // Sum of squared differences

                    // Loop over 3x3 window around pixel
                    for (int wy = -1; wy <= 1; ++wy) {
                        for (int wx = -1; wx <= 1; ++wx) {
                            int r1 = clampi(r + wy, 0, static_cast<int>(height) - 1);
                            int c1 = clampi(c + wx, 0, static_cast<int>(width) - 1);
                            int r2 = clampi(r + wy + shift[s][0], 0, static_cast<int>(height) - 1);
                            int c2 = clampi(c + wx + shift[s][1], 0, static_cast<int>(width) - 1);

                            float d = img[b][r1 * width + c1] - img[b][r2 * width + c2];
                            ssd += d * d;
                        }
                    }

                    // Keep the minimum SSD over all shifts
                    best = std::min(best, ssd);
                }

                // Store result for this pixel
                out[b][r * width + c] = best;
            }
        }
    }

    return out;
}

// ----------------------------------------------------------------------------
// computeEnergy – build total energy W(x,y) = α*Ws + β*Wi + γ*Wr
// Ws: squared difference between images
// Wi: Moravec informativeness (cornerness)
// Wr: optional forbidden zone raster (penalty)
// Output is a single-channel image with mean energy over RGB channels
// ----------------------------------------------------------------------------
Channel computeEnergy(
    const ImageRGB &I1, const ImageRGB &I2,
    const vector<double> &nodata_values1, 
    const vector<double> &nodata_values2,
    size_t width, size_t height,
    const EnergyParams &p,
    const ImageRGB *Wr) 
{
    size_t N = width * height;

    // Allocate output energy channel
    Channel W(N, std::numeric_limits<float>::infinity());

    // Compute Moravec informativeness for both images
    ImageRGB Wi1 = moravec(I1, width, height);
    ImageRGB Wi2 = moravec(I2, width, height);

    // Loop over pixels
    for (size_t i = 0; i < N; i++) {
        bool invalid = false;
        float sumEnergy = 0.0;
        
        // Sum over RGB channels
        for (size_t  b = 0; b < 3; ++b) {
            if((std::isnan(nodata_values1[b]) == false && I1[b][i] == static_cast<float>(nodata_values1[b])) ||
                 (std::isnan(nodata_values2[b]) == false && I2[b][i] == static_cast<float>(nodata_values2[b]))) {
                invalid = true;          // Make sure the seamline will not go through nodata pixels
                break;
            }
            if (std::isnan(I1[b][i]) || std::isnan(I2[b][i])) { invalid = true; break; }
            
            float diff = I1[b][i] - I2[b][i];
            float Ws = diff * diff;   // squared difference
            float Wi = Wi1[b][i] + Wi2[b][i];                           // informativeness
            float Wrval = (Wr ? (*Wr)[b][i] : 0.0f);                    // optional raster
            sumEnergy += p.alpha * Ws + p.beta * Wi + p.gamma * Wrval;
        }

        if (invalid) {
            W[i] = std::numeric_limits<float>::infinity();
        } else {
            // Mean over R,G,B
            W[i] = sumEnergy / 3.0f;
        }
    }

    return W;
}


// ----------------------------------------------------------------------------
// shortestPath – Dijkstra (8-connected grid)
// ----------------------------------------------------------------------------
vector<Point> shortestPath(
    const Channel &W, size_t width, size_t height, 
    Point Vb, Point Ve) 
{
    const size_t N = width * height;
    const float INF = std::numeric_limits<float>::infinity();

    vector<float> dist(N, INF);
    vector<int> par(N, -1);
    
    struct Node { float d; size_t idx; };
    struct Cmp { bool operator()(const Node&a, const Node&b)const{return a.d > b.d;} };
    
    std::priority_queue<Node, vector<Node>, Cmp> pq;
    
    size_t s = Vb.r * width + Vb.c;
    size_t t = Ve.r * width + Ve.c;
    
    if (!std::isfinite(W[s]) || !std::isfinite(W[t])) {
        std::cout << "Weight map at start or at end is infinite. \n";
        return {};
    }

    dist[s] = 0.0f;
    pq.push({0.0f, s});

    while(!pq.empty()){
        auto cur = pq.top(); pq.pop();
        float d = cur.d;
        size_t u = cur.idx;
        if(d != dist[u]) continue;
        if(u == t) break;
        size_t r = u / width, c = u % width;
        for(size_t k = 0; k < 8; k++){
            int rr = static_cast<int>(r) + dr8[k], cc = static_cast<int>(c) + dc8[k];
            if(rr < 0 || rr >= static_cast<int>(height) || cc < 0 || cc >= static_cast<int>(width)) continue;
            size_t v = idx(rr, cc, width);
            if (!std::isfinite(W[v])) continue;
            float w = 0.5f * (W[u] + W[v]) * edgeLen8[k];
            float nd = d + w;
            if(nd < dist[v]){
                dist[v] = nd; 
                par[v] = k;
                pq.push({nd, v});
            }
        }
    }

    if(!std::isfinite(dist[t])) {
        std::cout << "Path not found \n";
        return {};
    }

    vector<Point> path;
    size_t cur = t;
    while(true){
        path.push_back({static_cast<int>(cur / width), static_cast<int>(cur % width)});
        if(cur == s) break;
        int k = par[cur];
        if (k < 0) {
            return {};
        }
        int pr = (cur / width) - dr8[k];
        int pc = (cur % width) - dc8[k];
        cur = idx(pr, pc, width);
    }
    std::reverse(path.begin(),path.end());
    return path;
}

// ----------------------------------------------------------------------------
// splitPathToBlocks – divide coarse seam line into blocks
// ----------------------------------------------------------------------------
vector<vector<Point>> splitPathToBlocks(const vector<Point> &S, int blockLen)
{
    vector<vector<Point>> B;
    for(size_t i=0; i < S.size(); i += blockLen){
        size_t j = std::min(S.size(),i+blockLen);
        B.push_back(vector<Point>(S.begin() + i,S.begin() + j));
    }
    return B;
}

// ----------------------------------------------------------------------------
// computeBoundingBox – pad around segment
// ----------------------------------------------------------------------------
BBox computeBoundingBox(
    const vector<Point>& seg, 
    int pad, size_t H, size_t W)
{
    int rmin = INT_MAX, cmin = INT_MAX, rmax = INT_MIN, cmax = INT_MIN;

    for(auto &p : seg){
        rmin = std::min(rmin, p.r); cmin = std::min(cmin, p.c);
        rmax = std::max(rmax, p.r); cmax = std::max(cmax, p.c);
    }

    rmin = clampi(rmin-pad, 0, static_cast<int>(H)-1); cmin = clampi(cmin-pad, 0, static_cast<int>(W)-1);
    rmax = clampi(rmax+pad, 0, static_cast<int>(H)-1); cmax = clampi(cmax+pad, 0, static_cast<int>(W)-1);

    return {rmin,cmin,rmax,cmax};
}

// ----------------------------------------------------------------------------
// computeEndPoints – first and last point of a segment
// ----------------------------------------------------------------------------
std::pair<Point,Point> computeEndPoints(const vector<Point> &seg)
{
    return {seg.front(), seg.back()};
}

// ----------------------------------------------------------------------------
// readFineImage – extract a subregion (box) from the input RGB images
// ----------------------------------------------------------------------------
std::pair<ImageRGB, ImageRGB> readFineImage(
    const ImageRGB &I1, const ImageRGB &I2,
    size_t width, size_t height, const BBox &box,
    size_t &outW, size_t &outH)
{
    // Compute output size based on the bounding box
    outW = static_cast<size_t>(box.c1 - box.c0 + 1);
    outH = static_cast<size_t>(box.r1 - box.r0 + 1);

    // Allocate output RGB channels
    ImageRGB F1, F2;
    for (size_t b = 0; b < 3; ++b) {
        F1[b].resize(outW * outH, 0.0f);
        F2[b].resize(outW * outH, 0.0f);
    }

    // Copy pixels from input images to output subregion
    for (size_t  r = 0; r < outH; r++) {
        for (size_t  c = 0; c < outW; c++) {
            for (size_t  b = 0; b < 3; b++) {
                // Map subregion coordinates to full image coordinates
                F1[b][r * outW + c] = I1[b][(box.r0 + r) * width + (box.c0 + c)];
                F2[b][r * outW + c] = I2[b][(box.r0 + r) * width + (box.c0 + c)];
            }
        }
    }

    return {F1, F2};
}

// ----------------------------------------------------------------------------
// constructPathFromSegments – merge refined segments
// ----------------------------------------------------------------------------
vector<Point> constructPathFromSegments(const vector<vector<Point>> &Sfine)
{
    vector<Point>S;
    for(size_t i = 0; i < Sfine.size(); ++i){
        if(Sfine[i].empty()) continue;
        if(!S.empty() && (S.back().r == Sfine[i].front().r && S.back().c == Sfine[i].front().c))
            S.insert(S.end(), Sfine[i].begin() + 1, Sfine[i].end());
        else
            S.insert(S.end(), Sfine[i].begin(), Sfine[i].end());
    }
    return S;
}

// ----------------------------------------------------------------------------
// Algorithm 2: Parallel Hierarchical Seam Line Optimization
// ----------------------------------------------------------------------------
vector<Point> hierarchicalSeamLineOptimization(
    const ImageRGB &I1_overlap, const ImageRGB &I2_overlap,
    const vector<double> &nodata_values1, const vector<double> &nodata_values2,
    Point Vb, Point Ve,
    size_t width_overlap, size_t height_overlap,
    const EnergyParams &p)
{
    Channel W = computeEnergy(I1_overlap, I2_overlap, nodata_values1, nodata_values2, width_overlap, height_overlap, p);

    vector<Point> S = shortestPath(W, width_overlap, height_overlap, Vb, Ve);

    /*
    // Step 1: Overview images (coarse)
    int factor = std::max(1, static_cast<int>(std::min(width_overlap, height_overlap) / 100));
    size_t oW = 0, oH = 0;
    ImageRGB O1_overlap = readImageOverview(I1_overlap, width_overlap, height_overlap, factor, oW, oH);
    ImageRGB O2_overlap = readImageOverview(I2_overlap, width_overlap, height_overlap, factor, oW, oH);

    // Map endpoints to coarse grid (integer coordinates)
    Point Vb_coarse = { Vb.r / factor, Vb.c / factor };
    Point Ve_coarse = { Ve.r / factor, Ve.c / factor };

    // clamp to coarse dims
    Vb_coarse.r = clampi(Vb_coarse.r, 0, static_cast<int>(oH) - 1); Vb_coarse.c = clampi(Vb_coarse.c, 0, static_cast<int>(oW) - 1);
    Ve_coarse.r = clampi(Ve_coarse.r, 0, static_cast<int>(oH) - 1); Ve_coarse.c = clampi(Ve_coarse.c, 0, static_cast<int>(oW) - 1);

    // Compute coarse energy
    Channel Wcoarse = computeEnergy(O1_overlap, O2_overlap, nodata_values1, nodata_values2, oW, oH, p);

    // Compute coarse seam line
    vector<Point> Scoarse = shortestPath(Wcoarse, oW, oH, Vb_coarse, Ve_coarse);
    if (Scoarse.empty()) return {};

    vector<Point>Sglobal; Sglobal.reserve(Scoarse.size());
    for(auto &pnt : Scoarse) {
        int gr = clampi(pnt.r * factor + factor/2, 0, height_overlap - 1);
        int gc = clampi(pnt.c * factor + factor/2, 0, width_overlap - 1);
        Sglobal.push_back({gr, gc});
    }

    // split into blocks and refine hierarchically
    auto B = splitPathToBlocks(Sglobal);
    vector<Point>S = Sglobal;

    for(int step = 2; step <= 3; step++){
        if(step > 2) B = splitPathToBlocks(S);
        vector<vector<Point>> Sfine(B.size());

        #pragma omp parallel for schedule(dynamic)
        for(size_t i = 0; i < B.size(); i++){
            auto &Bi = B[i];
            if(Bi.empty()) continue;
            BBox box = computeBoundingBox(Bi, 16, height_overlap, width_overlap); //padding of 16 pixels
            size_t fW = 0, fH = 0;
            // read local patch (read-only), safe in parallel
            auto Fpair = readFineImage(I1_overlap, I2_overlap, width_overlap, height_overlap, box, fW, fH);
            ImageRGB F1 = std::move(Fpair.first);
            ImageRGB F2 = std::move(Fpair.second);

            Channel Wfine = computeEnergy(F1, F2, nodata_values1, nodata_values2, fW, fH, p);

            auto [Vb_seg, Ve_seg] = computeEndPoints(Bi);
            Point localVb = { Vb_seg.r - box.r0, Vb_seg.c - box.c0 };
            Point localVe = { Ve_seg.r - box.r0, Ve_seg.c - box.c0 };
            // clamp local endpoints
            localVb.r = clampi(localVb.r, 0, static_cast<int>(fH) - 1); localVb.c = clampi(localVb.c, 0, static_cast<int>(fW) - 1);
            localVe.r = clampi(localVe.r, 0, static_cast<int>(fH) - 1); localVe.c = clampi(localVe.c, 0, static_cast<int>(fW) - 1);

            vector<Point> Slocal = shortestPath(Wfine, fW, fH, localVb, localVe);
            for(auto &q : Slocal){ q.r += box.r0; q.c += box.c0; }
            Sfine[i] = move(Slocal);
        }

        S = constructPathFromSegments(Sfine);
        if(step == 3) break;
    }
    */

    return S;
}

void floodFillMask (
    vector<uint8_t>& diffMask,
    const vector<vector<Point>>& seamlines,
    const size_t unionW, const size_t unionH,
    const Point shift1, const Point shift1_union,
    const size_t nb1, const size_t nb2)
{
    std::queue<Point> q;

    // Debugging mask
    //saveMaskAsPNGG(
    //    "beforeSeamlined_"+std::string("debug_diffmask_") + std::to_string(nb1) + "_" + std::to_string(nb2) + ".png", diffMask, unionH, unionW);

    for (size_t i = 0; i < seamlines.size(); i++) {
        for (auto& p : seamlines[i]) {
            size_t p_c = p.c + shift1.c + shift1_union.c;
            size_t p_r = p.r + shift1.r + shift1_union.r;
            diffMask[p_r * unionW + p_c] = 255; // seam barrier
        }
    }

    // Debugging mask
    //saveMaskAsPNGG(
    //    "Seamlined_"+std::string("debug_diffmask_") + std::to_string(nb1) + "_" + std::to_string(nb2) + ".png", diffMask, unionH, unionW);

    for (int r = 0; r < unionH; r++) {
        for (int c = 0; c < unionW; c++) {
            if (diffMask[r * unionW + c] == 1 || diffMask[r * unionW + c] == 2){
                q.push({r, c});
            }
        }
    }
    
    // BFS expansion
    while (!q.empty()) {
        auto [r, c] = q.front(); q.pop();
        uint8_t val = diffMask[r * unionW + c];

        for (size_t k = 0; k < 4; ++k) {
            int rr = r + dr4[k];
            int cc = c + dc4[k];
            if (rr < 0 || rr >= static_cast<int>(unionH) ||
                cc < 0 || cc >= static_cast<int>(unionW))
                continue;

            uint8_t neighbor = diffMask[rr * unionW + cc];

            // Skip barriers
            if (neighbor == 255 || neighbor == 3)
                continue;

            // Fill 4’s and 5's only if not already 1/2
            if (neighbor == 4 || neighbor == 5) {
                diffMask[rr * unionW + cc] = val; // overwrite with 1 or 2
                q.push({rr, cc});
            }
        }
    }

    
    // Fill all isolated overlapping or candidates (i.e. pixel = 4 or 5)
    for (int r = 0; r < unionH; r++) {
        for (int c = 0; c < unionW; c++) {
            size_t count1 = 0;
            size_t count2 = 0;
            if (diffMask[r * unionW + c] == 4 || diffMask[r * unionW + c] == 5){
                for (size_t k = 0; k < 4; ++k) {
                    int rr = r + dr4[k];
                    int cc = c + dc4[k];
                    if (rr < 0 || rr >= static_cast<int>(unionH) ||
                        cc < 0 || cc >= static_cast<int>(unionW))
                        continue;
                    if (diffMask[rr * unionW + cc] == 1) {
                        count1 ++;
                    } else if(diffMask[rr * unionW + cc] == 2) {
                        count2 ++;
                    }
                }
                if(count1 == count2) {
                    diffMask[r * unionW + c] = 1; //arbitrary
                    continue;
                }
                diffMask[r * unionW + c] = (count1 > count2) ? 1 : 2;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// cutAlongSeam - "cut" an image using a seam line by setting values to nodata
// ----------------------------------------------------------------------------
void cutAlongSeam(
    ImageRGB& I,
    const vector<uint8_t>& diffMask,
    const size_t width, const size_t height,
    const size_t unionW,
    const Point shift_union,
    const vector<double>& nodata_values,
    const size_t imgIdx,
    const bool keepBoarder,
    const bool visibleSeam) 
{
    auto maskIdx = [&](size_t r, size_t c) { 
        return (r + shift_union.r) * unionW + (c + shift_union.c); 
    };

    // Convert nodata values to float and handle NaN
    std::array<float, 3> nodata;
    for (size_t ch = 0; ch < 3; ++ch) {
        if (std::isnan(nodata_values[ch])) nodata[ch] = 0.0f;
        else nodata[ch] = static_cast<float>(nodata_values[ch]);

        // Set to nodata all pixel that have value 4 in mask
        for(size_t r = 0; r < height; ++r) {
            for(size_t c = 0; c < width; ++c) {

                
                // to visualise the seamline
                if(keepBoarder && visibleSeam) {
                    if(diffMask[maskIdx(r, c)] == 255 && ch==0) I[ch][r * width + c] = 60000.0f;
                    if(diffMask[maskIdx(r, c)] == 255 && ch!=0) I[ch][r * width + c] = 1.0f;
                }
                
               
                if(diffMask[maskIdx(r, c)] != imgIdx) {
                    if(diffMask[maskIdx(r, c)] == 255 && keepBoarder) continue;
                    I[ch][r * width + c] = nodata[ch];
                }

            }
        }
    }
        
}


// ----------------------------------------------------------------------------
// Main seam line computation function
// ----------------------------------------------------------------------------
int Seamline(
    const std::string &img1, const std::string &img2, 
    size_t nb1, size_t nb2,
    const bool visibleSeam)
{
    size_t nx1, nx2;
    size_t ny1, ny2;
    size_t overlap_w, overlap_h;
    ImageRGB I1, I2;
    ImageRGB I1_overlap, I2_overlap; 
    vector<uint8_t> diffMask;
    size_t unionW, unionH;
    vector<Point> Vb, Ve;
    Point shift1; Point shift1_union;
    Point shift2; Point shift2_union;
    vector<double> nodata_values1, nodata_values2;

    EnergyParams p; p.alpha = 1; p.beta = 0.5; p.gamma = 0;

    
    try {
        if (!io_tiff_read_overlap<float>(
                img1, img2,
                nx1, ny1,
                nx2, ny2,
                I1, I2,
                I1_overlap, I2_overlap,
                diffMask, unionW, unionH,
                Vb, Ve,
                overlap_w, overlap_h,
                shift1, shift1_union,
                shift2, shift2_union,
                nodata_values1, nodata_values2)) {
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    Channel W = computeEnergy(I1_overlap, I2_overlap, nodata_values1, nodata_values2, overlap_w, overlap_h, p);

    vector<vector<Point>> seamlines;
    for (size_t i = 0; i < Vb.size(); i++) {
        seamlines.push_back(shortestPath(W, overlap_w, overlap_h, Vb[i], Ve[i]));
    }
    
    floodFillMask (diffMask, seamlines, unionW, unionH, shift1, shift1_union, nb1, nb2);

    cutAlongSeam(I1, diffMask, nx1, ny1, unionW, shift1_union, nodata_values1, 1, true, visibleSeam);
    cutAlongSeam(I2, diffMask, nx2, ny2, unionW, shift2_union, nodata_values2, 2, false, visibleSeam);

    //saveMaskAsPNGG(
    //    std::string("debug_diffmask_") + std::to_string(nb1) + "_" + std::to_string(nb2) + ".png",
    //    diffMask, unionH, unionW);

    io_tiff_modify<uint16_t>(img1.c_str(), I1, nx1, ny1);
    io_tiff_modify<uint16_t>(img2.c_str(), I2, nx2, ny2);

    return 0;
}
