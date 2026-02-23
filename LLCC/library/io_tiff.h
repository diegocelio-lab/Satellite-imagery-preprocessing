/**
 * Copyright (c) 2026, Diego Celio, diego.celio@me.com
 *
 * This file is part of the LLCC project.
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
 
#ifndef _IO_TIFF_H
#define _IO_TIFF_H

#include <tiffio.h>
#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <type_traits>
#include "gdal_priv.h"


// --- TIFF read ---
template <typename T>
inline T* io_tiff_read(const char* fname, size_t* nx, size_t* ny, size_t* nc, std::vector<double>* nodata_values = nullptr) {
    GDALAllRegister();

    GDALDataset* ds = (GDALDataset*)GDALOpen(fname, GA_ReadOnly);
    if (!ds) {
        fprintf(stderr, "Error: cannot open TIFF file %s\n", fname);
        return nullptr;
    }

    *nx = ds->GetRasterXSize();
    *ny = ds->GetRasterYSize();
    *nc = ds->GetRasterCount();

    T* buffer = new T[(*nx) * (*ny) * (*nc)];

    GDALDataType gdal_dtype = std::is_same<T, unsigned char>::value ? GDT_Byte :
                              std::is_same<T, uint16_t>::value ? GDT_UInt16 :
                              std::is_same<T, float>::value    ? GDT_Float32 :
                              GDT_Unknown;

    if (gdal_dtype == GDT_Unknown) {
        fprintf(stderr, "Unsupported data type\n");
        GDALClose(ds);
        delete[] buffer;
        return nullptr;
    }

    if (!buffer) {
        GDALClose(ds);
        return nullptr;
    }

    // --- Read each band ---
    for (size_t i = 0; i < *nc; i++) {
        GDALRasterBand* band = ds->GetRasterBand(i + 1);
        if (!band) {
            delete[] buffer;
            GDALClose(ds);
            return nullptr;
        }

        int success = 0;
        double nodata = band->GetNoDataValue(&success);
        if (nodata_values) {
            if (success)
                nodata_values->push_back(nodata);
            else
                nodata_values->push_back(std::numeric_limits<double>::quiet_NaN());
        }

        CPLErr err = band->RasterIO(GF_Read, 0, 0, static_cast<int>(*nx), static_cast<int>(*ny),
                                    buffer + i * (*nx) * (*ny),
                                    static_cast<int>(*nx), static_cast<int>(*ny),
                                    gdal_dtype, 0, 0);
        if (err != CE_None) {
            fprintf(stderr, "Error reading band %zu\n", i + 1);
            delete[] buffer;
            GDALClose(ds);
            return nullptr;
        }
    }

    GDALClose(ds);
    return buffer;
}

// --- TIFF write ---
template <typename T>
inline int io_tiff_write(const char* fname, const T* data, size_t nx, size_t ny, size_t nc, const char* fileref) {
    GDALAllRegister();

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!driver) return -1;

    GDALDataset* ref = (GDALDataset*)GDALOpen(fileref, GA_ReadOnly);

    GDALDataType gdal_dtype = std::is_same<T, unsigned char>::value ? GDT_Byte :
                              std::is_same<T, uint16_t>::value ? GDT_UInt16 :
                              std::is_same<T, float>::value    ? GDT_Float32 :
                              GDT_Unknown;
    
    if (gdal_dtype == GDT_Unknown) {
        fprintf(stderr, "Unsupported data type\n");
        if (ref) GDALClose(ref);
        return -1;
    }
    

    char** papszOptions = nullptr;

    GDALDataset* ds = driver->Create(fname, static_cast<int>(nx), static_cast<int>(ny),
                                     static_cast<int>(nc), gdal_dtype, papszOptions);
    if (!ds) {
        if (ref) GDALClose(ref);
        return -1;
    }

    if (ref) {
        ds->SetProjection(ref->GetProjectionRef());
        double geoTransform[6];
        if (ref->GetGeoTransform(geoTransform) == CE_None)
            ds->SetGeoTransform(geoTransform);
        GDALClose(ref);
    }

    for (size_t i = 0; i < nc; i++) {
        GDALRasterBand* band = ds->GetRasterBand(static_cast<int>(i + 1));
        if (!band) {
            GDALClose(ds);
            return -1;
        }

        CPLErr err = band->RasterIO(GF_Write, 0, 0, static_cast<int>(nx), static_cast<int>(ny),
                                    (void*)(data + i * nx * ny),
                                    static_cast<int>(nx), static_cast<int>(ny),
                                    gdal_dtype, 0, 0);
        if (err != CE_None) {
            GDALClose(ds);
            return -1;
        }

        band->SetNoDataValue(0.0);
    }

    GDALClose(ds);
    return 0;
}

#endif // _IO_TIFF_H

