/**
* This program is free software; you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the Free
* Software Foundation; either version 3 of the License, or (at your option)
* any later version.
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
* more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, see http://www.gnu.org/licenses/
 */

#pragma once



#include "sar/sar_image.hpp"

#include <gdal/gdal_priv.h>

inline void WriteOut(const char* path, SARImage& img)
{
    GDALRegister_GTiff();
    int h = img.YSize();
    int w = img.XSize();
    auto ds = GetGDALDriverManager()->GetDriverByName("gtiff")->Create(path, w, h, 1, GDT_Float32, nullptr);
    auto b = ds->GetRasterBand(1);


    std::vector<float> row(w);
    auto data = img.Get();
    const int x_stride = img.XStride();
    for(int y = 0; y < h; y++)
    {
        for(int x = 0; x < w; x++)
        {
            auto f = data[y * x_stride + x];
            row[x] = f.i * f.i + f.q * f.q;
        }
        b->RasterIO(GF_Write, 0, y, w, 1, row.data(), w, 1, GDT_Float32, 0, 0);

    }
    GDALClose(ds);
    printf("file = %s done!\n", path);
}
