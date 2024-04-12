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

#include "range_compression.hpp"

#include "util/util.hpp"
#include "util/checks.hpp"
#include "fft_helper.hpp"

void RunRangeCompression(SARImage& img, const std::vector<IQ<float>>& fd_chirp, int chirp_size)
{
    SARF_CHECK(img.XStride() == fd_chirp.size());
    const int range_fft_size = img.XStride();
    const int range_size = img.XSize();
    const int az_size = img.YSize();

    IQ<float>* data = img.Get();
    auto fft_forward = TimeStart();
    RunRangeFFT(data, range_fft_size, az_size, true);
    TimeStop(fft_forward, "forward range FFT");

    std::vector<IQ<float>> conj_fd_chirp = fd_chirp;
    for(auto& el : conj_fd_chirp)
    {
        el.q = -el.q;
    }

    for (int y = 0; y < az_size; y++) {
        for (int x = 0; x < range_fft_size; x++) {
            auto& d = data[x + range_fft_size * y];

            d = complex_mul(d, conj_fd_chirp[x]);
        }
    }


    RunRangeFFT(data, range_fft_size, az_size, false);

    float scaling = 1.0 / range_fft_size;
    int new_range_size = range_size - chirp_size;

    for (int y = 0; y < az_size; y++) {
        for (int x = 0; x < range_fft_size; x++) {
            auto& d = data[x + range_fft_size * y];
            if(x < new_range_size)
            {
                d.i *= scaling;
                d.q *= scaling;
            }
            else
            {
                d = {0, 0};
            }
        }
    }

    img.ResizeX(new_range_size);

}