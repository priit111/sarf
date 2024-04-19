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

#include "range_doppler_algorithm.hpp"

#include <cmath>
#include <thread>

#include "util/util.hpp"

#include "fft_helper.hpp"

namespace
{
constexpr int RCMC_SIZE = 16;
const float WINDOW[16] = {-1.38777878e-17, 1.67577197e-02, 7.70724198e-02, 2.00770143e-01,
                                       3.94012424e-01,  6.30000000e-01, 8.49229857e-01, 9.82157437e-01,
                                       9.82157437e-01,  8.49229857e-01, 6.30000000e-01, 3.94012424e-01,
                                       2.00770143e-01,  7.70724198e-02, 1.67577197e-02, -1.38777878e-17};

float sinc(float x) {
    const float pi = M_PI;
    if (x == 0) return 1.0f;
    return sinf(pi * x) / (pi * x);
}



void RunRCMC(const SARMetadata& metadata, int range_size, int range_stride, int azimuth_fft_size, int start_y, int end_y, const IQ<float>* data_in, IQ<float>* data_out)
{
    double slant_first = metadata.slant_range_first_sample;
    double range_spacing = metadata.range_spacing;
    double Vr = metadata.results.Vr;
    double prf = metadata.pulse_repetition_frequency;
    double lambda = metadata.wavelength;
    double fft_bin_step = prf / azimuth_fft_size;



    for(int y = start_y; y < end_y; y++)
    {
        double fn = 0.0f;
        if (y < azimuth_fft_size / 2) {
            fn = y * fft_bin_step;
        } else {
            fn = (y - azimuth_fft_size) * fft_bin_step;
        }
        for(int x = 0; x < range_size; x++)
        {
            double R0 = slant_first + x * range_spacing;
            double dR = (lambda * lambda * fn * fn) / (8 * Vr * Vr);
            dR *= R0;
            float shift = dR / range_spacing;
            const int range_walk_pixels = (int)std::round(shift);
            const int dst_idx = x + y * range_stride;
            const int N = RCMC_SIZE;
            float i_sum = 0.0f;
            float q_sum = 0.0f;
            float norm_sum = 0.0f;

            const int calc_x = x + range_walk_pixels;
            const int first_idx = calc_x - N / 2 + y * range_stride;

            for (int i = 0; i < N; i++) {
                if (calc_x - N / 2 < 0 || calc_x + N / 2 >= range_size) {
                    continue;
                }

                const IQ<float> data = data_in[first_idx + i];

                const float sinc_offset = shift - range_walk_pixels + N / 2;

                const float mult = WINDOW[i] * sinc(sinc_offset - i);
                i_sum += data.i * mult;
                q_sum += data.q * mult;
                norm_sum += mult;
            }

            if (norm_sum == 0.0f) {
                data_out[dst_idx] = {};
            } else {
                data_out[dst_idx] = {i_sum / norm_sum, q_sum / norm_sum};
            }
        }
    }
}

void RunReferenceMultiply(const SARMetadata& metadata, int range_size, int range_stride, int azimuth_fft_size, int y_start, int y_end, IQ<float>* data)
{
    const double slant_first = metadata.slant_range_first_sample;
    const double range_spacing = metadata.range_spacing;
    const double Vr = metadata.results.Vr;
    const double prf = metadata.pulse_repetition_frequency;
    const double lambda = metadata.wavelength;
    const double fft_bin_step = prf / azimuth_fft_size;

    auto ref_start = TimeStart();
    for(int y = y_start; y < y_end; y++)
    {
        double fn = 0.0f;
        if (y < azimuth_fft_size / 2) {
            fn = y * fft_bin_step;
        } else {
            fn = (y - azimuth_fft_size) * fft_bin_step;
        }
        for(int x = 0; x < range_size; x++)
        {
            double R0 = slant_first + x * range_spacing;
            const float Ka = (2 * Vr * Vr) / (lambda * R0);
            float phase = (-M_PI * fn * fn) / Ka;

            float sin_val = sin(phase);
            float cos_val = cos(phase);

            auto& d = data[x + y * range_stride];

            d = complex_mul(d, {cos_val, sin_val});
        }
    }
}
}

SARImage RangeDopplerAlgorithm(const SARMetadata& metadata, SARImage& in)
{

    IQ<float>* data_in = in.Get();
    int range_sz = in.XSize();
    int range_stride = in.XStride();
    int azimuth_sz = in.YStride();
    int azimuth_fft_size = in.YStride();

#if 1
    auto fft_forward = TimeStart();
    RunAzimuthFFT((fftwf_complex*)data_in, range_sz, range_stride, azimuth_sz, true);
    TimeStop(fft_forward, "forward Azimuth FFT");
#else
    auto fft_forward = TimeStart();
    std::unique_ptr<IQ<float>[]> tmp_buf(new IQ<float>[range_stride * azimuth_fft_size]);
    Transpose(data_in, tmp_buf.get(), range_stride, azimuth_fft_size);
    TimeStop(fft_forward, "First transpose");
    fft_forward = TimeStart();
    RunRangeFFT(tmp_buf.get(), azimuth_fft_size, range_sz, true);
    TimeStop(fft_forward, "transposed Azimuth FFT");
    fft_forward = TimeStart();
    Transpose(tmp_buf.get(), data_in, azimuth_fft_size, range_stride);
    TimeStop(fft_forward, "Second transpose");

#endif


    SARImage out;
    out.Init(range_sz, in.YSize(), range_stride, azimuth_sz);

    IQ<float>* data_out = out.Get();

    // slant range at range pixel

    auto rcmc_start = TimeStart();

    auto n_threads =  std::thread::hardware_concurrency();

    fmt::print("N THREADS = {}\n", n_threads);
    std::vector<std::thread> thread_vec;
    uint32_t y_step = (azimuth_fft_size / n_threads) + 1;

    int start_y = 0;
    for(int i = 0; i < n_threads; i++)
    {
        start_y = i * y_step;
        int end_y = start_y + y_step;
        end_y = std::min(end_y, azimuth_fft_size);
        std::thread t(RunRCMC, metadata, range_sz, range_stride, azimuth_fft_size, start_y, end_y, data_in, data_out);
        thread_vec.push_back(std::move(t));
    }


    for(auto& t : thread_vec)
    {
        t.join();
    }

    TimeStop(rcmc_start, "RCMC");

    auto ref_start = TimeStart();
    thread_vec.clear();
    start_y = 0;
    for(int i = 0; i < n_threads; i++)
    {
        start_y = i * y_step;
        int end_y = start_y + y_step;
        end_y = std::min(end_y, azimuth_fft_size);
        std::thread t(RunReferenceMultiply, metadata, range_sz, range_stride, azimuth_fft_size, start_y, end_y,data_out);
        thread_vec.push_back(std::move(t));
    }

    for(auto& t : thread_vec){
        t.join();
    }


    TimeStop(ref_start, "reference multiply");

    auto fft_back = TimeStart();
    RunAzimuthFFT((fftwf_complex*)data_out, range_sz, range_stride, azimuth_sz, false);
    TimeStop(fft_back, "reverse azimuth FFT");

    float mult = 1.0 / azimuth_sz;

    int cutoff = CalcAperturePixels(metadata, range_sz - 1);

    for(int y = 0; y < azimuth_sz; y++)
    {
        for(int x = 0; x < range_sz; x++)
        {
            auto& d = data_out[x + range_stride * y];
            if(y < cutoff / 2 || (y > azimuth_sz - cutoff / 2))
            {
                d = {};
            }
            else
            {
                d.i *= mult;
                d.q *= mult;
            }
        }
    }

    return out;
}