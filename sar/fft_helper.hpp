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

#include <map>

#include <fftw3.h>
#include <fmt/format.h>

#include "iq32f.hpp"

inline void RunRangeFFT(IQ<float>* data, int range_size, int azimuth_size, bool forward) {
    int rank = 1;
    int n[] = {range_size};
    int how_many = azimuth_size;
    int inembed[1] = {range_size};
    int istride = 1;
    int idist = range_size;
    int onembed[1] = {range_size};
    int ostride = 1;
    int odist = range_size;
    int sign = forward ? FFTW_FORWARD : FFTW_BACKWARD;

    static_assert(sizeof(fftwf_complex) == sizeof(IQ<float>));
    fftwf_complex* in = reinterpret_cast<fftwf_complex*>(data);
    fftwf_complex* out = reinterpret_cast<fftwf_complex*>(data);

    fftwf_plan plan = fftwf_plan_many_dft(rank, n, how_many, in, inembed, istride, idist, out, onembed, ostride, odist,
                                          sign, FFTW_ESTIMATE);

    fftwf_execute_dft(plan, in, out);

    fftwf_destroy_plan(plan);
};

template <class T>
void Transpose(const T* in, T* out, int x_size, int y_size) {
    for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
            int in_idx = x + y * x_size;
            int out_idx = y + x * y_size;
            out[out_idx] = in[in_idx];
        }
    }
}

#if 0
inline void Transpose(IQ<float>* data, int range_size, int azimuth_size)
{
    //http://agentzlerich.blogspot.com/2010/01/using-fftw-for-in-place-matrix.html

    int flags = FFTW_ESTIMATE;
    fftw_iodim howmany_dims[2];
    howmany_dims[0].n  = azimuth_size;
    howmany_dims[0].is = range_size;
    howmany_dims[0].os = 1;
    howmany_dims[1].n  = range_size;
    howmany_dims[1].is = 1;
    howmany_dims[1].os = azimuth_size;
    int howmany_rank = 2;
    fftwf_complex* in = reinterpret_cast<fftwf_complex*>(data);
    fftwf_complex* out = reinterpret_cast<fftwf_complex*>(data);
    auto plan = fftwf_plan_guru_dft(/*rank*/0, /*dims*/NULL,
                                    howmany_rank, howmany_dims,
                                    in, out, /*kind*/NULL, flags);

    fftwf_execute_dft(plan, in, out);
    fftwf_destroy_plan(plan);

    double* in = reinterpret_cast<double*>(data);
    double* out = reinterpret_cast<double*>(data);
    int flags = FFTW_ESTIMATE;
    fftw_iodim howmany_dims[2];
    howmany_dims[0].n  = azimuth_size;
    howmany_dims[0].is = range_size;
    howmany_dims[0].os = 1;
    howmany_dims[1].n  = range_size;
    howmany_dims[1].is = 1;
    howmany_dims[1].os = azimuth_size;
    int howmany_rank = 2;
    auto plan = fftw_plan_guru_r2r(/*rank*/0, /*dims*/NULL,
                                    howmany_rank, howmany_dims,
                                    in, out, /*kind*/NULL, flags);

    fftw_execute_r2r(plan, in, out);
    fftw_destroy_plan(plan);

};
#endif

inline void RunAzimuthFFT(fftwf_complex* data, int range_size, int range_stride, int azimuth_size, bool forward) {
    int rank = 1;
    int n[] = {azimuth_size};
    int how_many = range_size;
    int inembed[1] = {azimuth_size};
    int istride = range_stride;
    int idist = 1;
    int onembed[1] = {azimuth_size};
    int ostride = range_stride;
    int odist = 1;
    int sign = forward ? FFTW_FORWARD : FFTW_BACKWARD;

    fftwf_complex* in = data;
    fftwf_complex* out = data;

    fftwf_plan plan = fftwf_plan_many_dft(rank, n, how_many, in, inembed, istride, idist, out, onembed, ostride, odist,
                                          sign, FFTW_ESTIMATE);

    fftwf_execute_dft(plan, in, out);

    fftwf_destroy_plan(plan);
};

inline int IntPow(int base, int power) {
    int r = 1;
    for (int i = 0; i < power; i++) {
        r *= base;
    }
    return r;
}

inline std::map<int, std::array<int, 4>> GetOptimalFFTSizes() {
    std::map<int, std::array<int, 4>> result;

    // cuFFT documentation states that the most optimal FFT size is expressed with the following: 2^a×3^b×5^c×7^d
    // find all values up to a reasonable limit that makes sense for PALSAR
    const int limit = 131072;
    for (int i = 0;; i++) {
        int p1 = IntPow(2, i);
        if (p1 > limit) {
            break;
        }
        for (int j = 0;; j++) {
            int p2 = IntPow(3, j);
            if (p1 * p2 > limit) {
                break;
            }
            for (int k = 0;; k++) {
                int p3 = IntPow(5, k);
                if (p1 * p2 * p3 > limit) {
                    break;
                }
                for (int l = 0;; l++) {
                    int p4 = IntPow(7, l);
                    int val = p1 * p2 * p3 * p4;
                    if (val <= limit) {
                        result[val] = {i, j, k, l};
                    } else {
                        break;
                    }
                }
            }
        }
    }
    return result;
}
