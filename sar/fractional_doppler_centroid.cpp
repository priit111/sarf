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


#include "fractional_doppler_centroid.hpp"

#include <cmath>
#include <thread>
#include <vector>

#include <fmt/format.h>

namespace {
struct Params {
    const IQ<float>* data;
    int x_stride;
    int x_size;
    int y_start;
    int y_end;
    int y_max;
};

void CalcPhaseDifference(const IQ<float>* data, int x_size, int x_stride, int y_start, int y_end, IQ<double>* result) {
    IQ<double> sum = {0, 0};
    for (int y = y_start; y < y_end; y++) {
        IQ<float> tmp = {0, 0};
        for (int x = 0; x < x_size; x++) {
            IQ<float> cur = data[x + y * x_stride];
            if (cur.i == 0.0 && cur.q == 0.0) {
                continue;
            }
            IQ<float> prev = data[x + (y - 1) * x_stride];
            if (prev.i == 0.0 && prev.q == 0.0) {
                continue;
            }
            prev.q = -prev.q;
            IQ<float> res = complex_mul(prev, cur);
            tmp.i += res.i;
            tmp.q += res.q;
        }
        sum.i += tmp.i;
        sum.q += tmp.q;
    }
    *result = sum;
}

}  // namespace

double EstimateFractionDopplerCentroid(const SARImage& img, double prf) {
    const size_t n_threads = std::thread::hardware_concurrency();
    std::vector<IQ<double>> results(n_threads);
    std::vector<std::thread> threads;

    const int x_size = img.XSize();
    const int x_stride = img.XStride();
    const int y_max = img.YSize();

    int y_step = (y_max / n_threads) + 1;
    const IQ<float>* data = img.Get();
    for (size_t i = 0; i < n_threads; i++) {
        int y_start = i * y_step;
        if (y_start == 0) {
            y_start = 1;
        }
        int y_end = y_start + y_step;
        y_end = std::min(y_end, y_max - 1);
        threads.push_back(std::thread(&CalcPhaseDifference, data, x_size, x_stride, y_start, y_end, &results[i]));
    }

    for (auto& t : threads) {
        t.join();
    }
    IQ<double> tot = {};
    for (const auto& el : results) {
        tot.i += el.i;
        tot.q += el.q;
    }

    double angle = atan2(tot.q, tot.i);
    double fraction = angle / (2 * M_PI);
    double doppler_centroid = prf * fraction;

    fmt::print("Doppler centroid = {}\n", doppler_centroid);

    return doppler_centroid;
}