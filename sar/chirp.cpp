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

#include "chirp.hpp"

#include "fft_helper.hpp"

namespace
{
void ApplyHammingWindow(std::vector<IQ<float>>& data) {
    const size_t N = data.size();
    for (size_t i = 0; i < N; i++) {
        size_t n = i;
        double term = (2 * M_PI * n) / N;
        double m = 0.54 - 0.46 * cos(term);
        data[i].i *= m;
        data[i].q *= m;
    }
}
}


void GenerateChirpData(ChirpInfo& chirp, size_t range_fft_size)
{
    double Kr = chirp.Kr;

    const int64_t n_samples = chirp.n_samples;
    const double dt = 1.0 / chirp.range_sampling_rate;

    auto& td = chirp.time_domain_data;
    td.resize(n_samples);


    for (int64_t i = 0; i < n_samples; i++) {
        const double t = (i - n_samples / 2) * dt;
        const double phase = M_PI * Kr * t * t;

        IQ<float> iq = {};
        iq.i = cos(phase);
        iq.q = sin(phase);
        td[i] = iq;
    }

    ApplyHammingWindow(td);

    float scaling = sqrt(n_samples);
    for(auto& el : td)
    {
        el.i /= scaling;
        el.q /= scaling;
    }


    auto& fd = chirp.freq_domain_data;
    fd = td;
    fd.resize(range_fft_size, {0, 0});
    RunRangeFFT(fd.data(), fd.size(), 1, true);

}
