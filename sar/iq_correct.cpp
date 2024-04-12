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

#include "iq_correct.hpp"

#include <fmt/format.h>

namespace {
void DCBiasCorrection(SARImage& img, SARMetadata& metadata) {
    const int az_size = img.YSize();
    const int range_stride = img.XStride();
    size_t total_samples = 0;
    IQ<float>* data = img.Get();
    double i_total = 0.0;
    double q_total = 0.0;
    for (int y = 0; y < az_size; y++) {
        const int range_offset = metadata.left_range_offsets[y];
        const int n_samples = metadata.echo_samples[y];
        double i_sum = 0.0;
        double q_sum = 0.0;
        for (int x = 0; x < n_samples; x++) {
            const auto& d = data[range_offset + x + y * range_stride];
            i_sum += d.i;
            q_sum += d.q;
        }

        i_total += i_sum;
        q_total += q_sum;
        total_samples += n_samples;
    }

    float i_bias = i_total / total_samples;
    float q_bias = q_total / total_samples;
    fmt::print("I/Q DC bias = {} {}\n", i_bias, q_bias);

    for (int y = 0; y < az_size; y++) {
        const int range_offset = metadata.left_range_offsets[y];
        const int n_samples = metadata.echo_samples[y];
        for (int x = 0; x < n_samples; x++) {
            auto& d = data[range_offset + x + y * range_stride];
            d.i -= i_bias;
            d.q -= q_bias;
        }
    }
}

void GainCorrection(SARImage& img, SARMetadata& metadata) {
    const int az_size = img.YSize();
    const int range_stride = img.XStride();
    size_t total_samples = 0;
    IQ<float>* data = img.Get();
    double i_total = 0.0;
    double q_total = 0.0;
    for (int y = 0; y < az_size; y++) {
        const int range_offset = metadata.left_range_offsets[y];
        const int n_samples = metadata.echo_samples[y];
        double i_sum = 0.0;
        double q_sum = 0.0;
        for (int x = 0; x < n_samples; x++) {
            const auto& d = data[range_offset + x + y * range_stride];
            i_sum += d.i * d.i;
            q_sum += d.q * d.q;
        }

        i_total += i_sum;
        q_total += q_sum;
    }

    fmt::print("I/Q sq tot = {} {}\n", i_total, q_total);
    float correction = sqrt(i_total / q_total);
    fmt::print("gain bias = {}\n", correction);

    for (int y = 0; y < az_size; y++) {
        const int range_offset = metadata.left_range_offsets[y];
        const int n_samples = metadata.echo_samples[y];
        for (int x = 0; x < n_samples; x++) {
            auto& d = data[range_offset + x + y * range_stride];
            d.q *= correction;
        }
    }
}

void PhaseCorrection(SARImage& img, SARMetadata& metadata) {
    const int az_size = img.YSize();
    const int range_stride = img.XStride();
    size_t total_samples = 0;
    IQ<float>* data = img.Get();
    double iq_total = 0.0;
    double i_sq_total = 0.0;
    for (int y = 0; y < az_size; y++) {
        const int range_offset = metadata.left_range_offsets[y];
        const int n_samples = metadata.echo_samples[y];
        double iq_sum = 0.0;
        double i_sq_sum = 0.0;
        for (int x = 0; x < n_samples; x++) {
            const auto& d = data[range_offset + x + y * range_stride];
            iq_sum += d.i * d.q;
            i_sq_sum += d.i * d.i;
        }

        iq_total += iq_sum;
        i_sq_total += i_sq_sum;
    }

    double delta_theta = iq_total / i_sq_total;
    float sin_corr = delta_theta;
    float cos_corr = (1 - delta_theta * delta_theta);
    double quad_misconf = asin(delta_theta);

    double deg_misconf = 360 * quad_misconf / (2 * M_PI);

    fmt::print("I/Q phase error = {}\n", deg_misconf);

    for (int y = 0; y < az_size; y++) {
        const int range_offset = metadata.left_range_offsets[y];
        const int n_samples = metadata.echo_samples[y];
        for (int x = 0; x < n_samples; x++) {
            auto& d = data[range_offset + x + y * range_stride];
            float i = d.i;
            float q = d.q;

            float corr_q = (q - i * sin_corr) / cos_corr;
            d.q = corr_q;
        }
    }
}
}  // namespace

void CorrectIQ(SARImage& img, SARMetadata& metadata) {
    DCBiasCorrection(img, metadata);
    GainCorrection(img, metadata);
    PhaseCorrection(img, metadata);
}
