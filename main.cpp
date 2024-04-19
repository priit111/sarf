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

#include <iostream>

#include <algorithm>
#include <complex>
#include <thread>

#include <fmt/format.h>

#include "palsar/palsar_file.hpp"

#include "sar/sar_metadata.hpp"
#include "util/checks.hpp"

#include "sar/iq_correct.hpp"
#include "sar/sar_image.hpp"
#include "util/img_out.hpp"
#include "util/plot.hpp"
#include "util/util.hpp"

#include "sar/chirp.hpp"
#include "sar/fft_helper.hpp"
#include "sar/fractional_doppler_centroid.hpp"
#include "sar/processing_velocity_estimation.h"
#include "sar/range_compression.hpp"
#include "sar/range_doppler_algorithm.hpp"
#include "sar/terrain_correction.hpp"

void PlotChirp(const std::vector<std::complex<float>>& chirp) {
    PlotArgs a = {};
    a.out_path = "/tmp/chirp.html";
    a.graph_name = "chirp";

    Scatter i;
    Scatter q;

    i.line_name = "I";
    q.line_name = "Q";
    int x = 0;
    for (auto e : chirp) {
        i.y.push_back(e.real());
        q.y.push_back(e.imag());
        i.x.push_back(x);
        q.x.push_back(x);
        x++;
    }
    a.data = {i, q};

    Plot(a);
}

int main(int argc, const char* argv[]) {
    auto t = TimeStart();
    fftwf_init_threads();
    fftwf_plan_with_nthreads(std::thread::hardware_concurrency());


    if (argc != 2) {
        return 1;
    }

    const char* path = argv[1];
    palsar::FileSystem fs;
    fs.InitFromPath(path, "HH");

    SARMetadata metadata = {};

    auto parse_start = TimeStart();
    SARImage data;
    fs.ReadAllData(data, metadata);
    TimeStop(parse_start, "File parse");

#if 0
    const auto& osv = metadata.orbit_state_vectors;
    {
        std::vector<double> x;
        std::vector<double> y;
        for(size_t i = 0; i < osv.size(); i++)
        {
            x.push_back(osv[i].time);
            y.push_back(osv[i].x_pos);
        }
        std::vector<double> xi;
        std::vector<double> yi;
        for(int i = 0; i < 100; i++)
        {
            auto osv = InterpolateOrbit(metadata.orbit_state_vectors, CalcAzimuthTime(metadata, i));
            xi.push_back(osv.time);
            yi.push_back(osv.x_pos);
        }
        PlotArgs a = {};
        a.out_path = "./osv_test.html";
        a.graph_name = "test";
        a.x_axis_title = "time";
        a.y_axis_title = "x pos";
        a.data = {{.line_name = "OSV", .x = x, .y = y}, {.line_name="interp", .x=xi, .y=yi}};
        Plot(a);
        return 1;
    }
#endif

    BasicVrEstimation(metadata);
    fmt::print("Basic Vr = {}\n", metadata.results.Vr);
    metadata.results.Vr = EstimateProcessingVelocity(metadata, data.XSize(), data.YSize());

    GenerateChirpData(metadata.chirp, data.XStride());

    auto corr_time = TimeStart();
    CorrectIQ(data, metadata);
    TimeStop(corr_time, "I/Q correction");

    auto dc_time = TimeStart();
    EstimateFractionDopplerCentroid(data, metadata.pulse_repetition_frequency);
    TimeStop(dc_time, "Fract Doppler Centroid");

    WriteOut("/tmp/raw.tiff", data);

    auto rc_time = TimeStart();
    RunRangeCompression(data, metadata.chirp.freq_domain_data, metadata.chirp.n_samples);
    TimeStop(rc_time, "Range compression");

    WriteOut("/tmp/rc.tiff", data);
    

    auto az_start = TimeStart();
    auto az = RangeDopplerAlgorithm(metadata, data);
    TimeStop(az_start, "Azimuth compression");

    data.Clear();

    WriteOut("/tmp/az.tiff", az);

    TimeStop(t, "Total");

    RangeDopplerTerrainCorrection(az, metadata);

    return 0;
}
