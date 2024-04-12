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

#include <vector>

#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include "iq32f.hpp"

#include "util/geo_tools.hpp"
#include "osv.hpp"

struct ChirpInfo {
    double range_sampling_rate;
    double pulse_duration;
    int n_samples;
    double pulse_bandwidth;
    double Kr;
    std::vector<IQ<float>> time_domain_data;
    std::vector<IQ<float>> freq_domain_data;
};

struct ImgFormat {
    int range_size;
    int azimuth_size;
    int data_line_offset;
    int record_length;
};

struct SARResults {
    size_t total_samples;
    double dc_i;
    double dc_q;
    double Vr;
    double doppler_centroid;
};


inline double PTimeToDouble(boost::posix_time::ptime ptime)
{
    std::string ts("2000-01-01 00:00:00");
    boost::posix_time::ptime t(boost::gregorian::from_simple_string(ts));
    return (ptime - t).total_microseconds() * 1e-6;
}

struct SARMetadata {
    ChirpInfo chirp;
    std::string polarisation;
    std::vector<uint32_t> left_range_offsets;
    std::vector<uint32_t> echo_samples;
    double orbit_interval;
    std::vector<OrbitStateVector> orbit_state_vectors;
    boost::posix_time::ptime first_line_time;
    double pulse_repetition_frequency;
    double line_time_interval;
    double azimuth_bandwidth_fraction;
    double carrier_frequency;
    double wavelength;
    double platform_velocity;
    double range_spacing;
    double azimuth_spacing;
    double slant_range_first_sample;
    GeoPos3D center_point;
    boost::posix_time::ptime center_time;

    SARResults results;
};

inline double CalcSlantRange(const SARMetadata& metadata, int range_pixel)
{
    return metadata.slant_range_first_sample + metadata.range_spacing * range_pixel;
}

inline double CalcKa(const SARMetadata& metadata, int range_pixel)
{
    const double R0 = CalcSlantRange(metadata, range_pixel);
    const double Vr = metadata.results.Vr;
    constexpr double SOL = 299792458;
    return (2 * metadata.carrier_frequency * Vr * Vr) / (SOL * R0);
}

inline double CalcAperturePixels(const SARMetadata& metadata, int range_pixel)
{
    const double Ka = CalcKa(metadata, range_pixel);
    const double prf = metadata.pulse_repetition_frequency;
    return prf * prf / Ka;
}

inline double CalcAzimuthTime(const SARMetadata& metadata, int azimuth_idx)
{
    return metadata.line_time_interval * azimuth_idx;
}

inline void BasicVrEstimation(SARMetadata& metadata)
{
    auto osv = metadata.orbit_state_vectors.front();
    double Vs = sqrt(osv.x_vel * osv.x_vel + osv.y_vel * osv.y_vel + osv.z_vel * osv.z_vel);
    metadata.results.Vr = Vs * 0.94;
}
