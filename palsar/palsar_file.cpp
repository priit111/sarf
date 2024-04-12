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

#include "palsar_file.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/date_time/posix_time/ptime.hpp>

#include <fmt/format.h>

#include "sar/fft_helper.hpp"

namespace {
constexpr int SIGNAL_OFFSET = 412;
}

namespace palsar {
void MetadataFile::Open(const std::string& path) {
    FILE* fp = fopen(path.c_str(), "rb");

    size_ = std::filesystem::file_size(path);

    data_ = std::unique_ptr<char[]>(new char[size_]);

    if (fp) {
        size_t n_read = fread(data_.get(), 1, size_, fp);
        path_ = path;
        fclose(fp);
    } else {
        throw std::runtime_error("Failed to open file " + path);
    }
}

MetadataFile::~MetadataFile() {}

void ImgFile::Open(const std::string& path) {
    file_size_ = std::filesystem::file_size(path);
    FILE* fp = fopen(path.c_str(), "rb");
    data_ = std::unique_ptr<char[]>(new char[file_size_]);
    if (fp) {
        size_t read_n = fread(data_.get(), 1, file_size_, fp);
        path_ = path;
        fclose(fp);
    } else {
        throw std::runtime_error("Failed to open file " + path);
    }
}

void FileSystem::InitFromPath(const char* path, const char* polarization) {
    std::filesystem::path p(path);
    for (const auto& entry : std::filesystem::directory_iterator{p}) {
        auto path_str = entry.path().string();
        auto filename = entry.path().filename().string();
        std::vector<std::string> tokens;
        boost::algorithm::split(tokens, filename, boost::is_any_of("-"));
        if (tokens.size() < 3U) {
            continue;
        }

        const auto& file_type = tokens.at(0);

        if (!volume_directory_file_.IsOpen() && file_type == "VOL") {
            volume_directory_file_.Open(path_str);
        } else if (!trailer_file_.IsOpen() && file_type == "TRL") {
            trailer_file_.Open(path_str);
        } else if (!leader_file_.IsOpen() && file_type == "LED") {
            leader_file_.Open(path_str);
            scene_id_ = tokens.at(1);
            product_id_ = tokens.at(2);
            fmt::print("Scene ID = {}\n", scene_id_);
        } else if (!image_file_.IsOpen() && file_type == "IMG" && tokens.size() == 4U) {
            const auto& pol = tokens.at(1);
            if (pol == polarization) {
                image_file_.Open(path_str);
            }
        }
    }

    std::string err_msg;
    if (!volume_directory_file_.IsOpen()) {
        err_msg += " VOL ";
    }
    if (!trailer_file_.IsOpen()) {
        err_msg += " TRL ";
    }
    if (!leader_file_.IsOpen()) {
        err_msg += " LED ";
    }
    if (!image_file_.IsOpen()) {
        err_msg += " IMG ";
    }

    if (!err_msg.empty()) {
        throw std::invalid_argument("Missing the following files: " + err_msg);
    }
}

void FileSystem::ReadAllData(SARImage& img, SARMetadata& metadata) {
    Record data_set_summary = GetDataSetSummaryRecord();

    std::string processing_level = data_set_summary.ReadAn(1095, 16);
    processing_level.erase(std::remove(processing_level.begin(), processing_level.end(), ' '), processing_level.end());
    if (processing_level != "1.0") {
        std::string err_msg = "Expected file processing level to be 1.0, actually is = " + processing_level;
        throw std::runtime_error(err_msg);
    }

    metadata.wavelength = data_set_summary.ReadF16(501);
    metadata.carrier_frequency = data_set_summary.ReadF16(493) * 1e9;

    auto& chirp = metadata.chirp;

    chirp.Kr = data_set_summary.ReadF16(551);
    /*
    chirp.coefficient[0] = data_set_summary.ReadF16(535);  // Hz
    chirp.coefficient[1] = data_set_summary.ReadF16(551);  // Hz/s
    chirp.coefficient[2] = data_set_summary.ReadF16(567);  // quad term
    chirp.coefficient[3] = data_set_summary.ReadF16(583);  // quad term
    chirp.coefficient[4] = data_set_summary.ReadF16(599);  // quad term
     */

    chirp.range_sampling_rate = data_set_summary.ReadF16(711) * 1e6;  // MHz -> Hz;
    chirp.pulse_duration = data_set_summary.ReadF16(743) * 1e-6;      // useconds -> s
    chirp.n_samples = static_cast<int>(std::round(chirp.pulse_duration * chirp.range_sampling_rate));

    chirp.pulse_bandwidth = chirp.pulse_duration * fabs(chirp.Kr);

    constexpr double SOL = 299792458;
    metadata.range_spacing = SOL / (2 * chirp.range_sampling_rate);

    metadata.pulse_repetition_frequency = data_set_summary.ReadF16(935) * 1e-3;  // mHz -> Hz;
    metadata.azimuth_bandwidth_fraction = 0.8;

    if (data_set_summary.ReadI4(389) == 4) {
        // in polarimetry mode PRF needs to be divided by 2
        metadata.pulse_repetition_frequency /= 2.0;
    }

    metadata.line_time_interval = 1 / metadata.pulse_repetition_frequency;

    double center_lat = data_set_summary.ReadF16(117);
    double center_lon = data_set_summary.ReadF16(133);

    fmt::print("Center point lat lon = {} {}\n", center_lat, center_lon);

    metadata.center_point = Geo2xyzWgs84(center_lat, center_lon, 0);

    auto center_time_str = data_set_summary.ReadAn(69, 32);

    // format YYYYMMDDhhmmssttt
    int year = std::stoi(center_time_str.substr(0, 4));
    int month = std::stoi(center_time_str.substr(4, 2));
    int day = std::stoi(center_time_str.substr(6, 2));
    int hour = std::stoi(center_time_str.substr(8, 2));
    int minute = std::stoi(center_time_str.substr(10, 2));
    int second = std::stoi(center_time_str.substr(12, 2));
    int millisecond = std::stoi(center_time_str.substr(14, center_time_str.size() - 14));

    boost::posix_time::ptime center_time(boost::gregorian::date(year, month, day));

    center_time += boost::posix_time::hours(hour);
    center_time += boost::posix_time::minutes(minute);
    center_time += boost::posix_time::seconds(second);
    center_time += boost::posix_time::milliseconds(millisecond);

    metadata.center_time = center_time;

    size_t azimuth_size = 0;
    size_t echo_samples = 0;
    size_t record_length = 0;

    {
        Record file_descriptor_record = GetFileDescriptorRecord();
        azimuth_size = file_descriptor_record.ReadI6(181);
        record_length = file_descriptor_record.ReadI6(187);
        echo_samples = file_descriptor_record.ReadI8(281) / 2;

        if (record_length != echo_samples * 2 + SIGNAL_OFFSET) {
            // TODO(priit)
            throw std::runtime_error("IMG files with differing signal lengths not yet supported");
        }
    }

    {
        auto r = GetSignalDataRecord(0, record_length);
        auto year = r.ReadB4(37);
        auto day_of_year = r.ReadB4(41);
        auto ms_of_day = r.ReadB4(45);
        fmt::print("FIRST LINE Y = {} {} {}\n", year, day_of_year, ms_of_day);

        auto d = boost::gregorian::date(year, 1, 1);
        d += boost::gregorian::days(day_of_year - 1);
        metadata.first_line_time = boost::posix_time::ptime(d, boost::posix_time::milliseconds(ms_of_day));
    }

    {
        OrbitStateVector first_pos{};
        Record platform_position_data = GetPlatformPositionDataRecord();
        first_pos.x_pos = platform_position_data.ReadF16(45);
        first_pos.y_pos = platform_position_data.ReadF16(45 + 1 * 16);
        first_pos.z_pos = platform_position_data.ReadF16(45 + 2 * 16);
        first_pos.x_vel = platform_position_data.ReadF16(45 + 3 * 16);
        first_pos.y_vel = platform_position_data.ReadF16(45 + 4 * 16);
        first_pos.z_vel = platform_position_data.ReadF16(45 + 5 * 16);

        const int n_points = platform_position_data.ReadI4(141);

        const int year = platform_position_data.ReadI4(145);
        const int month = platform_position_data.ReadI4(149);
        const int day = platform_position_data.ReadI4(153);
        const double second_of_day = platform_position_data.ReadD22(161);
        metadata.orbit_interval = platform_position_data.ReadD22(183);

        auto first_osv = boost::posix_time::ptime(boost::gregorian::date(year, month, day),
                                                  boost::posix_time::seconds(static_cast<int>(second_of_day)));

        double first_osv_time = PTimeToDouble(first_osv) - PTimeToDouble(metadata.first_line_time);

        std::cout << "first osv = " << first_osv << "\n";
        std::cout << "first line time = " << metadata.first_line_time << "\n";

        for (int i = 0; i < n_points; i++) {
            OrbitStateVector osv = {};
            osv.time = first_osv_time + i * metadata.orbit_interval;
            osv.x_pos = platform_position_data.ReadD22(387 + i * 132);
            osv.y_pos = platform_position_data.ReadD22(387 + 22 + i * 132);
            osv.z_pos = platform_position_data.ReadD22(387 + 44 + i * 132);
            osv.x_vel = platform_position_data.ReadD22(453 + i * 132);
            osv.y_vel = platform_position_data.ReadD22(453 + 22 + i * 132);
            osv.z_vel = platform_position_data.ReadD22(453 + 44 + i * 132);

            metadata.orbit_state_vectors.push_back(osv);
        }

        double acc_vel = 0.0;
        for (const auto& osv : metadata.orbit_state_vectors) {
            acc_vel += sqrt(osv.x_vel * osv.x_vel + osv.y_vel * osv.y_vel + osv.z_vel * osv.z_vel);
        }
        metadata.platform_velocity = acc_vel / metadata.orbit_state_vectors.size();

        // TODO(priit) - better azimuth spacing estimate?
        metadata.azimuth_spacing = (0.88 * metadata.platform_velocity) / metadata.pulse_repetition_frequency;
    }

    uint32_t min_slant_range = UINT32_MAX;
    uint32_t prev_slant_range = UINT32_MAX;
    std::vector<uint32_t> slant_ranges;
    std::vector<uint32_t> offsets;
    for (int i = 0; i < azimuth_size; i++) {
        auto r = GetSignalDataRecord(i, record_length);
        const uint32_t slant_range = r.ReadB4(117);
        if (slant_range < min_slant_range) {
            min_slant_range = slant_range;
        }
        if (slant_range != prev_slant_range) {
            fmt::print("Azimuth idx = {} slant range = {}\n", i, slant_range);
            prev_slant_range = slant_range;
        }

        slant_ranges.push_back(slant_range);
    }

    metadata.slant_range_first_sample = min_slant_range;

    uint32_t prev_offset = UINT32_MAX;
    offsets.resize(slant_ranges.size());
    for (size_t i = 0; i < slant_ranges.size(); i++) {
        const uint32_t slant_range = slant_ranges[i];
        uint32_t offset = 0;
        if (slant_range != min_slant_range) {
            offset = std::round((slant_range - min_slant_range) / metadata.range_spacing);
        }
        if (offset != prev_offset) {
            fmt::print("Azimuth idx = {} range left offset = {}\n", i, offset);
            prev_offset = offset;
        }
        offsets[i] = offset;
    }

    metadata.left_range_offsets = offsets;

    const uint32_t max_offset = *std::max_element(offsets.begin(), offsets.end());

    const int range_size = echo_samples + max_offset;
    fmt::print("echo samples = {} max offset = {} range size = {}\n", echo_samples, max_offset, range_size);
    metadata.echo_samples.resize(azimuth_size, echo_samples);
    const int range_padded = range_size + metadata.chirp.n_samples;

    auto fft_sizes = GetOptimalFFTSizes();

    auto it_range = fft_sizes.upper_bound(range_padded);
    const int range_stride = it_range->first;
    fmt::print("chirp padding = {} range FFT padding = {}\n", metadata.chirp.n_samples, range_stride);
    fmt::print("Range FFT size = 2^{} * 3^{} * 5^{} * 7^{}\n", it_range->second[0], it_range->second[1],
               it_range->second[2], it_range->second[3]);

    auto it_azimuth = fft_sizes.upper_bound(azimuth_size);
    int azimuth_stride = it_azimuth->first;

    fmt::print("azimuth size = {} azimuth stride = {}\n", azimuth_size, azimuth_stride);
    fmt::print("Azimuth FFT size = 2^{} * 3^{} * 5^{} * 7^{}\n", it_azimuth->second[0], it_azimuth->second[1],
               it_azimuth->second[2], it_azimuth->second[3]);

    img.Init(range_size, azimuth_size, range_stride, azimuth_stride);

    IQ<float>* data = img.Get();

    int y = 0;
    for (; y < azimuth_size; y++) {
        int x = 0;
        const int near_offset = offsets[y];

        for (; x < near_offset; x++) {
            data[x + y * range_stride] = {};
        }

        const char* iq8 = image_file_.SignalDataRow(y, record_length) + SIGNAL_OFFSET;
        for (int i = 0; i < echo_samples; i++) {
            auto& d = data[x + y * range_stride];
            d.i = iq8[2 * i];
            d.q = iq8[2 * i + 1];
            x++;
        }
        for (; x < range_stride; x++) {
            data[x + y * range_stride] = {};
        }
    }

    for (; y < azimuth_stride; y++) {
        for (int x = 0; x < range_stride; x++) {
            data[x + y * range_stride] = {};
        }
    }
}

}  // namespace palsar