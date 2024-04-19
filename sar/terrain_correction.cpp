#include "terrain_correction.hpp"

#include <algorithm>
#include <thread>
#include <atomic>
#include <vector>

#include "fmt/format.h"

#include "util/geo_tools.hpp"
#include "util/util.hpp"
#include "util/img_out.hpp"


namespace {
struct Vec3d {
    double x;
    double y;
    double z;
};
double getDopplerFrequency(Vec3d earthPoint, Vec3d sensorPosition, Vec3d sensorVelocity, double wavelength) {
    const auto xDiff = earthPoint.x - sensorPosition.x;
    const auto yDiff = earthPoint.y - sensorPosition.y;
    const auto zDiff = earthPoint.z - sensorPosition.z;
    const auto distance = sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff);

    return 2.0 * (sensorVelocity.x * xDiff + sensorVelocity.y * yDiff + sensorVelocity.z * zDiff) /
           (distance * wavelength);
}

struct Args {
    double lon_start;
    double lat_start;
    double pixel_spacing_y;
    double pixel_spacing_x;
    double line_time_interval;
    double wavelength;

    double range_spacing;
    double slant_range_first_sample;
    int in_azimuth_size;
    const OrbitStateVector* osv;
    int n_osv;
    const Vec3d* sat_pos_az;
    const Vec3d* sat_vel_az;
};

inline double GetEarthPointZeroDopplerTime(double line_time_interval, double wavelength, Vec3d earth_point,
                                               int n_azimuth, const Vec3d* sensor_position,
                                               const Vec3d* sensor_velocity) {
    const double first_line_utc = 0.0;
    // binary search is used in finding the zero doppler time
    int lower_bound = 0;
    int upper_bound = n_azimuth - 1;
    auto lower_bound_freq =
        getDopplerFrequency(earth_point, sensor_position[lower_bound], sensor_velocity[lower_bound], wavelength);
    auto upper_bound_freq =
        getDopplerFrequency(earth_point, sensor_position[upper_bound], sensor_velocity[upper_bound], wavelength);

    if (std::abs(lower_bound_freq) < 1.0) {
        return first_line_utc + lower_bound * line_time_interval;
    } else if (std::abs(upper_bound_freq) < 1.0) {
        return first_line_utc + upper_bound * line_time_interval;
    } else if (lower_bound_freq * upper_bound_freq > 0.0) {
        return NAN;
    }

    // start binary search
    double mid_freq;
    while (upper_bound - lower_bound > 1) {
        const auto mid = (int)((lower_bound + upper_bound) / 2.0);
        mid_freq = sensor_velocity[mid].x * (earth_point.x - sensor_position[mid].x) +
                   sensor_velocity[mid].y * (earth_point.y - sensor_position[mid].y) +
                   sensor_velocity[mid].z * (earth_point.z - sensor_position[mid].z);

        if (mid_freq * lower_bound_freq > 0.0) {
            lower_bound = mid;
            lower_bound_freq = mid_freq;
        } else if (mid_freq * upper_bound_freq > 0.0) {
            upper_bound = mid;
            upper_bound_freq = mid_freq;
        } else if (mid_freq == 0.0) {
            return first_line_utc + mid * line_time_interval;
        }
    }

    const auto y0 =
        lower_bound - lower_bound_freq * (upper_bound - lower_bound) / (upper_bound_freq - lower_bound_freq);
    return first_line_utc + y0 * line_time_interval;
}

inline Vec3d GetPosition(double time, const OrbitStateVector* vectors, int n_osv) {
    int i0 = 0;
    int iN = n_osv - 1;

    Vec3d result{0, 0, 0};
    for (int i = 0; i <= n_osv; ++i) {
        auto const orbI = vectors[i];

        double weight = 1;
        for (int j = i0; j <= iN; ++j) {
            if (j != i) {
                double const time2 = vectors[j].time;
                weight *= (time - time2) / (orbI.time - time2);
            }
        }
        result.x += weight * orbI.x_pos;
        result.y += weight * orbI.y_pos;
        result.z += weight * orbI.z_pos;
    }
    return result;
}


bool myIsnan(double v) {
    std::uint64_t i;
    memcpy(&i, &v, 8);
    return ((i&0x7ff0000000000000)==0x7ff0000000000000)&&(i&0xfffffffffffff);
}

void RunTC(const IQ<float>* in, float* out, int out_x_size, int y_start, int y_end, int in_range_size, int in_range_stride, Args args) {
    for (int y = y_start; y < y_end; y++) {
        for (int x = 0; x < out_x_size; x++) {

            const double lat = args.lat_start + y * args.pixel_spacing_y + 0.5 * args.pixel_spacing_y;
            const double lon = args.lon_start + x * args.pixel_spacing_x + 0.5 * args.pixel_spacing_x;
            auto tmp = Geo2xyzWgs84(lat, lon, 0);
            Vec3d earth_point = {.x = tmp.x, .y = tmp.y, .z =tmp.z};
            double az_time =
                GetEarthPointZeroDopplerTime(args.line_time_interval, args.wavelength, earth_point,
                                                 args.in_azimuth_size, args.sat_pos_az, args.sat_vel_az);
            const int out_idx = x + y * out_x_size;
            if(myIsnan(az_time))
            {
                out[out_idx] = 0.0f;
                continue;
            }
            double az_idx = az_time / args.line_time_interval;

            if (az_idx < 0.0 || az_idx >= args.in_azimuth_size) {
                out[out_idx] = 0.0f;
                continue;
            }

            // slant range


            //Vec3d sat_pos = GetPosition(az_time, args.osv, args.n_osv);

            Vec3d sat_pos = args.sat_pos_az[static_cast<int>(std::round(az_idx))];

            // range idx
            double dx = earth_point.x - sat_pos.x;
            double dy = earth_point.y - sat_pos.y;
            double dz = earth_point.z - sat_pos.z;
            double slant_range = sqrt(dx * dx + dy * dy + dz * dz);
            double rg_idx = (slant_range - args.slant_range_first_sample) / args.range_spacing;

            if(rg_idx < 0 || rg_idx >= in_range_size)
            {
                out[out_idx] = 0.0f;
                continue;
            }

            int x_i = (int)rg_idx;
            int y_i = (int)az_idx;

            float val0 = ToIntens(in[y_i * in_range_stride + x_i]);
            float val1 = ToIntens(in[y_i * in_range_stride + x_i + 1]);
            float val2 = ToIntens(in[(y_i + 1) * in_range_stride + x_i]);
            float val3 = ToIntens(in[(y_i + 1) * in_range_stride + x_i + 1]);

            float interp_x = (float)rg_idx - x_i;
            float interp_y = (float)az_idx - y_i;

            float int_r0 = val0 + interp_x * (val1 - val0);
            float int_r1 = val2 + interp_x * (val3 - val2);
            float data = int_r0 + interp_y * (int_r1 - int_r0);

            // write result

            // out[out_idx] = 10.0f * logf(data);
            out[out_idx] = data;  // * 0.0f + 1.0f;
        }
    }
}
}  // namespace

void RangeDopplerTerrainCorrection(const SARImage& img, const SARMetadata& metadata) {
    int range_size = img.XSize();
    int range_stride = img.XStride();
    int azimuth_size = img.YSize();
    auto UL = xyz2geoWGS84(RangeDopplerGeoLocate(metadata, 0, 0));
    auto UR = xyz2geoWGS84(RangeDopplerGeoLocate(metadata, range_size - 1, 0));
    auto LL = xyz2geoWGS84(RangeDopplerGeoLocate(metadata, 0, azimuth_size - 1));
    auto LR = xyz2geoWGS84(RangeDopplerGeoLocate(metadata, range_size - 1, azimuth_size - 1));

    double lon_max = std::max({UL.longitude, UR.longitude, LL.longitude, LR.longitude});
    double lon_min = std::min({UL.longitude, UR.longitude, LL.longitude, LR.longitude});
    double lat_max = std::max({UL.latitude, UR.latitude, LL.latitude, LR.latitude});
    double lat_min = std::min({UL.latitude, UR.latitude, LL.latitude, LR.latitude});

    fmt::print("geobox = {} {} {} {}\n", lat_max, lat_min, lon_max, lon_min);

    fmt::print("az spacing = {} rg spacing = {}\n", metadata.azimuth_spacing, metadata.range_spacing);
    double pixel_spacing_in_meter = std::max({metadata.azimuth_spacing, metadata.range_spacing, 1.0});
    double pixel_spacing_in_degree = pixel_spacing_in_meter / WGS84::A * (180.0 / M_PI);
    double dif_lat = std::fabs(lat_max - lat_min);
    double dif_lon = std::fabs(lon_max - lon_min);

    int y_size = dif_lat / pixel_spacing_in_degree;
    int x_size = dif_lon / pixel_spacing_in_degree;

    size_t total = (size_t)x_size * y_size;
    fmt::print("Terrain correction output dimensions = {} {}, pixel spacing = {} m\n", x_size, y_size,
               pixel_spacing_in_meter);


    float* out = new float[total];
    fmt::print("out allocated...\n");
    Args args = {};
    args.lat_start = lat_max;
    args.lon_start = lon_min;
    args.pixel_spacing_x = pixel_spacing_in_degree;
    args.pixel_spacing_y = -pixel_spacing_in_degree;

    std::vector<Vec3d> pos;
    std::vector<Vec3d> vel;
    for(int i = 0; i < azimuth_size; i++)
    {
        auto osv = InterpolateOrbit(metadata.orbit_state_vectors, CalcAzimuthTime(metadata, i));
        pos.push_back({osv.x_pos, osv.y_pos, osv.z_pos});
        vel.push_back({osv.x_vel, osv.y_vel, osv.z_vel});
    }

    auto osv = metadata.orbit_state_vectors;

    osv.clear();

    for(auto& e : metadata.orbit_state_vectors)
    {
        if(e.time >= -180 || e.time <= 180)
        {
            osv.push_back(e);
        }
    }

    args.line_time_interval = metadata.line_time_interval;
    args.wavelength = metadata.wavelength;
    args.slant_range_first_sample = metadata.slant_range_first_sample;
    args.range_spacing = metadata.range_spacing;
    args.sat_pos_az = pos.data();
    args.sat_vel_az = vel.data();
    args.osv = osv.data();
    args.n_osv = osv.size();
    args.in_azimuth_size = azimuth_size;

    auto tc_start = TimeStart();
    const int n_threads =  std::thread::hardware_concurrency();
    fmt::print("TC N THREADS = {}\n", n_threads);
    std::vector<std::thread> thread_vec;
    uint32_t y_step = (y_size / n_threads) + 1;

    int start_y = 0;
    for(int i = 0; i < n_threads; i++)
    {
        start_y = i * y_step;
        int end_y = start_y + y_step;
        end_y = std::min(end_y, y_size);
        std::thread t(RunTC, img.Get(), out, x_size, start_y, end_y, range_size, range_stride, args);

        thread_vec.push_back(std::move(t));
    }


    for(auto& t : thread_vec)
    {
        t.join();
    }

    TimeStop(tc_start, "TC time");

    double gt[6] = {};
    gt[0] = lon_min;
    gt[1] = pixel_spacing_in_degree;
    gt[3] = lat_max;
    gt[5] = -pixel_spacing_in_degree;

    WriteOut("/tmp/tc.tiff", out, x_size, y_size, gt);
}