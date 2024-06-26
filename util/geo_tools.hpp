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

#include <cmath>

struct GeoPosLLH {
    double latitude;
    double longitude;
    double height;
};

struct GeoPos3D {
    double x;
    double y;
    double z;
};

namespace WGS84 {                                        // NOLINT
constexpr double A = 6378137.0;                          // m
constexpr double B = 6356752.3142451794975639665996337;  // 6356752.31424518; // m
constexpr double FLAT_EARTH_COEF = 1.0 / ((A - B) / A);  // 298.257223563;
constexpr double E2 = 2.0 / FLAT_EARTH_COEF - 1.0 / (FLAT_EARTH_COEF * FLAT_EARTH_COEF);
constexpr double EP2 = E2 / (1 - E2);
}  // namespace WGS84

constexpr double DTOR = M_PI / 180.0;
constexpr double RTOD = 180 / M_PI;

inline double Distance(GeoPos3D xyz1, GeoPos3D xyz2) {
    double dx = xyz1.x - xyz2.x;
    double dy = xyz1.y - xyz2.y;
    double dz = xyz1.z - xyz2.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

inline GeoPos3D Geo2xyzWgs84(double latitude, double longitude, double altitude) {
    double const lat = latitude * DTOR;
    double const lon = longitude * DTOR;

    double sin_lat;
    double cos_lat;
    sincos(lat, &sin_lat, &cos_lat);

    double const sinLat = sin_lat;

    double const N = (WGS84::A / sqrt(1.0 - WGS84::E2 * sinLat * sinLat));
    double const NcosLat = (N + altitude) * cos_lat;

    double sin_lon;
    double cos_lon;
    sincos(lon, &sin_lon, &cos_lon);

    GeoPos3D r = {};
    r.x = NcosLat * cos_lon;
    r.y = NcosLat * sin_lon;
    r.z = (N + altitude - WGS84::E2 * N) * sinLat;
    return r;
}

inline GeoPosLLH xyz2geoWGS84(GeoPos3D xyz) {
    double x = xyz.x;
    double y = xyz.y;
    double z = xyz.z;
    double s = sqrt(x * x + y * y);
    double theta = atan(z * WGS84::A / (s * WGS84::B));

    GeoPosLLH llh = {};
    llh.longitude = atan(y / x) * RTOD;

    if (llh.longitude < 0.0 && y >= 0.0) {
        llh.longitude += 180.0;
    } else if (llh.longitude > 0.0 && y < 0.0) {
        llh.longitude -= 180.0;
    }

    llh.latitude =
        atan((z + WGS84::EP2 * WGS84::B * pow(sin(theta), 3)) / (s - WGS84::E2 * WGS84::A * pow(cos(theta), 3))) * RTOD;
    return llh;
}

struct Velocity3D {
    double x;
    double y;
    double z;
};

GeoPos3D RangeDopplerGeoLocate(Velocity3D sensor_vel, GeoPos3D sensor_pos, GeoPos3D init_guess, double slant_range);

double CalcIncidenceAngle(GeoPos3D earth_point, GeoPos3D sensor_pos);