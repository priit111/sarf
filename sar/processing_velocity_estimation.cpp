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
#include "processing_velocity_estimation.h"

#include <sstream>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "fmt/format.h"
#include "sar_metadata.hpp"
#include "util/geo_tools.hpp"

namespace {
double SquareFitVr(const Eigen::VectorXd& xvals, const Eigen::VectorXd& yvals) {
    Eigen::MatrixXd A(xvals.size(), 2);

    A.setZero();

    for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++) {
        A(j, 1) = xvals(j) * xvals(j);
    }

    auto Q = A.fullPivHouseholderQr();
    auto result = Q.solve(yvals);

    return sqrt(result[1]);
}

double CalcDistance(OrbitStateVector pos, GeoPos3D xyz) {
    double dx = pos.x_pos - xyz.x;
    double dy = pos.y_pos - xyz.y;
    double dz = pos.z_pos - xyz.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

double EstimateVrAtRange(const SARMetadata& metadata, int az_idx, int range_idx) {
    double R0 = metadata.slant_range_first_sample + range_idx * metadata.range_spacing;
    const int aperture_size = CalcAperturePixels(metadata, range_idx) / 3;
    const int N = 9;

    auto center_osv = InterpolateOrbit(metadata.orbit_state_vectors, CalcAzimuthTime(metadata, az_idx));

    auto earth_point =
        RangeDopplerGeoLocate({center_osv.x_vel, center_osv.y_vel, center_osv.z_vel},
                              {center_osv.x_pos, center_osv.y_pos, center_osv.z_pos}, metadata.center_point, R0);

    auto test = xyz2geoWGS84(earth_point);
    fmt::print("mid point = {} {}\n", test.latitude, test.longitude);

    const int step = aperture_size / (N - 1);

    // calculate distance between center point and positions on the aperture
    // goal is to find data points from real orbit state vector for the hyperbolic range function
    // R^2(n) = R0^2 + Vr^2 * (n)
    // R - slant range across azimuth time points
    // R0 - slant range and closes point
    // Vr - processing / effective radar velocity
    // (n) - relative azimuth time
    Eigen::VectorXd y_vals(N);  // Vr
    Eigen::VectorXd x_vals(N);  // t
    for (int i = 0; i < N; i++) {
        const int aperture_az_idx = az_idx + (i - (N / 2)) * step;
        ;
        auto idx_time = CalcAzimuthTime(metadata, aperture_az_idx);
        auto osv = InterpolateOrbit(metadata.orbit_state_vectors, idx_time);
        double dx = osv.x_pos - earth_point.x;
        double dy = osv.y_pos - earth_point.y;
        double dz = osv.z_pos - earth_point.z;
        double R_square = dx * dx + dy * dy + dz * dz;
        double R0_square = R0 * R0;
        double dt = center_osv.time - osv.time;

        y_vals[i] = R_square - R0_square;
        x_vals[i] = dt;
    }

    // Now we have calculated data, where slant range to center point varies with azimuth time

    // mathematically the vectors now contain data points for the equation
    // y = ax^2 + c, the best data fit for a gives us an estimate for Vr^2
    // Now we have calculated data, where slant range to center point varies with azimuth time

    return SquareFitVr(x_vals, y_vals);
}

}  // namespace

double EstimateProcessingVelocity(const SARMetadata& metadata, int range_size, int azimuth_size) {
    const int az_idx = azimuth_size / 2;  // TODO How much does Vr vary in azimuth direction?
    const int rg_idx = (range_size - metadata.chirp.n_samples) / 2;

    double Vr = EstimateVrAtRange(metadata, az_idx, rg_idx);

    fmt::print("Estimated Vr = {}\n", Vr);

    return Vr;
}