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

/**
* ENVISAT and ERS ASAR instrument focusser for QA4EO activity (c) by CGI Estonia AS
*
* ENVISAT and ERS ASAR instrument focusser for QA4EO activity is licensed under a
* Creative Commons Attribution-ShareAlike 4.0 International License.
*
* You should have received a copy of the license along with this
* work. If not, see http://creativecommons.org/licenses/by-sa/4.0/
*/
#include "geo_tools.hpp"

#include "Eigen/Dense"

namespace {
Eigen::Vector3d ToVec3(GeoPos3D geo_pos) { return Eigen::Vector3d{geo_pos.x, geo_pos.y, geo_pos.z}; }

// Eigen::Vector3d ToVec3(Velocity3D vel) { return Eigen::Vector3d{vel.x, vel.y, vel.z}; }
}  // namespace

GeoPos3D RangeDopplerGeoLocate(Velocity3D velocity, GeoPos3D position, GeoPos3D init_guess, double slant_range) {
   constexpr double a = 6378137.0;
   constexpr double b = 6356752.314245;


   Eigen::Vector3d result = {init_guess.x, init_guess.y, init_guess.z};
   Eigen::Vector3d sensor_vel = {velocity.x, velocity.y, velocity.z};
   Eigen::Vector3d sensor_pos = {position.x, position.y, position.z};

   int iter = 0;
   while (1) {
       Eigen::Matrix3d A(3, 3);

       Eigen::Vector3d pos_vec = result - sensor_pos;

       /*
        * x,y,z => point on the ground we are trying to find
        *
        *
        * 1. [(x - x_sat), (y - y_sat), (z - z_sat)] dot [vel_x_sat, vel_y_sat, vel_z_sat] = 0
        *  This is from the Doppler equation, if satellite velocity is perpendicular to ground point, then the Doppler
        * effect is zero
        *
        *
        * 2. |pos_vec|^2 = R^2
        *
        * distance between satellite and ground point must be equal to the range distance
        *
        * 3. (x/a)^2 + (y/a)^2 + (z/b)^2 = 1
        *
        * The point must be on the WGS84 ellipsoid
        */

       // partial derivatives of the three equations x,y,z
       A.row(0) = sensor_vel;
       A.row(1) = 2.0 * pos_vec;  //
       A(2, 0) = (2 * result[0]) / (a * a);
       A(2, 1) = (2 * result[1]) / (a * a);
       A(2, 2) = (2 * result[2]) / (b * b);

       Eigen::Vector3d B;
       B[0] = sensor_vel.dot(pos_vec);
       B[1] = pos_vec.dot(pos_vec) - slant_range * slant_range;  // distance between
       B[2] = pow(result[0] / a, 2.0) + pow(result[1] / a, 2.0) + pow(result[2] / b, 2.0) - 1.0;

       // TODO better comments
       Eigen::Vector3d dx = A.fullPivHouseholderQr().solve(-B);

       result = result + dx;
       iter++;

       if (dx.norm() < 1e-6) {
           return {result[0], result[1], result[2]};
       }

       if (iter > 100) {
           throw std::runtime_error("Failed to converge for " + std::string(__FUNCTION__));
       }
   }
}

double CalcIncidenceAngle(GeoPos3D earth_point, GeoPos3D sensor_pos) {
   // find surface normal equation to wgs84
   double a_sq = WGS84::A * WGS84::A;
   double b_sq = WGS84::B * WGS84::B;

   // should be gradient of the ellipsoid equation, x^2/a^2 + y^2/a^2 + z^2 / b^2 = 1
   // -> partial derivative with x, y, z
   // -> e.g partial derivative by x would be 2x / a^2, no point multiplying by 2 since we only care about direction
   Eigen::Vector3d v1 = {earth_point.x / a_sq, earth_point.y / a_sq, earth_point.z / b_sq};

   Eigen::Vector3d v2 = ToVec3(sensor_pos) - ToVec3(earth_point);

   // angle between the ground point surface normal vector and vector from satellite to earth point should be the
   // incidence angle
   double cos_theta = v1.dot(v2) / (v1.norm() * v2.norm());

   // acos + rad2deg
   const double degree = (360.0 / (2 * M_PI)) * acos(cos_theta);

   return degree;
}