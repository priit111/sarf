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

struct OrbitStateVector {
    double time;
    double x_pos;
    double y_pos;
    double z_pos;
    double x_vel;
    double y_vel;
    double z_vel;
};

inline OrbitStateVector InterpolateOrbit(const std::vector<OrbitStateVector>& osv, double time) {

    OrbitStateVector r = {};
    r.time = time;
    const int n = osv.size();
    const double first_time = osv.front().time;
    for (int i = 0; i < n; i++) {
        double mult = 1;
        for (int j = 0; j < n; j++) {
            if (i == j) continue;

            double xj = osv.at(j).time;
            double xi = osv.at(i).time;
            mult *= (time - xj) / (xi - xj);
        }

        r.x_pos += mult * osv[i].x_pos;
        r.y_pos += mult * osv[i].y_pos;
        r.z_pos += mult * osv[i].z_pos;
        r.x_vel += mult * osv[i].x_vel;
        r.y_vel += mult * osv[i].y_vel;
        r.z_vel += mult * osv[i].z_vel;
    }

    return r;
}