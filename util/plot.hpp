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

#include <string>
#include <vector>

struct Scatter{
    std::string line_name;
    std::vector<double> x;
    std::vector<double> y;
};






struct PlotArgs{
    std::string out_path;
    std::string graph_name;
    std::string x_axis_title;
    std::string y_axis_title;
    std::vector<Scatter> data;
};

void Plot(const PlotArgs& graph);