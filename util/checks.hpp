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


#include <cstdlib>
#include <cstdio>



inline void VerifyFunc(bool condition, const char* filename, const char* function, int line)
{
    if(!condition)
    {
        printf("SARF_CHECK fail! file = %s  function = %s() line = %d\n", filename, function, line);
        exit(1);
    }
}


#define SARF_CHECK(a) VerifyFunc(a, __FILE__, __FUNCTION__, __LINE__)