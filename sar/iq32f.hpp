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


template<class T>
struct IQ
{
    T i;
    T q;
};

template<class T>
inline IQ<T> complex_mul(IQ<T> a, IQ<T> b)
{
    T real = a.i * b.i  - a.q * b.q;
    T imag = a.i * b.q + a.q * b.i;
    return {real, imag};
}


