#pragma once

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

#include <cstring>
#include <memory>

#include "iq32f.hpp"

class SARImage {
public:
    void Init(int x_size, int y_size, int x_stride, int y_stride) {
        x_size_ = x_size;
        y_size_ = y_size;
        x_stride_ = x_stride;
        y_stride_ = y_stride;
        data_ = std::unique_ptr<IQ<float>[]>(new IQ<float>[x_stride * y_stride]);
    }

    int XSize() const { return x_size_; }

    int XStride() const { return x_stride_; }

    int YSize() const { return y_size_; }

    int YStride() const { return y_stride_; }

    void ResizeX(int x_size) { x_size_ = x_size; }

    void ZeroFill() { memset(data_.get(), 0, ByteSize()); }

    size_t ByteSize() const { return static_cast<size_t>(x_stride_) * y_stride_ * sizeof(data_[0]); }

    IQ<float>* Get() { return data_.get(); }

    const IQ<float>* Get() const { return data_.get(); }

private:
    std::unique_ptr<IQ<float>[]> data_;
    int x_size_;
    int y_size_;
    int x_stride_;
    int y_stride_;
};