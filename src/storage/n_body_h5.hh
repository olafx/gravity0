/*
    MIT License

    Copyright (c) 2021 Olaf Willocx

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#pragma once
#include <H5Cpp.h>
#include <string>

namespace Storage
{

using namespace H5;

//  time dependent 3D n-body storage directly via HDF5
//  TODO
//    - support what the Python initial condition writer now writes

struct N_Body_h5
{
    H5File file;

    N_Body_h5(const std::string& name)
        : file {name, H5F_ACC_RDWR}
    {
    }

    void write(const double *const buffer, const hsize_t size, const std::string& name) const
    {   hsize_t shape[2] {size, 3};
        DataSpace dataspace {2, shape};
        FloatType datatype {PredType::NATIVE_DOUBLE};
        DataSet dataset {file.createDataSet(name, datatype, dataspace)};
        dataset.write(buffer, PredType::NATIVE_DOUBLE);
    }

    void read(double *const buffer, const std::string& name) const
    {   DataSet dataset {file.openDataSet(name)};
        dataset.read(buffer, PredType::NATIVE_DOUBLE);
    }

    hsize_t read_size(const std::string& name) const
    {   DataSet dataset {file.openDataSet(name)};
        DataSpace dataspace {dataset.getSpace()};
        hsize_t shape[2];
        dataspace.getSimpleExtentDims(shape);
        return shape[0] / 2;
    }

    hsize_t size() const
    {   Group root {file.openGroup("/")};
        return root.getNumObjs();
    }
};

}
