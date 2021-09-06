#pragma once
#include <H5Cpp.h>
#include <string>

using namespace H5;

struct Storage
{
    H5File file;

    Storage() : file {"gravity0.h5", H5F_ACC_TRUNC}
    {
    }

    void add(double *s, size_t n, size_t t)
    {   hsize_t dims[2] {n, 2};
        DataSpace dataspace {2, dims};
        FloatType datatype {PredType::NATIVE_DOUBLE};
        DataSet dataset {file.createDataSet(std::to_string(t), datatype, dataspace)};
        dataset.write(s, PredType::NATIVE_DOUBLE);
    }
};
