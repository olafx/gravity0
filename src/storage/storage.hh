#pragma once
#include <string>
#include <H5Cpp.h>

using namespace H5;

struct Storage
{
    H5File file;

    void new_dataset(void *buffer, hsize_t size, std::string name)
    {   hsize_t shape[2] {size, 2};
        DataSpace dataspace {2, shape};
        FloatType datatype {PredType::NATIVE_DOUBLE};
        DataSet dataset {file.createDataSet(name, datatype, dataspace)};
        dataset.write(buffer, PredType::NATIVE_DOUBLE);
    }

    void read_dataset(void *buffer, std::string name)
    {   DataSet dataset {file.openDataSet(name)};
        dataset.read(buffer, PredType::NATIVE_DOUBLE);
    }

    hsize_t read_dataset_size(std::string name)
    {   DataSet dataset {file.openDataSet(name)};
        DataSpace dataspace {dataset.getSpace()};
        hsize_t shape[2];
        dataspace.getSimpleExtentDims(shape);
        return shape[0];
    }

    hsize_t size()
    {   Group root {file.openGroup("/")};
        return root.getNumObjs();
    }
};
