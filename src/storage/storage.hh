#pragma once
#include <string>
#include <H5Cpp.h>

//  This homebrew HDF5 usage is temporary, will update to using Xdmf3 standard through VTK with SZip time compression.

using namespace H5;

struct Storage
{
    H5File file;

    void new_dataset(const void *buffer, hsize_t size, const std::string& name) const
    {   hsize_t shape[2] {size, 3};
        DataSpace dataspace {2, shape};
        FloatType datatype {PredType::NATIVE_DOUBLE};
        DataSet dataset {file.createDataSet(name, datatype, dataspace)};
        dataset.write(buffer, PredType::NATIVE_DOUBLE);
    }

    void read_dataset(void *buffer, const std::string& name) const
    {   DataSet dataset {file.openDataSet(name)};
        dataset.read(buffer, PredType::NATIVE_DOUBLE);
    }

    hsize_t read_dataset_size(const std::string& name) const
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
