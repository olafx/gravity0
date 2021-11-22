#pragma once
#include <H5Cpp.h>
#include <string>

namespace Storage
{

using namespace H5;

//  time dependent n-body storage directly via HDF5
//  TODO
//    - store a few large 3D arrays rather than tons of 2D arrays
//    - store times
//    - format attributes in root

struct n_hdf5
{
    H5File file;

    n_hdf5(const std::string& name)
        : file {name, H5F_ACC_RDWR}
    {
    }

    void new_dataset(const double *const buffer, const hsize_t size, const std::string& name) const
    {   hsize_t shape[2] {size, 3};
        DataSpace dataspace {2, shape};
        FloatType datatype {PredType::NATIVE_DOUBLE};
        DataSet dataset {file.createDataSet(name, datatype, dataspace)};
        dataset.write(buffer, PredType::NATIVE_DOUBLE);
    }

    void read_dataset(double *const buffer, const std::string& name) const
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

}
