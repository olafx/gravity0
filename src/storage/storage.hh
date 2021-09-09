#pragma once
#include <string>
#include <H5Cpp.h>

using namespace H5;

struct Storage
{
    H5File file;

    void new_dataset(void *buffer, hsize_t size, uint64_t name)
    {   hsize_t shape[2] {size, 2};
        DataSpace dataspace {2, shape};
        FloatType datatype {PredType::NATIVE_DOUBLE};
        DataSet dataset {file.createDataSet(std::to_string(name), datatype, dataspace)};
        dataset.write(buffer, PredType::NATIVE_DOUBLE);
    }

    /*
    void new_size(uint64_t size)
    {   DataSpace attr_dataspace {H5S_SCALAR};
        IntType attr_datatype {PredType::NATIVE_UINT64};
        Attribute attr {file.createAttribute("size", attr_datatype, attr_dataspace)};
        attr.write(PredType::NATIVE_UINT64, &size);
    }

    void write_size(uint64_t size)
    {   Attribute attr {file.openAttribute("size")};
        attr.write(PredType::NATIVE_UINT64, &size);
    }

    uint64_t read_size()
    {   Attribute attr {file.openAttribute("size")};
        uint64_t buffer;
        attr.read(PredType::NATIVE_UINT64, &buffer);
        return buffer;
    }
    */

    void read_dataset(void *buffer, uint64_t name)
    {   DataSet dataset {file.openDataSet(std::to_string(name))};
        dataset.read(buffer, PredType::NATIVE_DOUBLE);
    }

    hsize_t read_dataset_size(uint64_t name)
    {   DataSet dataset {file.openDataSet(std::to_string(name))};
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
