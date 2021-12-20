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
#include <algorithm>

namespace Storage
{

using namespace H5;

//  time dependent 3D n-body storage directly via HDF5
//  TODO
//    - Write SZip compressed.

struct N_Body_h5
{
    H5File      file;
    std::size_t current_group_name;
    Group       current_group;
    DataSet     current_pos_set;
    DataSet     current_vel_set;
    DataSet     current_time_set;
    hsize_t     time_group_size;
    hsize_t     time_chunk_size;


    explicit N_Body_h5(std::string name, hsize_t time_group_size = 256, hsize_t time_chunk_size = 16)
        : file               {std::move(name) + ".h5", H5F_ACC_RDWR},
          current_group_name {file.openGroup("/").getNumObjs() - 2},
          current_group      {file.openGroup(std::to_string(current_group_name))},
          current_pos_set    {current_group.openDataSet("pos")},
          current_vel_set    {current_group.openDataSet("vel")},
          current_time_set   {current_group.openDataSet("time")},
          time_group_size    {std::move(time_group_size)},
          time_chunk_size    {std::move(time_chunk_size)}
    {}


    hsize_t n_objects(std::size_t group_name = 0) const
    {
        Group     group   = file.openGroup(std::to_string(group_name));
        DataSet   pos_set = group.openDataSet("pos");
        DataSpace space   = pos_set.getSpace();

        hsize_t shape[3];
        space.getSimpleExtentDims(shape);

        return shape[1];
    }


    void read(double *const pos_buffer,
              const std::size_t group_name = 0, const hsize_t n_time_steps = 1) const
    {
        Group     group   = file.openGroup(std::to_string(group_name));
        DataSet   pos_set = group.openDataSet("pos");
        DataSpace space   = pos_set.getSpace();

        hsize_t shape[3];
        space.getSimpleExtentDims(shape);

        hsize_t   mem_shape[3] {std::min(n_time_steps, shape[0]), shape[1], shape[2]}, mem_offset[3] {};
        DataSpace mem_space    {3, mem_shape};

        space.selectHyperslab(H5S_SELECT_SET, mem_shape, mem_offset);
        pos_set.read(pos_buffer, PredType::NATIVE_DOUBLE, mem_space, space);

        if (n_time_steps > shape[0])
        {   hsize_t offset_size = shape[0] * shape[1] * shape[2];
            read(pos_buffer + offset_size, group_name + 1, n_time_steps - shape[0]);
        }
    }


    void read(double *const pos_buffer, double *const vel_buffer,
              const std::size_t group_name = 0, const hsize_t n_time_steps = 1) const
    {
        Group     group   = file.openGroup(std::to_string(group_name));
        DataSet   pos_set = group.openDataSet("pos");
        DataSet   vel_set = group.openDataSet("vel");
        DataSpace space   = pos_set.getSpace();

        hsize_t shape[3];
        space.getSimpleExtentDims(shape);

        hsize_t   mem_shape[3]  {std::min(n_time_steps, shape[0]), shape[1], shape[2]};
        hsize_t   mem_offset[3] {};
        DataSpace mem_space     {3, mem_shape};

        space.selectHyperslab(H5S_SELECT_SET, mem_shape, mem_offset);
        pos_set.read(pos_buffer, PredType::NATIVE_DOUBLE, mem_space, space);
        vel_set.read(vel_buffer, PredType::NATIVE_DOUBLE, mem_space, space);

        if (n_time_steps > shape[0])
        {   hsize_t size = shape[0] * shape[1] * shape[2];
            read(pos_buffer + size, vel_buffer + size, group_name + 1, n_time_steps - shape[0]);
        }
    }


    void write(const double *const pos_buffer, const double *const vel_buffer, const double time)
    {
        DataSpace current_pos_space  = current_pos_set .getSpace();
        DataSpace current_time_space = current_time_set.getSpace();
        hsize_t current_pos_shape[3];
        current_pos_space.getSimpleExtentDims(current_pos_shape);

        if (current_pos_shape[0] == time_group_size)
        {   
            current_group = file.createGroup(std::to_string(++current_group_name));

            current_pos_shape[0] = 0;
            hsize_t max_pos_shape[3] {H5S_UNLIMITED, current_pos_shape[1], 3};
            current_pos_space.setExtentSimple(3, current_pos_shape, max_pos_shape);

            DSetCreatPropList properties;
            hsize_t chunk_shape[3] {time_chunk_size, current_pos_shape[1], 3};
            properties.setChunk(3, chunk_shape);
            constexpr double filler = 0;
            properties.setFillValue(PredType::NATIVE_DOUBLE, &filler);

            current_pos_set = current_group.createDataSet("pos", PredType::NATIVE_DOUBLE, current_pos_space, properties);
            current_vel_set = current_group.createDataSet("vel", PredType::NATIVE_DOUBLE, current_pos_space, properties);

            current_time_space.setExtentSimple(1, current_pos_shape, max_pos_shape);

            properties.setChunk(1, chunk_shape);

            current_time_set = current_group.createDataSet("time", PredType::NATIVE_DOUBLE, current_time_space, properties);
        }

        current_pos_shape[0]++;

        current_pos_set. extend(current_pos_shape);
        current_vel_set. extend(current_pos_shape);
        current_time_set.extend(current_pos_shape);

        current_pos_space  = current_pos_set .getSpace();
        current_time_space = current_time_set.getSpace();

        hsize_t slab_shape [3] {1, current_pos_shape[1], 3};
        hsize_t slab_offset[3] {current_pos_shape[0] - 1, 0, 0};
        current_pos_space .selectHyperslab(H5S_SELECT_SET, slab_shape, slab_offset);
        current_time_space.selectHyperslab(H5S_SELECT_SET, slab_shape, slab_offset);

        DataSpace pos_mem_space  {3, slab_shape};
        DataSpace time_mem_space {1, slab_shape};

        current_pos_set .write(pos_buffer, PredType::NATIVE_DOUBLE, pos_mem_space,  current_pos_space);
        current_vel_set .write(vel_buffer, PredType::NATIVE_DOUBLE, pos_mem_space,  current_pos_space);
        current_time_set.write(&time     , PredType::NATIVE_DOUBLE, time_mem_space, current_time_space);
    }


    void write(const double *const pos_buffer, const double time)
    {
        DataSpace current_pos_space  = current_pos_set .getSpace();
        DataSpace current_time_space = current_time_set.getSpace();
        hsize_t current_pos_shape[3];
        current_pos_space.getSimpleExtentDims(current_pos_shape);

        if (current_pos_shape[0] == time_group_size)
        {   
            current_group = file.createGroup(std::to_string(++current_group_name));

            current_pos_shape[0] = 0;
            hsize_t max_pos_shape[3] {H5S_UNLIMITED, current_pos_shape[1], 3};
            current_pos_space.setExtentSimple(3, current_pos_shape, max_pos_shape);

            DSetCreatPropList properties;
            hsize_t chunk_shape[3] {time_chunk_size, current_pos_shape[1], 3};
            properties.setChunk(3, chunk_shape);
            constexpr double filler = 0;
            properties.setFillValue(PredType::NATIVE_DOUBLE, &filler);

            current_pos_set = current_group.createDataSet("pos", PredType::NATIVE_DOUBLE, current_pos_space, properties);

            current_time_space.setExtentSimple(1, current_pos_shape, max_pos_shape);

            properties.setChunk(1, chunk_shape);

            current_time_set = current_group.createDataSet("time", PredType::NATIVE_DOUBLE, current_time_space, properties);
        }

        current_pos_shape[0]++;

        current_pos_set. extend(current_pos_shape);
        current_time_set.extend(current_pos_shape);

        current_pos_space  = current_pos_set .getSpace();
        current_time_space = current_time_set.getSpace();

        hsize_t slab_shape [3] {1, current_pos_shape[1], 3};
        hsize_t slab_offset[3] {current_pos_shape[0] - 1, 0, 0};
        current_pos_space .selectHyperslab(H5S_SELECT_SET, slab_shape, slab_offset);
        current_time_space.selectHyperslab(H5S_SELECT_SET, slab_shape, slab_offset);

        DataSpace pos_mem_space  {3, slab_shape};
        DataSpace time_mem_space {1, slab_shape};

        current_pos_set .write(pos_buffer, PredType::NATIVE_DOUBLE, pos_mem_space,  current_pos_space);
        current_time_set.write(&time     , PredType::NATIVE_DOUBLE, time_mem_space, current_time_space);
    }
};

}
