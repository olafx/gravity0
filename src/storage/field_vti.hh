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
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkAbstractArray.h>
#include <string>
#include <type_traits>
#include <utility>
#include <algorithm>
#include <fmt/format.h>

#include <iostream> // TODO temp

namespace Storage
{

//  time dependent 2D or 3D uniform grid storage via VTK image format (multiple files, time stored)
//  TODO
//    - Figure out how to read TimeValues (written by WriteNextTime) rather than TimeStep, which is an integer count.
//      Can then implement multi time step and custom time initial condition reading.

template <std::size_t time_steps_per_file = 256, std::size_t order = 2,
          typename = std::enable_if_t<order == 2 || order == 3>>
class Field_vti
{
    vtkXMLImageDataWriter *writer   = NULL;
    vtkXMLImageDataReader *reader   = NULL;
    vtkImageData          *image    = NULL;
    vtkDoubleArray        *rho      = NULL;
    vtkDoubleArray        *vel      = NULL;
    std::array<int, order> res;
    std::size_t time_count = 0;
    std::string name_no_suffix;


    void pre_read()
    {
        image = vtkImageData::SafeDownCast(reader->GetOutputDataObject(0));

        rho = vtkDoubleArray::SafeDownCast(image->GetPointData()->GetAbstractArray("Density"));
        vel = vtkDoubleArray::SafeDownCast(image->GetPointData()->GetAbstractArray("Velocity"));

        int extent[6];
        image->GetExtent(extent);

        res[0] = extent[1] + 1;
        res[1] = extent[3] + 1;
        if constexpr (order == 3)
            res[2] = extent[5] + 1;

        // reader->GetTimeStep();
        // reader->SetTimeStep();
        // reader->GetNumberOfTimeSteps();
    }


    void pre_write()
    {
        writer = vtkXMLImageDataWriter::New();
        image  = vtkImageData::New();
        rho    = vtkDoubleArray::New();
        vel    = vtkDoubleArray::New();

        rho->SetName("Density");
        rho->SetNumberOfComponents(1);

        vel->SetName("Velocity");
        vel->SetNumberOfComponents(3);

        image->GetPointData()->SetScalars(rho);
        image->GetPointData()->AddArray(vel);
        image->SetExtent(0, res[0] - 1, 0, res[1] - 1, 0, order == 3 ? res[2] - 1 : 0);

        writer->SetInputData(image);
        writer->SetNumberOfTimeSteps(time_steps_per_file);
    }


public:

    explicit Field_vti(std::string name)
        : name_no_suffix {std::move(name)}
    {
        reader = vtkXMLImageDataReader::New();
        reader->SetFileName((name_no_suffix + ".vti").c_str());
        reader->Update();
    }


    ~Field_vti()
    {
        if (writer != NULL)
        {
            writer->Stop();

            rho->Delete();
            if (vel != NULL)
                vel->Delete();
            image->Delete();
            writer->Delete();
        }
        if (reader != NULL)
        {
            reader->Delete();
        }
    }


    std::array<int, order> resolution()
    {
        if (image == NULL)
            pre_read();

        return res;
    }


    void read(double *const rho_buffer, double *const vel_buffer)
    {
        if (image == NULL)
            pre_read();
        
        double *rho_raw = rho->GetPointer(0);
        double *vel_raw = vel->GetPointer(0);

        for (std::size_t i = 0; i < res[0] * res[1] * (order == 3 ? res[2] : 1); i++)
            rho_buffer[i] = rho_raw[i];
        for (std::size_t i = 0; i < res[0] * res[1] * (order == 3 ? res[2] : 1) * 3; i++)
            vel_buffer[i] = vel_raw[i];

        reader->Delete();
        reader = NULL;

        pre_write();
        set_buffers(rho_buffer, vel_buffer);
        write(0);
    }


    void set_buffer(double *const rho_buffer)
    {
        rho->SetVoidArray(rho_buffer, res[0] * res[1] * (order == 3 ? res[2] : 1), 1);
    }


    void set_buffers(double *const rho_buffer, double *const vel_buffer)
    {
        rho->SetVoidArray(rho_buffer, res[0] * res[1] * (order == 3 ? res[2] : 1)    , 1);
        vel->SetVoidArray(vel_buffer, res[0] * res[1] * (order == 3 ? res[2] : 1) * 3, 1);
    }


    void no_vel()
    {
        vel->Delete();
        vel = NULL;
    }


    void write(const double time)
    {
        if (time_count % time_steps_per_file == 0)
        {
            if (time_count != 0)
                writer->Stop();

            writer->SetFileName(fmt::format("{}_{:04d}.vti", name_no_suffix, time_count / time_steps_per_file).c_str());
            writer->Start();
        }

        rho->Modified();
        if (vel != NULL)
            vel->Modified();
        writer->WriteNextTime(time);
        time_count++;
    }
};

}
