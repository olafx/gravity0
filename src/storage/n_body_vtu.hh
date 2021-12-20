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
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkInformation.h>
#include <vtkCompositeDataPipeline.h>
#include <vtkDoubleArray.h>
#include <string>
#include <utility>
#include <fmt/format.h>

namespace Storage
{

//  time dependent 3D n-body storage via VTK unstructured grid format (multiple files, time stored)

template <std::size_t time_steps_per_file = 16384>
struct N_Body_vtu
{
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkPoints> points;
    vtkNew<vtkDoubleArray> points_data;
    std::size_t time_count;
    std::string name_no_suffix;

    N_Body_vtu(std::string name, double *const data, const int n)
        : time_count {0}, name_no_suffix {std::move(name)}
    {
        points_data->SetVoidArray(data, 3 * n, 1);
        points_data->SetNumberOfComponents(3);

        points->SetDataType(VTK_DOUBLE);
        points->SetNumberOfPoints(n);
        points->SetData(points_data);

        grid->SetPoints(points);

        writer->SetInputData(grid);
        writer->SetNumberOfTimeSteps(time_steps_per_file);
    }

    void write(const double time)
    {
        if (time_count % time_steps_per_file == 0)
        {   if (time_count != 0)
                writer->Stop();
            writer->SetFileName(fmt::format("{}_{:04d}.vtu", name_no_suffix, time_count / time_steps_per_file).c_str());
            writer->Start();
        }
        points->Modified();
        writer->WriteNextTime(time);
        time_count++;
    }

    ~N_Body_vtu()
    {   writer->Stop();
    }
};

}
