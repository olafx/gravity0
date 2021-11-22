#pragma once
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkInformation.h>
#include <vtkCompositeDataPipeline.h>
#include <vtkDoubleArray.h>
#include <string>
#include <fmt/format.h>

namespace Storage
{

//  time dependent n-body storage via VTK unstructured grid format (multiple files, time not stored)
//  TODO
//    - points needs to use data; need to find how to set void ptr

template <std::size_t time_steps_per_file = 16384>
struct n_vtu
{
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkPoints> points;
    std::size_t time_count;
    std::string name_no_suffix;

    //  VTK uses integers for counts and doesn't mark data pointer for writing const
    n_vtu(const std::string& name, double *const data, const int n)
        : time_count {0}, name_no_suffix {name}
    {
        points->SetDataType(VTK_DOUBLE);
        points->SetNumberOfPoints(n);
        //  TODO temp nonsense
        for (std::size_t i = 0; i < n; i++)
            points->SetPoint(i, data[3*i], data[3*i+1], data[3*i+2]);

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
        time_count++;
        points->Modified();
        writer->WriteNextTime(time);
    }

    ~n_vtu()
    {   writer->Stop();
    }
};

}
