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

    //  VTK uses integers for counts and doesn't mark data pointer for writing const
    N_Body_vtu(const std::string& name, double *const data, const int n)
        : time_count {0}, name_no_suffix {name}
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
        time_count++;
        points->Modified();
        writer->WriteNextTime(time);
    }

    ~N_Body_vtu()
    {   writer->Stop();
    }
};

}
