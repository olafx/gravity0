#pragma once
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXdmf3Writer.h>
#include <vtkDoubleArray.h>

struct Storage2
{
    vtkSmartPointer<vtkDoubleArray>      points_data = decltype(points_data)::New();
    vtkSmartPointer<vtkPoints>           points      = decltype(points)     ::New();
    vtkSmartPointer<vtkUnstructuredGrid> grid        = decltype(grid)       ::New();
    vtkSmartPointer<vtkXdmf3Writer>      writer      = decltype(writer)     ::New();

    Storage2(const std::string& name)
    {
        points_data->SetNumberOfComponents(3);

        points->SetDataType(VTK_DOUBLE);
        points->SetData(points_data);

        grid->SetPoints(points);

        writer->SetFileName(name.c_str());
        writer->SetLightDataLimit(0);
        writer->SetInputData(grid);
    }

    void set(double *const data, const std::size_t n)
    {
        points_data->SetVoidArray(data, 3 * n, 1);
    }

    void write()
    {
        writer->Write();
    }
};
