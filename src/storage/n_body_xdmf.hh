#pragma once
#include <vtkNew.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXdmf3Writer.h>

namespace Storage
{

//  time dependent 3D n-body storage via VTK unstructured grid format
//  TODO
//    - time dependence
//    - constructor needs to be more like the others; no set function

struct N_Body_xdmf
{
    vtkNew<vtkDoubleArray>      points_data;
    vtkNew<vtkPoints>           points;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkXdmf3Writer>      writer;

    N_Body_xdmf(const std::string& name)
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

}
