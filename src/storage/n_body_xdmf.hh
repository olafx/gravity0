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
    vtkNew<vtkDoubleArray> points_data;
    vtkNew<vtkPoints> points;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkXdmf3Writer> writer;

    N_Body_xdmf(const std::string& name)
    {   points_data->SetNumberOfComponents(3);
        points->SetDataType(VTK_DOUBLE);
        points->SetData(points_data);
        grid->SetPoints(points);
        writer->SetFileName(name.c_str());
        writer->SetLightDataLimit(0);
        writer->SetInputData(grid);
    }

    void set(double *const data, const std::size_t n) const
    {   points_data->SetVoidArray(data, 3 * n, 1);
    }

    void write()
    {   writer->Write();
    }
};

}
