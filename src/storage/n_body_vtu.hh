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
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkAbstractArray.h>
#include <string>
#include <cstring>
#include <utility>

namespace Storage
{

//  time dependent 3D n-body storage via VTK unstructured grid format (multiple files, time stored)
//  TODO
//    - Figure out how to read TimeValues (written by WriteNextTime) rather than TimeStep, which is an integer count.
//      Can then implement multi time step and custom time initial condition reading.

template <size_t time_steps_per_file = 256>
class N_Body_vtu
{
    vtkXMLUnstructuredGridWriter *writer = NULL;
    vtkXMLUnstructuredGridReader *reader = NULL;
    vtkUnstructuredGrid *grid = NULL;
    vtkPoints *points = NULL;
    vtkDoubleArray *pos = NULL;
    vtkDoubleArray *vel = NULL;
    vtkIdType n_objs;
    size_t time_count = 0;
    std::string name_no_suffix;

    void pre_read()
    {   grid = vtkUnstructuredGrid::SafeDownCast(reader->GetOutputDataObject(0));
        pos = vtkDoubleArray::SafeDownCast(grid->GetPoints()->GetData());
        vel = vtkDoubleArray::SafeDownCast(grid->GetPointData()->GetAbstractArray("Velocity"));
        n_objs = grid->GetNumberOfPoints();
        // reader->GetTimeStep();
        // reader->SetTimeStep();
        // reader->GetNumberOfTimeSteps();
    }

    void pre_write()
    {   writer = vtkXMLUnstructuredGridWriter::New();
        grid = vtkUnstructuredGrid::New();
        points = vtkPoints::New();
        pos = vtkDoubleArray::New();
        vel = vtkDoubleArray::New();
        pos->SetNumberOfComponents(3);
        vel->SetName("Velocity");
        vel->SetNumberOfComponents(3);
        points->SetDataType(VTK_DOUBLE);
        points->SetNumberOfPoints(n_objs);
        points->SetData(pos);
        grid->SetPoints(points);
        grid->GetPointData()->AddArray(vel);
        writer->SetInputData(grid);
        writer->SetNumberOfTimeSteps(time_steps_per_file);
    }

public:

    explicit N_Body_vtu(std::string name)
        : name_no_suffix {std::move(name)}
    {   reader = vtkXMLUnstructuredGridReader::New();
        reader->SetFileName((name_no_suffix + ".vtu").c_str());
        reader->Update();
    }

    ~N_Body_vtu()
    {   if (writer != NULL)
        {   writer->Stop();
            pos->Delete();
            if (vel != NULL)
                vel->Delete();
            points->Delete();
            grid->Delete();
            writer->Delete();
        }
        if (reader != NULL)
        {   reader->Delete();
        }
    }

    vtkIdType n_objects()
    {   if (grid == NULL)
            pre_read();
        return n_objs;
    }

    void read(double *const pos_buffer, double *const vel_buffer)
    {   if (grid == NULL)
            pre_read();
        double *pos_raw = pos->GetPointer(0);
        double *vel_raw = vel->GetPointer(0);
        memcpy(pos_buffer, pos_raw, 3 * n_objs * sizeof(double));
        memcpy(vel_buffer, vel_raw, 3 * n_objs * sizeof(double));
        reader->Delete();
        reader = NULL;
        pre_write();
        set_buffers(pos_buffer, vel_buffer);
        write(0);
    }

    void set_buffer(double *const pos_buffer)
    {   pos->SetVoidArray(pos_buffer, 3 * n_objs, 1);
    }

    void set_buffers(double *const pos_buffer, double *const vel_buffer)
    {   pos->SetVoidArray(pos_buffer, 3 * n_objs, 1);
        vel->SetVoidArray(vel_buffer, 3 * n_objs, 1);
    }

    void no_vel()
    {   vel->Delete();
        vel = NULL;
    }

    void write(const double time)
    {   if (time_count % time_steps_per_file == 0)
        {   if (time_count != 0)
                writer->Stop();
            writer->SetFileName((name_no_suffix + std::to_string(time_count / time_steps_per_file)).c_str());
            writer->Start();
        }
        pos->Modified();
        if (vel != NULL)
            vel->Modified();
        writer->WriteNextTime(time);
        time_count++;
    }
};

}
