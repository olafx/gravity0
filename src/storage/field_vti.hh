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
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkXMLImageDataWriter.h>
#include <string>
#include <array>
#include <type_traits>
#include <utility>
#include <fmt/format.h>

namespace Storage
{

//  time dependent 2D or 3D uniform grid storage via VTK image format (multiple files, time stored)

template <std::size_t time_steps_per_file = 16384, std::size_t order = 2,
          typename = std::enable_if_t<order == 2 || order == 3>>
struct Field_vti
{
    vtkNew<vtkXMLImageDataWriter> writer;
    vtkNew<vtkImageImport> importer;
    std::size_t time_count = 0;
    std::string name_no_suffix;

    Field_vti(std::string name, double *const data, const std::array<int, order>& dims)
        : name_no_suffix {std::move(name)}
    {
        importer->SetDataSpacing(1, 1, 1);
        importer->SetDataOrigin(0, 0, 0);
        importer->SetWholeExtent(0, dims[0] - 1, 0, dims[1] - 1, 0, order == 2 ? 0 : dims[2] - 1);
        importer->SetDataExtentToWholeExtent();
        importer->SetDataScalarType(VTK_DOUBLE);
        importer->SetNumberOfScalarComponents(1);
        importer->SetImportVoidPointer(data);

        writer->SetInputConnection(importer->GetOutputPort());
        writer->SetNumberOfTimeSteps(time_steps_per_file);
    }

    void write(const double time)
    {
        if (time_count % time_steps_per_file == 0)
        {   if (time_count != 0)
                writer->Stop();
            writer->SetFileName(fmt::format("{}_{:04d}.vti", name_no_suffix, time_count / time_steps_per_file).c_str());
            writer->Start();
        }
        importer->Modified();
        writer->WriteNextTime(time);
        time_count++;
    }

    ~Field_vti()
    {
        writer->Stop();
    }
};

}
