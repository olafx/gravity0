#pragma once
#include <vtkNew.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDoubleArray.h>
#include <string>
#include <array>
#include <type_traits>
#include <fmt/format.h>

namespace Storage
{

//  time dependent uniform grid storage via VTK image format (multiple files, time not stored)

template <std::size_t order = 2, typename = std::enable_if_t<order == 2 || order == 3>>
struct f_vti
{
    vtkNew<vtkXMLImageDataWriter> writer;
    vtkNew<vtkImageImport> image_import;
    std::size_t time_count;
    std::string name_no_suffix;

    //  VTK uses integers for counts and doesn't mark data pointer for writing const
    f_vti(const std::string& name, const std::array<int, order>& dims, double *const data)
        : time_count {0}, name_no_suffix {name}
    {
        image_import->SetDataSpacing(1, 1, 1);
        image_import->SetDataOrigin(0, 0, 0);
        image_import->SetWholeExtent(0, dims[0] - 1,
                                     0, dims[1] - 1,
                                     0, order == 2 ? 0 : dims[2] - 1);
        image_import->SetDataExtentToWholeExtent();
        image_import->SetDataScalarType(VTK_DOUBLE);
        image_import->SetNumberOfScalarComponents(1);
        image_import->SetImportVoidPointer(data);

        writer->SetInputConnection(image_import->GetOutputPort());
    }

    void write()
    {
        writer->SetFileName(fmt::format("{}_{:08d}.vti", name_no_suffix, time_count++).c_str());
        writer->Write();
    }
};

}
