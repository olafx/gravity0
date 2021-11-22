#include "field_periodic.hh"
#include "../storage/f_vti.hh"
#include "../storage/n_vtu.hh"

int main()
{
    std::array<int, 2> dims {5, 5};
    auto *data_1 = new double[dims[0] * dims[1]];
    Storage::Field_vti<4, 2> storage_1 {"field", data_1, dims};
    for (std::size_t i = 0; i < 10; i++)
    {   for (std::size_t j = 0; j < dims[0] * dims[1]; j++)
            data_1[j] = .1 * i;
        storage_1.write(.1 * i);
    }

    int size = 5;
    auto *data_2 = new double[3 * 5];
    Storage::N_Body_vtu<4> storage_2 {"n_body", data_2, size};
    for (std::size_t i = 0; i < 10; i++)
    {   for (std::size_t j = 0; j < 3 * 5; j++)
            data_2[j] = .1 * i;
        storage_2.write(.1 * i);
    }

    delete[] data_1;
    delete[] data_2;
}
