#include "field_periodic.hh"
#include "../storage/field_vti.hh"
#include "../storage/n_body_vtu.hh"
#include <random>

int main()
{
    std::random_device rd;
    std::mt19937_64 mt {rd()};
    std::uniform_real_distribution<> dist {0, 1};

    std::array<int, 2> dims {5, 5};
    auto *data_1 = new double[dims[0] * dims[1]];
    Storage::Field_vti<4, 2> storage_1 {"field", data_1, dims};
    for (std::size_t i = 0; i < 10; i++)
    {   for (std::size_t j = 0; j < dims[0] * dims[1]; j++)
            data_1[j] = sqrt(.2 * i) + j;
        storage_1.write(.1 * i);
    }

    int size = 256;
    auto *data_2 = new double[3 * size];
    Storage::N_Body_vtu<512> storage_2 {"n_body", data_2, size};
    for (std::size_t i = 0; i < 512 * 8; i++)
    {   for (std::size_t j = 0; j < 3 * size; j++)
            data_2[j] = dist(mt);
        storage_2.write(.1 * i);
    }

    std::cout << data_2 << '\n';

    delete[] data_1;
    delete[] data_2;
}
