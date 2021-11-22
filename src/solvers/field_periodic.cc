#include "field_periodic.hh"
#include "../storage/f_vti.hh"

int main()
{
    //  temp test
    std::array<int, 2> dims {5, 5};
    auto *data = new double[dims[0] * dims[1]];
    Storage::f_vti storage {"field_periodic", dims, data};
    for (std::size_t i = 0; i < 10; i++)
    {   for (std::size_t j = 0; j < dims[0] * dims[1]; j++)
            data[j] = i;
        storage.write();
    }
    delete[] data;
}
