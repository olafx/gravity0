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

#include "field_periodic.hh"
#include "field_vti.hh"
#include "n_body_vtu.hh"
#include "n_body_h5.hh"
#include <random>

int main()
{
    std::random_device rd;
    std::mt19937_64 mt {rd()};
    std::uniform_real_distribution<> dist {0, 1};


    {   
        std::array<int, 2> dims {5, 5};
        constexpr std::size_t steps_per_file = 4;
        std::size_t n_files = 3;
        std::size_t offset  = 1;

        auto *data = new double[dims[0] * dims[1]];

        Storage::Field_vti<steps_per_file> storage {"field", data, dims};

        for (std::size_t i = 0; i < steps_per_file * n_files - offset; i++)
        {   for (std::size_t j = 0; j < dims[0] * dims[1]; j++)
                data[j] = sqrt(.2 * i) + j;
            storage.write(.1 * i);
        }

        delete[] data;
    }


    {   int size = 3;

        constexpr std::size_t steps_per_file = 5;
        std::size_t n_files = 9;
        std::size_t offset  = 1;

        auto *data = new double[3 * size];

        Storage::N_Body_vtu<steps_per_file> storage {"n_body", data, size};

        for (std::size_t i = 0; i < steps_per_file * n_files - offset; i++)
        {   for (std::size_t j = 0; j < 3 * size; j++)
                data[j] = dist(mt);
            storage.write(.1 * i);
        }

        delete[] data;
    }
}
