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
#include <array>

#include <iostream>

int main()
{
    // {
    //     Storage::Field_vti<5> storage {"0"};
    //     std::array<int, 2> res = storage.resolution();
    //     size_t size = res[0] * res[1];
    //     auto *data = new double[(1 + 3) * size];
    //     double *rho = data;
    //     double *vel = data + size;
    //     storage.read(rho, vel);
    //     for (size_t i = 1; i < 13; i++)
    //     {   for (size_t j = 0; j < size; j++)
    //             rho[j] = i;
    //         for (size_t j = 0; j < size * 3; j++)
    //             vel[j] = .1 * i;
    //         storage.write(.01 * i);
    //     }
    // }

    double rho_mean = 1e-4;
    double dt = 1e-3;
    double G = 1;
    size_t N = 100;
    size_t n = 10;

    Storage::Field_vti<16, 2> storage {"0"};
    std::array<int, 2> res = storage.resolution();
    size_t size = res[0] * res[1];

    auto *data = new fftw_complex[size * (1 + 2 + 1 + 2)];
    for (size_t i = 0; i < size * (1 + 2 + 1 + 2); i++)
    {   data[i][0] = 0;
        data[i][1] = 0;
    }

    storage.read(static_cast<double *>(*(data + 3 * size)),
                 static_cast<double *>(*(data + 4 * size))); // this also sets the buffer to where it will write from
    {   fftw_complex *rho = data + 3 * size;
        fftw_complex *vel = data + 4 * size;
        double *rho_ = static_cast<double *>(*rho);
        double *vel_ = static_cast<double *>(*vel);
        for (size_t i = size - 1; i != 0; i--)
            rho[i][0] = rho_[i];
        rho[0][0] = rho_[0];
        for (size_t i = 2 * size - 1; i != 0; i--)
            vel[i][0] = vel_[i];
        vel[0][0] = vel_[0];
    }

    Field_Periodic::init_rho(data, res, size, rho_mean);
    Field_Periodic::init(data, res, size);
    // since initial velocity was set to 0, its FT might be bad, although it seems fine in numpy's FFT

    for (size_t s = 1; s <= N; s++)
    {   printf("{%zu}/{%zu}\n", s, N);
        Field_Periodic::forward(data, res, size, G, dt);
        if (s % n == 0)
        {   Field_Periodic::prep_for_store(data, res, size, rho_mean);
            storage.write(s * dt);
        }
    }

    delete[] data;
}
