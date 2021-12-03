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

#include "direct_verlet.hh"
#include "n_body_h5.hh"
#include <fmt/core.h>

int main()
{
    /* number of time steps  */ std::size_t n_t = 1000000;
    /* steps between process */ std::size_t n_s = 1000;
    /* time step             */ double dt   = 1e-5;
    /* softening eps ^2      */ double eps2 = 3e-5;
    /* h5 file with init con */ Storage::N_Body_h5 storage {"0.h5"};
    /* number of objs        */ std::size_t n = storage.read_size("ic");
    /* integrator memory     */ auto *ic    = new double[6 * n];
                                auto *state = new double[6 * n];

    /* read and store ic     */ storage.read(ic, "ic");
                                storage.write(ic, n, "0");



    //  process state # s
    auto process = [&](std::size_t s)
    {   storage.write(state, n, std::to_string(s / n_s));
        fmt::print("{}/{}\n", s, n_t);
    };



    //  initialization step
    Direct_Verlet::forward_init(ic, state, n, dt, eps2);
    if (n_s == 1)
        process(1);

    delete[] ic;


    //  main loop
    for (std::size_t s = 2; s <= n_t; s++)
    {   Direct_Verlet::forward(state, n, dt, eps2);
        if (s % n_s == 0)
            process(s);
    }

    delete[] state;
}
