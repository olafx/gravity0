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

int main()
{
    /* number of time steps  */ constexpr size_t n_t = 300000;
    /* steps between process */ constexpr size_t n_s = 1000;
    /* time step             */ constexpr double dt   = 1e-5;
    /* softening eps ^2      */ constexpr double eps2 = 3e-5;
    /* h5 file with ic       */ Storage::N_Body_h5 storage {"king"};
                                size_t n = storage.n_objects();
    /* integrator memory     */ auto *ic    = new double[6 * n];
                                auto *state = new double[6 * n];
    /* store vel also?       */ constexpr bool store_velocities = false;

    //  read ic, pos and vel
    storage.read(ic, ic + 3 * n);

    //  process state # s
    auto process = [&](size_t s)
    {   if constexpr (store_velocities)
            storage.write(state, state + 3 * n, s * dt);
        else
            storage.write(state, s * dt);
        printf("{%zu}/{%zu}\n", s, n_t);
    };

    //  initialization step
    Direct_Verlet::forward_init(ic, state, n, dt, eps2);
    if (n_s == 1)
        process(1);

    delete[] ic;

    //  main loop
    for (size_t s = 2; s <= n_t; s++)
    {   Direct_Verlet::forward(state, n, dt, eps2);
        if (s % n_s == 0)
            process(s);
    }

    delete[] state;
}
