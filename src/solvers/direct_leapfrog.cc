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

#include "direct_leapfrog.hh"
#include "n_body_h5.hh"
#include "n_body_vtu.hh"
#include <fmt/core.h>

int main()
{
    /* number of time steps  */ constexpr std::size_t n_t = 1000000;
    /* steps between process */ constexpr std::size_t n_s = 1000;
    /* time step             */ constexpr double dt   = 1e-6;
    /* softening eps^2       */ constexpr double eps2 = 3e-5;
    /* h5 file with init con */ Storage::N_Body_vtu storage {"king"};
    /* number of objs        */ std::size_t n = storage.n_objects();
    /* integrator memory     */ auto *state = new double[6 * n];
    /* store vel also?       */ constexpr bool store_velocities = true;


    //  read ic, pos and vel
    storage.read(state, state + 3 * n);


    //  vel buffer was set by 'read', so need explicit 'no_vel' to not update
    if constexpr (!store_velocities)
        storage.no_vel();


    //  process state function
    auto process = [&](std::size_t s)
    {   storage.write(s * dt);
        fmt::print("{}/{}\n", s, n_t);
    };


    //  main loop
    for (std::size_t s = 1; s <= n_t; s++)
    {   Direct_Leapfrog::forward(state, n, dt, eps2);
        if (s % n_s == 0)
            process(s);
    }

    delete[] state;
}
