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
#include <cmath>

namespace Direct_Verlet
{

/*  Direct Verlet

    direct n-body solver
    constant time step
    3D
    2 positions
    0 velocities
    symplectic
    softening
    2nd order in space
    1st order in time
    Taylor series initialization from 1 position 1 velocity data
*/

void forward_init(const double *const ic, double *const state,
                  const std::size_t n, const double dt, const double eps2) noexcept
{
    for (std::size_t i = 0; i < 3 *n; i++)
        state[i] = 0;

    for (std::size_t i = 0; i < n; i++)
    {   state[3*(i+n)  ] = ic[3*i  ];
        state[3*(i+n)+1] = ic[3*i+1];
        state[3*(i+n)+2] = ic[3*i+2];
    }

    #pragma omp parallel for default(none) shared(n, state, ic, eps2, dt)
    for (std::size_t i = 0; i < n; i++)
    {   for (std::size_t j = 0; j < n; j++)
            if (i != j)
            {   const double b1 = state[3*(i+n)  ] - state[3*(j+n)  ];
                const double b2 = state[3*(i+n)+1] - state[3*(j+n)+1];
                const double b3 = state[3*(i+n)+2] - state[3*(j+n)+2];
                double c = 1 / sqrt(b1 * b1 +
                                    b2 * b2 +
                                    b3 * b3 + eps2);
                c = c * c * c;
                state[3*i  ] -= c * b1;
                state[3*i+1] -= c * b2;
                state[3*i+2] -= c * b3;
            }
        state[3*i  ] = state[3*(i+n)  ] + (.5 * dt * state[3*i  ] + ic[3*(i+n)  ]) * dt;
        state[3*i+1] = state[3*(i+n)+1] + (.5 * dt * state[3*i+1] + ic[3*(i+n)+1]) * dt;
        state[3*i+2] = state[3*(i+n)+2] + (.5 * dt * state[3*i+2] + ic[3*(i+n)+2]) * dt;
    }
}



void forward(double *const state,
             const std::size_t n, const double dt, const double eps2) noexcept
{
    #pragma omp parallel for default(none) shared(n, state, eps2, dt)
    for (std::size_t i = 0; i < n; i++)
    {   double a1 = 0;
        double a2 = 0;
        double a3 = 0;
        for (std::size_t j = 0; j < n; j++)
            if (i != j)
            {   const double b1 = state[3*i  ] - state[3*j  ];
                const double b2 = state[3*i+1] - state[3*j+1];
                const double b3 = state[3*i+2] - state[3*j+2];
                double c = 1 / sqrt(b1 * b1 +
                                    b2 * b2 +
                                    b3 * b3 + eps2);
                c = c * c * c;
                a1 -= c * b1;
                a2 -= c * b2;
                a3 -= c * b3;
            }
        const double d1 = state[3*(i+n)  ];
        const double d2 = state[3*(i+n)+1];
        const double d3 = state[3*(i+n)+2];
        state[3*(i+n)  ] = state[3*i  ];
        state[3*(i+n)+1] = state[3*i+1];
        state[3*(i+n)+2] = state[3*i+2];
        state[3*i  ] += state[3*i  ] - d1 + a1 * (dt * dt);
        state[3*i+1] += state[3*i+1] - d2 + a2 * (dt * dt);
        state[3*i+2] += state[3*i+2] - d3 + a3 * (dt * dt);
    }
}

};
