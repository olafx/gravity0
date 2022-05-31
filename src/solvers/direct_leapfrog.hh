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

namespace Direct_Leapfrog
{

void forward(double *const state, const size_t n,
             const double dt, const double eps2) noexcept
{
    #pragma omp parallel for default(none) shared(n, state, eps2, dt)
    for (size_t i = 0; i < n; i++)
    {   double a1 = 0;
        double a2 = 0;
        double a3 = 0;
        for (size_t j = 0; j < n; j++)
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
        state[3*(i+n)  ] += a1 * dt;
        state[3*(i+n)+1] += a2 * dt;
        state[3*(i+n)+2] += a3 * dt;
        state[3*i  ] += state[3*(i+n)  ] * dt;
        state[3*i+1] += state[3*(i+n)+1] * dt;
        state[3*i+2] += state[3*(i+n)+2] * dt;
    }
}

};
