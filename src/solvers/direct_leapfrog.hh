#pragma once
#include <cmath>

/*  Direct Leapfrog

    direct n-body solver
    constant time step
    3D
    1 position
    1 velocity
    symplectic
    softening
    2nd order in space
    1st order in time
    no initializaion
*/

namespace Direct_Leapfrog
{

void forward(double *const state, const std::size_t n,
             const double dt, const double eps2)
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
        state[3*(i+n)  ] += a1 * dt;
        state[3*(i+n)+1] += a2 * dt;
        state[3*(i+n)+2] += a3 * dt;
        state[3*i  ] += state[3*(i+n)  ] * dt;
        state[3*i+1] += state[3*(i+n)+1] * dt;
        state[3*i+2] += state[3*(i+n)+2] * dt;
    }
}

};
