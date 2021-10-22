#pragma once
#include <cstddef>
#include <cmath>

/*
initial condition storage format
    particle 1 pos x, particle 1 pos y, particle 1 pos z,
    ...,
    particle n pos x, particle n pos y, particle n pos z,
    particle 1 vel x, particle 1 vel y, particle 1 vel z,
    ...,
    particle n vel x, particle n vel y, particle n vel z

state storage format
    particle 1 pos 2 x, particle 1 pos 2 y, particle 1 pos 2 z,
    ...,
    particle n pos 2 x, particle n pos 2 y, particle n pos 2 z,
    particle 1 pos 1 x, particle 1 pos 1 y, particle 1 pos 1 z,
    ...,
    particle n pos 1 x, particle 1 pos 1 y, particle 1 pos 1 z
*/

/*
initialize state using ic of size n via 2nd order Taylor series with timestep dt
*/
void init(const double *const ic, double *const state, const size_t n, const double dt)
{
    // state pos 2 prep for init
    for (size_t i = 0; i < 3 * n; i++)
        state[i] = 0;

    // state pos 1 is given by ic
    for (size_t i = 0; i < n; i++)
    {   state[3 * (i + n)]     = ic[3 * i];
        state[3 * (i + n) + 1] = ic[3 * i + 1];
        state[3 * (i + n) + 2] = ic[3 * i + 2];
    }

    // state pos 2 via 2nd order Taylor series

    // PP loop
    for (size_t i = 0; i < n; i++)
    {   for (size_t j = 0; j < n; j++)
        {   // 2nd Taylor term
            if (i != j)
            {   // pos 1 distance
                double b1 = state[3 * (i + n)]     - state[3 * (j + n)];
                double b2 = state[3 * (i + n) + 1] - state[3 * (j + n) + 1];
                double b3 = state[3 * (i + n) + 2] - state[3 * (j + n) + 2];
                // pos 1 common acceleration factor
                double c = 1 / sqrt(b1 * b1 + b2 * b2 + b3 * b3);
                c = c * c * c;
                // add acceleration factor
                state[3 * i]     -= c * b1;
                state[3 * i + 1] -= c * b2;
                state[3 * i + 2] -= c * b3;
            }
        }
        // add velocity from ic, fix 1st and 2nd Taylor terms
        state[3 * i]     = (state[3 * i]     * (0.5 * dt) + ic[3 * (i + n)])     * dt;
        state[3 * i + 1] = (state[3 * i + 1] * (0.5 * dt) + ic[3 * (i + n) + 1]) * dt;
        state[3 * i + 2] = (state[3 * i + 2] * (0.5 * dt) + ic[3 * (i + n) + 2]) * dt;
        // 0th Taylor term
        state[3 * i]     += state[3 * (i + n)];
        state[3 * i + 1] += state[3 * (i + n) + 1];
        state[3 * i + 2] += state[3 * (i + n) + 2];
    }
}

/*
update state of size n via direct Verlet integration with timestep dt
*/
void update(double *const state, const size_t n, const double dt)
{
    // PP loop
    for (size_t i = 0; i < n; i++)
    {   // acceleration counters
        double a1 = 0;
        double a2 = 0;
        double a3 = 0;
        for (size_t j = 0; j < n; j++)
        {   if (i != j)
            {   // pos 1 distance
                double b1 = state[3 * i]     - state[3 * j];
                double b2 = state[3 * i + 1] - state[3 * j + 1];
                double b3 = state[3 * i + 2] - state[3 * j + 2];
                // pos 1 common acceleration factor
                double c = 1 / sqrt(b1 * b1 + b2 * b2 + b3 * b3);
                c = c * c * c;
                // add acceleration factor
                a1 -= c * b1;
                a2 -= c * b2;
                a3 -= c * b3;
            }
        }
        // save pos 1
        double d1 = state[3 * (i + n)];
        double d2 = state[3 * (i + n) + 1];
        double d3 = state[3 * (i + n) + 2];
        // set pos 1 to pos 2
        state[3 * (i + n)]     = state[3 * i];
        state[3 * (i + n) + 1] = state[3 * i + 1];
        state[3 * (i + n) + 2] = state[3 * i + 2];
        // new pos 2 via saved pos 1 and acceleration factor
        state[3 * i]     += state[3 * i]     - d1 + a1 * (dt * dt);
        state[3 * i + 1] += state[3 * i + 1] - d2 + a2 * (dt * dt);
        state[3 * i + 2] += state[3 * i + 2] - d3 + a3 * (dt * dt);
    }
}
