#pragma once
#include <cstddef>
#include <cmath>

/*
initial condition storage format
    particle 1 pos x, particle 1 pos y,
    ...,
    particle n pos x, particle n pos y,
    particle 1 vel x, particle 1 vel y,
    ...,
    particle n vel x, particle n vel y

state storage format
    particle 1 pos 2 x, particle 1 pos 2 y,
    ...,
    particle n pos 2 x, particle n pos 2 y
    particle 1 pos 1 x, particle 1 pos 1 y,
    ...,
    particle n pos 1 x, particle 1 pos 1 y,
*/

/*
initialize state using ic of size n via 2nd order Taylor series with timestep dt
*/
void ic_init(const double *const ic, double *const state, const size_t n, const double dt)
{
    // state pos 2 prep for init
    for (size_t i {}; i < 2*n; i++)
        state[i] = 0;

    // state pos 1 is given by ic
    for (size_t i {}; i < n; i++)
    {   state[2*(i+n)] = ic[2*i];
        state[2*(i+n)+1] = ic[2*i+1];
    }

    // state pos 2 via 2nd order Taylor series

    // PP loop
    for (size_t i {}; i < n; i++)
    {   for (size_t j {}; j < n; j++)
        {   // 2nd Taylor term
            if (i != j)
            {   // pos 1 distance
                double b1 {state[2*(i+n)]-state[2*(j+n)]};
                double b2 {state[2*(i+n)+1]-state[2*(j+n)+1]};
                // pos 1 common acceleration factor
                double c {1/sqrt(b1*b1+b2*b2)};
                c = c*c*c;
                // add acceleration factor
                state[2*i] -= c*b1;
                state[2*i+1] -= c*b2;
            }
        }
        // add velocity from ic, fix 1st and 2nd Taylor terms
        state[2*i] = (state[2*i]*(0.5*dt)+ic[2*(i+n)])*dt;
        state[2*i+1] = (state[2*i+1]*(0.5*dt)+ic[2*(i+n)+1])*dt;
        // 0th Taylor term
        state[2*i] += state[2*(i+n)];
        state[2*i+1] += state[2*(i+n)+1];
    }
}

/*
update state of size n via direct Verlet integration with timestep dt
*/
void update(double *const state, const size_t n, const double dt)
{
    // PP loop
    for (size_t i {}; i < n; i++)
    {   // acceleration counters
        double a1 {};
        double a2 {};
        for (size_t j {}; j < n; j++)
        {   if (i != j)
            {   // pos 1 distance
                double b1 {state[2*i]-state[2*j]};
                double b2 {state[2*i+1]-state[2*j+1]};
                // pos 1 common acceleration factor
                double c {1/sqrt(b1*b1+b2*b2)};
                c = c*c*c;
                // add acceleration factor
                a1 -= c*b1;
                a2 -= c*b2;
            }
        }
        // save pos 1
        double d1 {state[2*(i+n)]};
        double d2 {state[2*(i+n)+1]};
        // set pos 1 to pos 2
        state[2*(i+n)] = state[2*i];
        state[2*(i+n)+1] = state[2*i+1];
        // new pos 2 via saved pos 1 and acceleration factor
        state[2*i] += state[2*i]-d1+a1*(dt*dt);
        state[2*i+1] += state[2*i+1]-d2+a2*(dt*dt);
    }
}
