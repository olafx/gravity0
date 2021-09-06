#include <cstddef>
#include <cmath>
#include <iostream>
#include "storage.hh"

int main()
{
    // number of particles
    size_t n {4};

    // time step
    double dt {0.00001};

    // time steps, including ic
    size_t n_t {40000001};

    // store data every how manieth step
    size_t n_s {1000};

    // name of current dataset ()
    size_t i_s {0};

    // set up h5 storage
    Storage storage;

    // initial condition
    // storage format
    //   particle 1 pos x, particle 1 pos y,
    //   ...,
    //   particle n pos x, particle n pos y,
    //   particle 1 vel x, particle 1 vel y,
    //   ...,
    //   particle n vel x, particle n vel y
    double *ic = new double[4*n]
    {   0.98, 0,
        0, 1,
        -1, 0,
        0, -1,

        0, 1,
        -1, 0,
        0, -1,
        1, 0
    };

    // store initial condition
    storage.add(ic, n, i_s++);

    // state
    // storage format
    //   particle 1 pos 2 x, particle 1 pos 2 y,
    //   ...,
    //   particle n pos 2 x, particle n pos 2 y
    //   particle 1 pos 1 x, particle 1 pos 1 y,
    //   ...,
    //   particle n pos 1 x, particle 1 pos 1 y,
    // pos 2 first since first requires less index math and pos 2 needs more computation
    double *s = new double[4*n];

    // state pos 1 is initially ic
    for (size_t i {0}; i < n; i++)
    {   // particle i pos 1 x
        s[2*(i+n)] = ic[2*i];
        // particle i pos 1 y
        s[2*(i+n)+1] = ic[2*i+1];
    }

    // state pos 2 calculated via 2nd order Taylor series
    // start by setting 0
    for (size_t i {0}; i < n; i++)
    {   // particle i pos 1 x
        s[2*i] = 0;
        // particle i pos 1 y
        s[2*i+1] = 0;
    }
    // precalculate recurring number
    {
    double a {0.5 * dt};
    // iterate over pairs
    for (size_t i {0}; i < n; i++)
    {   for (size_t j {0}; j < n; j++)
        {   // 2nd term Taylor
            if (i != j)
            {   // pos 1 distance
                double b {s[2*(i+n)]-s[2*(j+n)]};
                double c {s[2*(i+n)+1]-s[2*(j+n)+1]};
                // pos 1 common acceleration factor
                double d {1/sqrt(b*b+c*c)};
                d = d*d*d;
                // add acceleration to particle i pos 2, fix factor later
                s[2*i] -= d*b;
                s[2*i+1] -= d*c;
            }
        }
        // partially fix 2nd Taylor factor, add velocity from ic, fix both 1st and 2nd Taylor factors
        s[2*i] = (s[2*i]*a+ic[2*(i+n)])*dt;
        s[2*i+1] = (s[2*i+1]*a+ic[2*(i+n)+1])*dt;
        // 0th term Taylor, add pos 1 to pos 2
        s[2*i] += s[2*(i+n)];
        s[2*i+1] += s[2*(i+n)+1];
    }
    }

    // done with first integration step

    // Verlet integration
    {
    // precalculate recurring number
    double a {dt * dt};
    for (size_t t {1}; t < n_t; t++)
    {
        // iterate over pairs
        for (size_t i {0}; i < n; i++)
        {   // acceleration counters initially 0
            double b {0};
            double c {0};
            for (size_t j {0}; j < n; j++)
            {   if (i != j)
                {   // pos 1 distance
                    double d {s[2*i]-s[2*j]};
                    double e {s[2*i+1]-s[2*j+1]};
                    // pos 1 common acceleration factor
                    double f {1/sqrt(d*d+e*e+0.0001)};
                    f = f*f*f;
                    // add to acceleration counters, fix factor later
                    b -= f*d;
                    c -= f*e;
                }
            }
            // set pos 1 to pos 2 but save old pos 1
            double d {s[2*(i+n)]};
            double e {s[2*(i+n)+1]};
            s[2*(i+n)] = s[2*i];
            s[2*(i+n)+1] = s[2*i+1];
            // new pos 2, b and c factor fixed without writing to b and c
            s[2*i] += s[2*i]-d+b*a;
            s[2*i+1] += s[2*i+1]-e+c*a;
        }

        // store state
        if (t % n_s == 0)
            storage.add(s, n, i_s++);
    }
    }
}
