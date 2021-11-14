#pragma once
#include <Irulan/Static.h>
#include <Irulan/Dynamic.h>
#include <cmath>

struct Direct_Verlet
{
    Irulan::Dynamic::Array<Irulan::Static::Array<double[3]>[1],
                           Irulan::Allocate<false>,
                           Irulan::EfficientShape<true>> ic, state;
    double dt;
    double eps2;
    size_t n;

    Direct_Verlet(double dt_, double eps2_, size_t n_, double *ic_, double *state_)
        : dt {dt_}, eps2 {eps2_}, n {n_}
    {   ic()    = reinterpret_cast<decltype(ic())>   (ic_);
        state() = reinterpret_cast<decltype(state())>(state_);
    }

    //  initialize state using ic via 2nd order Taylor series
    void forward_init()
    {
        //  state pos 2 prep for init
        for (size_t i = 0; i < n; i++)
        {   state(i)(0) = 0;
            state(i)(1) = 0;
            state(i)(2) = 0;
        }

        //  state pos 1 is given by ic_
        for (size_t i = 0; i < n; i++)
        {   state(i + n)(0) = ic(i)(0);
            state(i + n)(1) = ic(i)(1);
            state(i + n)(2) = ic(i)(2);
        }

        //  state pos 2 via 2nd order Taylor series

        //  PP loop
        for (size_t i = 0; i < n; i++)
        {   for (size_t j = 0; j < n; j++)
            {   //  2nd Taylor term
                if (i != j)
                {   //  pos 1 distance
                    double b1 = state(i + n)(0) - state(j + n)(0);
                    double b2 = state(i + n)(1) - state(j + n)(1);
                    double b3 = state(i + n)(2) - state(j + n)(2);
                    //  pos 1 common acceleration factor
                    double c = 1 / sqrt(b1 * b1 + b2 * b2 + b3 * b3 + eps2);
                    c = c * c * c;
                    //  add acceleration factor
                    state(i)(0) -= c * b1;
                    state(i)(1) -= c * b2;
                    state(i)(2) -= c * b3;
                }
            }
            //  add velocity from ic, fix 1st and 2nd Taylor terms
            state(i)(0) = (state(i)(0) * (0.5 * dt) + ic(i + n)(0)) * dt;
            state(i)(1) = (state(i)(1) * (0.5 * dt) + ic(i + n)(1)) * dt;
            state(i)(2) = (state(i)(2) * (0.5 * dt) + ic(i + n)(2)) * dt;
            //  0th Taylor term
            state(i)(0) += state(i + n)(0);
            state(i)(1) += state(i + n)(1);
            state(i)(2) += state(i + n)(2);
        }
    }



    //  integrate state via direct Verlet integration
    void forward()
    {
        //  PP loop
        for (size_t i = 0; i < n; i++)
        {   //  acceleration counters
            double a1 = 0,
                   a2 = 0,
                   a3 = 0;
            for (size_t j = 0; j < n; j++)
            {   if (i != j)
                {   //  pos 1 distance
                    double b1 = state(i)(0) - state(j)(0),
                           b2 = state(i)(1) - state(j)(1),
                           b3 = state(i)(2) - state(j)(2);
                    //  pos 1 common acceleration factor
                    double c = 1 / sqrt(b1 * b1 + b2 * b2 + b3 * b3 + eps2);
                    c = c * c * c;
                    //  add acceleration factor
                    a1 -= c * b1;
                    a2 -= c * b2;
                    a3 -= c * b3;
                }
            }
            //  save pos 1
            double d1 = state(i + n)(0),
                   d2 = state(i + n)(1),
                   d3 = state(i + n)(2);
            //  set pos 1 to pos 2
            state(i + n)(0) = state(i)(0);
            state(i + n)(1) = state(i)(1);
            state(i + n)(2) = state(i)(2);
            //  new pos 2 via saved pos 1 and acceleration factor
            state(i)(0) += state(i)(0) - d1 + a1 * (dt * dt);
            state(i)(1) += state(i)(1) - d2 + a2 * (dt * dt);
            state(i)(2) += state(i)(2) - d3 + a3 * (dt * dt);
        }
    }
};
