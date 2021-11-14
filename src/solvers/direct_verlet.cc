#include "direct_verlet.hh"
#include "../storage/storage.hh"
#include <fmt/core.h>

int main()
{
    //  number of time steps, excluding initial condition
    std::size_t n_t = 10000000;

    //  store data every so many steps
    std::size_t n_s = 10000;

    //  open h5 file containing initial condition
    Storage storage {{"0.h5", H5F_ACC_RDWR}};

    //  reading number of objects
    std::size_t n = storage.read_dataset_size("ic");

    //  memory
    double *ic    = new double[3 * 2 * n],
           *state = new double[3 * 2 * n];

    //  setting up solver with time step, softening parameter, number of objects, and memory
    Direct_Verlet integrator {1e-6, 3e-5, n, ic, state};



    //  reading initial condition
    storage.read_dataset(integrator.ic(), "ic");

    //  storing initial condition pos as step 0
    storage.new_dataset(integrator.ic(), integrator.n, "0");



    //  things that happen when we care about state number s
    auto process = [&](std::size_t s)
    {   storage.new_dataset(integrator.state(), integrator.n, std::to_string(s / n_s));
        fmt::print("{}/{}\n", s, n_t);
    };



    //  initializing state to be ready for Verlet integration
    //  this counts as the first integration step
    integrator.forward_init();
    if (n_s == 1)
        process(1);

    //  initial condition no longer needed
    delete[] ic;

    //  time loop
    for (std::size_t s = 2; s <= n_t; s++)
    {   integrator.forward();
        if (s % n_s == 0)
            process(s);
    }

    delete[] state;
}
