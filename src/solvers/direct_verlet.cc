#include "direct_verlet.hh"
#include "../storage/storage.hh"

int main()
{
    // number of particles
    size_t n {4};

    // time step
    double dt {0.00001};

    // # time steps, including ic
    size_t n_t {10000001};

    // store data every so many steps
    size_t n_s {1000};



    // set up h5 storage
    Storage storage {{"direct_verlet.h5", H5F_ACC_TRUNC}};

    double *state = new double[4*n];

    // initial condition
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

    storage.new_dataset(ic, n, 0);

    ic_init(ic, state, n, dt);
    delete[] ic;

    for (size_t t {1}; t < n_t; t++)
    {   if (t % n_s == 0)
            storage.new_dataset(state, n, t/n_s);
        update(state, n, dt);
    }

    delete[] state;
}
