#include "direct_verlet.hh"
#include "../storage/storage.hh"

int main()
{
    //  time step
    double dt = 0.00001;

    //  number of time steps, including ic
    size_t n_t = 10000001;

    //  store data every so many steps
    size_t n_s = 1000;



    //  open h5 file with initial condition
    Storage storage {{"0.h5", H5F_ACC_RDWR}};

    //  get number of objects from initial condition
    size_t n = storage.read_dataset_size("ic") / 2;

    //  memory allocation
    double *state = new double[6 * n],
           *ic    = new double[6 * n];

    //  reading initial condition
    storage.read_dataset(ic, "ic");

    //  storing initial condition pos as step 0
    storage.new_dataset(ic, n, "0");



    //  initializing state to be ready for Verlet integration
    init(ic, state, n, dt);
    delete[] ic;

    //  time loop
    for (size_t s = 1; s < n_t; s++)
    {   //  store state
        if (s % n_s == 0)
            storage.new_dataset(state, n, std::to_string(s / n_s));
        update(state, n, dt);
    }

    delete[] state;
}
