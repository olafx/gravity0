#include "direct_leapfrog.hh"
#include "../storage/n_hdf5.hh"
#include <iostream>

int main()
{
    /* number of time steps  */ std::size_t n_t = 1000000;
    /* steps between process */ std::size_t n_s = 1000;
    /* time step             */ double dt   = 1e-5;
    /* softening eps^2       */ double eps2 = 3e-5;
    /* h5 file with init con */ Storage::n_hdf5 storage {{"0.h5", H5F_ACC_RDWR}};
    /* number of objs        */ std::size_t n = storage.read_dataset_size("ic");
    /* integrator memory     */ auto *state = new double[6*n];

    /* read and store ic     */ storage.read_dataset(state, "ic");
                                storage.new_dataset(state, n, "0");



    //  process state # s
    auto process = [&](std::size_t s)
    {   storage.new_dataset(state, n, std::to_string(s / n_s));
        std::cout << s << '/' << n_t << '\n';
    };



    //  main loop
    for (std::size_t s = 1; s <= n_t; s++)
    {   Direct_Leapfrog::forward(state, n, dt, eps2);
        if (s % n_s == 0)
            process(s);
    }

    delete[] state;
}
