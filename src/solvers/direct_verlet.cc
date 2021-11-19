#include "direct_verlet.hh"
#include "../storage/storage.hh"
#include <iostream>

int main()
{
    /* number of time steps  */ std::size_t n_t = 1000000;
    /* steps between process */ std::size_t n_s = 1000;
    /* time step             */ double dt   = 1e-5;
    /* softening eps ^2      */ double eps2 = 3e-5;
    /* h5 file with init con */ Storage storage {{"0.h5", H5F_ACC_RDWR}};
    /* number of objs        */ std::size_t n = storage.read_dataset_size("ic");
    /* integrator memory     */ auto *ic    = new double[6*n];
                                auto *state = new double[6*n];

    /* read and store ic     */ storage.read_dataset(ic, "ic");
                                storage.new_dataset(ic, n, "0");



    //  process state # s
    auto process = [&](std::size_t s)
    {   storage.new_dataset(state, n, std::to_string(s / n_s));
        std::cout << s << '/' << n_t << '\n';
    };



    using namespace Direct_Verlet;

    //  initialization step
    forward_init(ic, state, n, dt, eps2);
    if (n_s == 1)
        process(1);

    delete[] ic;


    //  main loop
    for (std::size_t s = 2; s <= n_t; s++)
    {   forward(state, n, dt, eps2);
        if (s % n_s == 0)
            process(s);
    }

    delete[] state;
}
