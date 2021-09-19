#include "ic_gen.hh"
#include "../storage/storage.hh"

int main()
{   using namespace std;

    size_t n {4};
    
    Storage storage {{"0.h5", H5F_ACC_EXCL}};

    mt19937_64 generator;
    normal_distribution<double> distribution {0, 10};

    double *ic = new double[6*n];

    for (size_t i {}; i < 6*n; i++)
        ic[i] = distribution(generator);

    storage.new_dataset(ic, 2*n, "ic");

    delete[] ic;
}
