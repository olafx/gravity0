#include <H5Cpp.h>

int main()
{   using namespace H5;
    H5File fp {"a.h5", H5F_ACC_TRUNC};
}
