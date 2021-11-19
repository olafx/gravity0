#pragma once
#include <cstddef>

namespace Kahan_Babuska
{
    double sum(const double *const data, const std::size_t n)
    {
        double sum = 0;
        double compensator = 0;
        for (std::size_t i = 0; i < n; i++)
        {   const double a = data[i] - c;
            const double b = sum + a;
            compensator = (b - sum) - a;
            sum = t;
        }
        return sum;
    }
}
