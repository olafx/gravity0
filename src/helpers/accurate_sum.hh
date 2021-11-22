#pragma once
#include <cstddef>

namespace Accurate_Sum
{
    double kahan_babuska(const double *const data, const std::size_t n)
    {
        double sum = data[0] + data[1];
        double compensator = (data[0] + data[1]) - data[1];
        for (std::size_t i = 2; i < n; i++)
        {   const double a = data[i] - compensator;
            const double b = sum + a;
            compensator = (b - sum) - a;
            sum = b;
        }
        return sum;
    }

    template <std::size_t n_direct = 128>
    double pairwise(const double *const data, const std::size_t n)
    {
        if (n <= n_direct)
        {   double sum = data[0];
            for (std::size_t i = 1; i < n; i++)
                sum += data[i];
            return sum;
        }
        else
            return sum<n_direct>(data, n / 2) + sum<n_direct>(data + n / 2, n / 2);
    }
}
