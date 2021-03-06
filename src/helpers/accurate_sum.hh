/*
    MIT License

    Copyright (c) 2021 Olaf Willocx

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#pragma once
#include <cstddef>
#include <numeric>

namespace Accurate_Sum
{
    double kahan_babuska(const double *const data, const size_t n)
    {
        double sum = data[0] + data[1];
        double com = sum - data[1];
        for (size_t i = 2; i < n; i++)
        {   const double a = data[i] - com;
            const double b = sum + a;
            com = b - sum - a;
            sum = b;
        }
        return sum;
    }

    template <size_t n_direct = 128>
    double pairwise(const double *const data, const size_t n)
    {
        if (n <= n_direct)
        {   double sum = 0;
            for (size_t i = 0; i < n; i++)
                sum += data[i];
            return sum;
        }
        else
        {   return pairwise<n_direct>(data, n / 2) + pairwise<n_direct>(data + n / 2, n / 2);
        }
    }
}
