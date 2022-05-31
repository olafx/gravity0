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

#include <fftw3.h>
#include <cmath>
#include <array>

namespace Field_Periodic
{
    // rho is in contrast format rn in data[3:4], turn it into absolute rho in data[3:4]
    void init_rho(fftw_complex *data, std::array<int, 2> res, size_t size, double rho_mean) noexcept
    {
        fftw_complex *rho = data + 3 * size;
        for (size_t i = 0; i < size; i++)
            rho[i][0] = rho_mean * (1 + rho[i][0]);
    }

    // turn rho in data[3:4] and vel in data[4:6] in their FT and store in data[0:1] and data[1:3]
    void init(fftw_complex *data, std::array<int, 2> res, size_t size) noexcept
    {
        fftw_plan plan_1 = fftw_plan_dft_2d(res[0], res[1], data + 3 * size, data           , FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan plan_2 = fftw_plan_dft_2d(res[0], res[1], data + 4 * size, data + 1 * size, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan plan_3 = fftw_plan_dft_2d(res[0], res[1], data + 5 * size, data + 2 * size, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan_1);
        fftw_execute(plan_2);
        fftw_execute(plan_3);
        fftw_destroy_plan(plan_1);
        fftw_destroy_plan(plan_2);
        fftw_destroy_plan(plan_3);
    }

    // at this point rho and vel in freq space are set in data[0:1] and data[1:3],
    // and we calculate their inverse FTs and store them in data[3:4] and data[4:6].
    // do complex to real so that it stores right.
    // also the density must become a density perturbation again
    void prep_for_store(fftw_complex *data, std::array<int, 2> res, size_t size, double rho_mean) noexcept
    {
        fftw_plan plan_1 = fftw_plan_dft_2d(res[0], res[1], data,            data + 3 * size, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan plan_2 = fftw_plan_dft_2d(res[0], res[1], data + 1 * size, data + 4 * size, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan plan_3 = fftw_plan_dft_2d(res[0], res[1], data + 2 * size, data + 5 * size, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_1);
        fftw_execute(plan_2);
        fftw_execute(plan_3);
        fftw_destroy_plan(plan_1);
        fftw_destroy_plan(plan_2);
        fftw_destroy_plan(plan_3);

        fftw_complex *rho = data + 3 * size;
        fftw_complex *vel = data + 4 * size;
        double *rho_ = static_cast<double *>(*rho);
        double *vel_ = static_cast<double *>(*vel);
        for (size_t i = 0; i < size; i++)
            rho_[i] = sqrt(rho[i][0] * rho[i][0] + rho[i][1] * rho[i][1]);
        for (size_t i = 0; i < 2 * size; i++)
            vel_[i] = sqrt(vel[i][0] * vel[i][0] + vel[i][1] * vel[i][1]);

        for (size_t i = 0; i < size; i++)
            rho_[i] = rho_[i] / rho_mean - 1;
    }

    void forward(fftw_complex *data, std::array<int, 2> res, size_t size, double G, double dt) noexcept
    {
        // update vel in data[1:3]
        //     multiplying by i means what? (a+bi)i = -b+ai, so 
        //     for adding to our real part, look at negatige imaginary part, and
        //     for adding to our imaginary part, look at real part
        fftw_complex *rho = data;
        fftw_complex *vel_x = data + 1 * size;
        fftw_complex *vel_y = data + 2 * size;
        for (size_t i = 0; i < res[0]; i++)
            for (size_t j = 0; j < res[1]; j++)
            {   // we interpret the space step to be 1. here we calculate what
                // the frequency should be under that assumption, as returned by FFT.
                double k_x = j <= .5 * res[0] ? static_cast<double>(j) / res[0] : static_cast<double>(j - res[0]) / res[0];
                double k_y = i <= .5 * res[1] ? static_cast<double>(i) / res[1] : static_cast<double>(i - res[1]) / res[1];
                double k = sqrt(k_x * k_x + k_y * k_y);
                if (k == 0) // dumb temp fix for when k = 0, idrk what do about this
                    k = 1;
                size_t p = i*res[0]+j;
                vel_x[p][0] += dt * 2 * (-rho[p][1]) * G * k_x / (k * k);
                vel_x[p][1] += dt * 2 * ( rho[p][0]) * G * k_x / (k * k);
                vel_y[p][0] += dt * 2 * (-rho[p][1]) * G * k_y / (k * k);
                vel_y[p][1] += dt * 2 * ( rho[p][0]) * G * k_y / (k * k);
            }

        // inverse FT of rho in freq in data[0:1], into data[3:4]
        // 2 inverse FTs of vel in freq in data[1:3], goes into data[4:6]
        fftw_plan plan_1 = fftw_plan_dft_2d(res[0], res[1], data,            data + 3 * size, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan plan_2 = fftw_plan_dft_2d(res[0], res[1], data + 1 * size, data + 4 * size, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan plan_3 = fftw_plan_dft_2d(res[0], res[1], data + 2 * size, data + 5 * size, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_1);
        fftw_execute(plan_2);
        fftw_execute(plan_3);
        fftw_destroy_plan(plan_1);
        fftw_destroy_plan(plan_2);
        fftw_destroy_plan(plan_3);
        // evaluate product of above 2 and put in data[4:6]
        //     an FFT returns frequencies multiplied by size. this isn't a problem anywhere because
        //     we stay in freq space anyway and don't use freq data directly, but since we calculate
        //     product of two FFTs we'll have multiplied twice by size, so need to divide this product
        //     by size after.
        rho = data + 3 * size;
        vel_x = data + 4 * size;
        vel_y = data + 5 * size;
        for (size_t i = 0; i < size; i++)
        {   // (a+bi)(c+di) = (ac-bd)+(ad+bc)i
            vel_x[i][0] = rho[i][0] * vel_x[i][0] - rho[i][1] * vel_x[i][1];
            vel_y[i][0] = rho[i][0] * vel_y[i][0] - rho[i][1] * vel_y[i][1];
            vel_x[i][0] /= size;
            vel_y[i][0] /= size;
        }
        // 2 FTs of that product in data[4:6] in place
        fftw_plan plan_4 = fftw_plan_dft_2d(res[0], res[1], data + 4 * size, data + 4 * size, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan plan_5 = fftw_plan_dft_2d(res[0], res[1], data + 5 * size, data + 5 * size, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan_4);
        fftw_execute(plan_5);
        fftw_destroy_plan(plan_4);
        fftw_destroy_plan(plan_5);
        // update rho in data[0:1]
        rho = data;
        vel_x = data + 1 * size;
        vel_y = data + 2 * size;
        for (size_t i = 0; i < res[0]; i++)
            for (size_t j = 0; j < res[1]; j++)
            {   double k_x = j <= .5 * res[0] ? static_cast<double>(j) / res[0] : static_cast<double>(j - res[0]) / res[0];
                double k_y = i <= .5 * res[1] ? static_cast<double>(i) / res[1] : static_cast<double>(i - res[1]) / res[1];
                size_t p = i*res[0]+j;
                rho[p][0] += -dt * 2 * M_PI * (-vel_x[p][1] * k_x - vel_y[p][1] * k_y);
                rho[p][1] += -dt * 2 * M_PI * ( vel_x[p][0] * k_x + vel_y[p][0] * k_y);
            }
    }
}
