'''
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
'''

'''
King model cluster sampling following King (1966).

The calculations don't follow the paper exactly, but the results (should be/are) identical.
-   King evaluates the differential equation describing density as a function of potential numerically.
    The differential equation indeed has no analytic solution, but the solution can be expressed using the error function, which
    can at least be approximated analytically.
-   King uses a boundary condition in the calculation of the partial differential equation describing density as a function
    of radius. This boundary condition was determined based on experimentel data. I found that that this Gaussian differential
    equation is extremely stable under this boundary condition; any remotely reasonable value gives indistinguishable results.
    This makes me question the accuracy of said boundary condition. For this reason and the stability reason, I just pick 0
    for said boundary condition.
-   I don't bother undimensionalizing the final Gaussian differential equation.

The sampling is done via uniform rejection sampling. Since the probability distribution function is numerically evaluated,
rejection sampling is likely the most accurate and fastest sampling method for this problem (although uniform definitely isn't).
'''

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, arccos, linspace, random, array, transpose
from scipy.special import erf
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import h5py

from n_body_vtu import N_Body_vtu
from n_body_h5 import N_Body_h5



#   Make RNG reproducible by creating a generator with seed 0.
generator = random.default_rng(0)



#   Same variable names as King (1966).
k  = .1
j  = 1
V0 = -1
G  = 1

#   The number of radii to be linearly interpolated between when sampling.
N = 10000

#   The number of stars to sample.
n = 8

#   In principle the previous variables define at what distance V hits 0, but the range of distances may be accidentally chosen
#   too small in the numerical evaluation to reach that distance. The r_max variable is the maximum distance of said range.
#   If the program exits with error 'could not continue since r0 > r_max', increase r_max here.
#   Best not to make it too large since this will increase the distance between steps in the range. Make sure there are enough
#   steps from the origin (which is printed during evaluation).
#   Could in principle use estimates and/or trial and error to automate this and not need this variable.
r_max = 10

#   File name.
filename = '0.h5'

#   The name of the dataset in which the data will be stored in the following format as a row major 2 * n by 3 doubles matrix:
#       particle 1 pos x, particle 1 pos y, particle 1 pos z,
#       ...,
#       particle n pos x, particle n pos y, particle n pos z,
#       particle 1 vel x, particle 1 vel y, particle 1 vel z,
#       ...,
#       particle n vel x, particle n vel y, particle n vel z
dataset_name = 'ic'

#   If the potential reaches 0 too fast, interpolation is inaccurate.
#   This is the minimum number of steps to reach the boundary that is acceptable, otherwise program gives a warning.
THRESHOLD_STEPS_TO_BOUNDARY = 16



file_format = filename.rsplit('.')[-1]
if file_format not in ['h5', 'vtu']:
    raise Exception('file format not supported')

#   The density as a function of potential is the result of a non analytic differential equation.
#   King solves it numerically. I solve it using the error function, which can be approximated analytically.
def rho(V):
    #   The piecewise nature of this function is not for numerical reasons, the density really is discontinuous due to galactic
    #   tidal forces carrying away stars that wonder out beyond the cluster radius, which is the key idea behind the King model.
    return 0 if V > 0 else sqrt(pi**3) * k / j**3 * exp(2 * j**2 * (V0 - V)) * erf(j * sqrt(-2 * V))\
                           -2 * pi * k * sqrt(-2 * V) * exp(2 * j**2 * V0) * (1 / j**2 - 4 / 3 * V)

#   This is the right hand side of the partial differential equation describing density as a function of radial distance.
#   This differential equation is the differential version of Gauss's law for gravity in radial coordinates, under spherical
#   symmetry.
def rhs(r, q):
    a = 4 * pi * G * rho(q[0])
    #   An inverse square field has a potential equal to zero at distance 0.
    if r != 0:
        a -= 2 * q[1] / r
    return [q[1], a]



#   Notice the 0 initial condition mentioned earlier. King here would use some undimensionalized experimental value (8 I
#   believe) in his undimensionalized version of the partial differential equation instead.
sol = solve_ivp(rhs, [0, r_max], [V0, 0], t_eval=linspace(0, r_max, N))



#   Plain linear interpolation. Anything else might make the boundary smooth, which is unphysical.
V_interp = interp1d(sol.t, sol.y[0], assume_sorted=True)

#   Check if the range was chosen sufficiently large.
i = 0
while sol.y[0, i] < 0:
    i += 1
    if i == N:
        raise RuntimeError("can't continue since r0 > r_max")

#   Since the interpolation is linear, the exact boundary of the cluster can easily be evaluated.
r0 = sol.t[i-1] - sol.y[0, i-1] * r_max / N / (sol.y[0, i] - sol.y[0, i-1])

if i < THRESHOLD_STEPS_TO_BOUNDARY:
    raise RuntimeWarning('steps to boundary too small; reduce r_max and/or increase N')



#   Rejection sampling the main probability distribution function.
#   The radial distance sampling is simple. The velocity sampling range depends on the radial distance that was just sampled
#   (because e.g. stars far from the cluster center with high gravitational energy must have low kinetic energy).
samples_r = np.empty(n)
samples_v = np.empty(n)
i = 0
attempts = 0
while i < n:
    r = generator.uniform(0, r0)
    v = generator.uniform(0, sqrt(-2 * V_interp(r)))
    p = generator.uniform(0, 1 / np.e**2)
    if p < r**2 * v**2 * exp(-2 * j**2 * (V_interp(r) - V0)) * (exp(-j**2 * v**2) - exp(j**2 * 2 * V_interp(r))):
        samples_r[i] = r
        samples_v[i] = v
        i += 1
    attempts += 1

print(f'sampling efficiency {i / attempts * 100:.2f}%')



#   The sampling above did not include directional data. Position and velocity are isotropically distributed.
samples_r_polar = arccos(generator.uniform(-1, 1, n))
samples_v_polar = arccos(generator.uniform(-1, 1, n))
samples_r_azimuthal = generator.uniform(0, 2 * pi, n)
samples_v_azimuthal = generator.uniform(0, 2 * pi, n)



#   The samples are complete now, but should be transformed from spherical to Cartesian coordinates for storage and direct use
#   by the solvers.
def spherical_to_Cartesian(r, polar, azimuthal):
    return array([r * sin(polar) * cos(azimuthal),
                  r * sin(polar) * sin(azimuthal),
                  r * cos(polar)])
positions  = spherical_to_Cartesian(samples_r, samples_r_polar, samples_r_azimuthal)
velocities = spherical_to_Cartesian(samples_v, samples_v_polar, samples_v_azimuthal)



#   A struct of array was calculated for sampling, but an array of structs is needed for simulation (x, y, z together).
positions  = transpose(positions)
velocities = transpose(velocities)


if file_format == 'vtu':
    N_Body_vtu('king').write(positions, velocities)
elif file_format == 'h5':
    N_Body_h5('king').write(positions, velocities)
