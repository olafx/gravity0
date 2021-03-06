"""
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
"""

import sys
import numpy
import scipy.fft

'''
Generate a square or cubic Gaussian random field with a negative power law as power spectrum.

In 3D, with power -2, this produces a roughly scale invariant spectrum, physically similar to quantum fluctuations, and
so is ideal as a field based cosmological initial condition depicting local density variation. For n-body solvers, this
can also be used, but need to follow with the Zel'dovich approximation to get particle position offsets from the density
variations (see doi: 10.1063/1.4822978).

e.g. python field_cosmological.py 2 4096 -1.3 plot store ic
e.g. python field_cosmological.py 3 256 -2 plot
'''

generator = numpy.random.default_rng(0)

def field_cosmological(size: int, n_dims: int, power: float):
    if n_dims == 2:
        k_i = numpy.fft.fftshift(numpy.mgrid[:size, :size] - (size + 1) // 2)
        amplitude = (k_i[0]**2 + k_i[1]**2) ** (.5 * power)
        del k_i
        amplitude[0,0] = 0
        noise = generator.normal(size=(size, size)) \
         + 1j * generator.normal(size=(size, size))
        return numpy.fft.ifft2(noise * amplitude).real
    elif n_dims == 3:
        k_i = numpy.fft.fftshift(numpy.mgrid[:size, :size, :size] - (size + 1) // 2)
        amplitude = (k_i[0]**2 + k_i[1]**2 + k_i[2]**2) ** (.5 * power)
        del k_i
        amplitude[0,0,0] = 0
        noise = generator.normal(size=(size, size, size)) \
         + 1j * generator.normal(size=(size, size, size))
        return numpy.fft.ifftn(noise * amplitude).real


if __name__ == '__main__':
    if len(sys.argv) < 3:
        exit(1)
    n_dims = int(sys.argv[1])
    size = int(sys.argv[2])
    power = float(sys.argv[3])
    density = field_cosmological(size, n_dims, power)
    if 'plot' in sys.argv:
        import matplotlib.pyplot as plt
        plt.imshow(density if n_dims == 2 else density[size // 2])
        plt.title(f'power {power}, region {size}^{n_dims}')
        plt.show()
    if 'store' in sys.argv:
        from field_vti import Field_vti
        velocity = numpy.zeros([size, size, size, n_dims] if n_dims == 3 else [size, size, n_dims]);
        Field_vti(sys.argv[sys.argv.index('store') + 1]).write(density, velocity)
