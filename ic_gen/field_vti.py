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

from vtk import *
from vtk.util.numpy_support import *
import numpy as np

'''
Create a vtk XML image file (.vti) from an initial condition.
An initial condition can optionally contain velocities and describe multiple times.

Unfortunately, when the file gets somewhat large, Python VTK kicks the bucket.
It can do up to 256^3 or 4096^2 on my machine, which is only a few hundred Mb.
'''


class Field_vti:

    def __init__(self, name: str, n_times: float = 1):
        self.writer = vtkXMLImageDataWriter()
        self.image = vtkImageData()
        self.densities = vtkDoubleArray()
        self.time_count = 0
        self.n_times = n_times
        self.densities.SetNumberOfComponents(1)
        self.densities.SetName('Density')
        #   Image contains point data (which itself contains densities).
        self.image.GetPointData().SetScalars(self.densities)
        self.writer.SetFileName(name + '.vti')
        self.writer.SetInputData(self.image)
        self.writer.SetNumberOfTimeSteps(self.n_times)

    def write(self, densities: np.ndarray, velocities: np.ndarray = None, time: float = 0):
        if len(densities.shape) not in [2, 3] or velocities is not None and len(velocities.shape) not in [3, 4]:
            raise BufferError('data has the wrong order')
        if velocities is not None and densities.shape != velocities.shape[:-1]:
            raise BufferError('densities and velocities should have a consistent shape')
        if self.time_count == 0:
            #   Need to set the image extent still.
            shape = densities.shape
            print('SHAPE', shape)
            self.image.SetExtent(0, shape[0] - 1, 0, shape[1] - 1, 0, 0 if len(shape) == 2 else shape[2] - 1)
            if velocities is not None:
                #   Can only know now the user wants velocities.
                #   Velocities are represented by point data.
                self.velocities = vtkDoubleArray()
                self.velocities.SetName('Velocity')
                self.velocities.SetNumberOfComponents(len(shape))
                #   Image stores velocities as point data. Need to add an array, can
                #   in general have multiple point data arrays, or 'attributes' as in Xdmf3.
                self.image.GetPointData().AddArray(self.velocities)
            self.writer.Start()
        #   Update data. The 1 signifies not to deallocate. (Makes more sense in C++.)
        #   Needs to be flattened because numpy_to_vtk only works up to order 2. This is because it's supposed to
        #   receive an array of vectors in general, which is in this case an array of scalars.
        self.densities.SetArray(numpy_to_vtk(densities.flatten()), densities.size, 1)
        #   Signal that data has been modified so that new time gets its own data
        #   rather than referencing the data from the previous time.
        self.densities.Modified()
        if velocities is not None:
            #   Needs to be flattened also, but this time to obtain an array of size 3 vectors.
            self.velocities.SetArray(numpy_to_vtk(velocities.reshape(velocities.size // len(shape), len(shape))), velocities.size, 1)
            self.velocities.Modified()
        self.writer.WriteNextTime(time)
        self.time_count += 1
        if self.time_count == self.n_times:
            self.writer.Stop()


if __name__ == '__main__':

    generator = np.random.default_rng(0)
    n = 16

    #   Test with velocities over multiple times.
    n_times = 4
    times = np.linspace(0, 1, n_times)
    writer = Field_vti('0', n_times)
    for time in times:
        densities = generator.normal(0, 1, (n, n))
        velocities = generator.normal(0, 1, (n, n, 3))
        writer.write(densities, velocities, time)
