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
Create a vtk XML unstructured grid file (.vtu) from an initial condition.
An initial condition can optionally contain velocities and describe multiple times.
'''


class N_Body_vtu:

    def __init__(self, name: str, n_times: float = 1):
        self.writer = vtkXMLUnstructuredGridWriter()
        self.grid = vtkUnstructuredGrid()
        self.points = vtkPoints()
        #   Points are represented by positions.
        self.positions = vtkDoubleArray()
        self.time_count = 0
        self.n_times = n_times
        self.positions.SetNumberOfComponents(3)
        #   Grid contains points (which themselves contain positions).
        self.points.SetData(self.positions)
        self.grid.SetPoints(self.points)
        self.writer.SetFileName(name + '.vtu')
        self.writer.SetInputData(self.grid)
        self.writer.SetNumberOfTimeSteps(n_times)

    def write(self, positions: np.ndarray, velocities: np.ndarray = None, time: float = 0):
        if velocities is not None:
            if positions.size != velocities.size:
                raise BufferError('positions and velocities should have equal size')
        if self.time_count == 0:
            if velocities is not None:
                #   Can only know now the user wants velocities.
                #   Velocities are represented by point data.
                self.velocities = vtkDoubleArray()
                self.velocities.SetName('Velocity')
                self.velocities.SetNumberOfComponents(3)
                #   Grid stores velocities as point data. Need to add an array, can
                #   in general have multiple point data arrays, or 'attributes' as in Xdmf3.
                self.grid.GetPointData().AddArray(self.velocities)
            self.writer.Start()
        #   Update data. The 1 signifies not to deallocate. (Makes more sense in C++.)
        self.positions.SetArray(numpy_to_vtk(positions), positions.size, 1)
        #   Signal that data has been modified so that new time gets its own data
        #   rather than referencing the data from the previous time.
        self.positions.Modified()
        if velocities is not None:
            self.velocities.SetArray(numpy_to_vtk(velocities), velocities.size, 1)
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
    writer = N_Body_vtu('1', n_times)
    for time in times:
        positions = generator.normal(0, 1, (n, 3))
        velocities = generator.normal(0, 1, (n, 3))
        writer.write(positions, velocities, time)

    #   Test without velocities and one time.
    positions = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]])
    N_Body_vtu('0').write(positions)
