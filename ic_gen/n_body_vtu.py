'''
Create a vtk XML unstructured grid file (.vtu) from an initial condition constituting positions and velocities.
'''

#   TODO multiple times

from vtk import *
from vtk.util.numpy_support import *
import numpy as np

class N_Body_vtu():
    def __init__(self, name):
        self.writer = vtkXMLUnstructuredGridWriter()
        self.grid = vtkUnstructuredGrid()
        self.points = vtkPoints()
        self.writer.SetInputData(self.grid)
        self.writer.SetFileName(name + '.vtu')
    def write(self, positions, velocities):
        self.points.SetData(numpy_to_vtk(positions, deep=True))
        self.grid.SetPoints(self.points)
        self.writer.Write()

def main():
    writer = N_Body_vtu('n_body_vtu')
    positions, velocities = np.random.rand(10, 3), np.random.rand(10, 3)
    writer.write(positions, velocities)

main()
