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
Create an HDF5 file (.h5) from an initial condition.
An initial condition can optionally contain velocities and describe multiple times.

These files don't have a standard format, but compress much better due to control
over data storage and having a suitable compressor (Szip).

Groups are named as increasing integers.
Position and velocity data are stored as separate datasets 'pos' and 'vel' in
the first group ('0').
These datasets contain data over multiple times, so are 3rd order tensors of
dimension (# times, # objects, 3).
A 'time' dataset is also stored in each group, containing the times.

Access speed from fast to slow: spatial dimensions, objects, times.

Root gets a string attribute describing much of what is described here.
'''

import h5py
import numpy as np

class N_Body_h5:
    def __init__(self, name: str, n_times: int = 1, time_group_size: int = 2**14, time_chunk_size: int = 2**8):
        self.fp              = h5py.File(name + '.h5', 'w', libver='latest')
        self.time_count      = 0
        #   The only reason # of times is used is because N_Body_h5's dtor can't be relied on to write
        #   the final lingering time data, since Python starts shutting down before HDF5 is done writing.
        #   Knowing # of times in advance, the write function can add the linger times.
        self.n_times         = n_times
        self.times           = np.empty([min(n_times, time_group_size)], dtype=np.float64)
        self.time_group_size = time_group_size
        self.time_chunk_size = time_chunk_size
        #   Number of objects in the previous time step, to check for consistency.
        self.n_prev          = None
        self.had_vel_prev    = None
        #   Writing info attributes.
        self.fp.attrs.create('time chunk size', time_group_size, dtype=np.uint32)
        self.fp.attrs['info'] = \
            '''
            Groups in root represent data in time chunks, and are named as increasing integers.
            The number of time steps per chunk is given by the 'time chunk size' attribute.
            Groups contain a position 'pos' dataset, the 'times' dataset, and optionally a 'vel' velocity dataset.
            '''
    def write(self, positions: np.ndarray, velocities: np.ndarray = None, time: float = 0):
        if velocities is not None:
            if positions.shape != velocities.shape:
                raise BufferError('positions and velocities should have the same shape')
        n = positions.shape[0]
        if self.n_prev is not None and self.n_prev != n:
            raise BufferError("number of objects can't be changing throughout the initial condition")
        if self.had_vel_prev and velocities is None:
            raise BufferError('either all or none of the times must have velocities')
        if self.time_count % self.time_group_size == 0:
            #   Save times.
            if self.time_count != 0:
                self.fp['/' + str(self.time_count // self.time_group_size - 1)].create_dataset('times', data=self.times)
            #   Time for a new group.
            group = self.fp.create_group(str(self.time_count // self.time_group_size))
            group.create_dataset('pos', shape=(0, n, 3), maxshape=(None, n, 3),
                                 chunks=(self.time_chunk_size, n, 3), dtype=np.float64, compression='szip')
            if velocities is not None:
                group.create_dataset('vel', shape=(0, n, 3), maxshape=(None, n, 3),
                                     chunks=(self.time_chunk_size, n, 3), dtype=np.float64, compression='szip')
        else:
            group = self.fp['/' + str(self.time_count // self.time_group_size)]
        self.times[self.time_count % self.time_group_size] = time
        group['pos'].resize(self.time_count % self.time_group_size + 1, axis=0)
        if velocities is not None:
            group['vel'].resize(self.time_count % self.time_group_size + 1, axis=0)
        group['pos'][self.time_count % self.time_group_size] = positions
        if velocities is not None:
            group['vel'][self.time_count % self.time_group_size] = velocities
        if self.time_count + 1 == self.n_times:
            #   Save lingering time steps.
            group.create_dataset('times', data=self.times[:(self.time_count % self.time_group_size + 1)])
        self.n_prev = n
        self.had_vel_prev = velocities is not None
        self.time_count += 1

    def __del__(self):
        self.fp.close()

if __name__ == '__main__':
    generator = np.random.default_rng(0)
    n         = 16
    #   Test without velocities and one time.
    writer    = N_Body_h5('0')
    positions = generator.normal(0, 1, (n, 3))
    writer.write(positions)
    #   Test with velocities over multiple times.
    n_times   = 7
    times     = np.linspace(0, 1, n_times)
    writer    = N_Body_h5('1', n_times, 5)
    for time in times:
        positions  = generator.normal(0, 1, (n, 3))
        velocities = generator.normal(0, 1, (n, 3))
        writer.write(positions, velocities, time)
