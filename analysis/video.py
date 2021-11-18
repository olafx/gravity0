import sys
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim

if len(sys.argv) != 3:
    exit()

data_name = sys.argv[1]
video_name = sys.argv[2]

n = 1000

ffmpeg_writer = anim.writers['ffmpeg']
writer = ffmpeg_writer(fps=60)

fig = plt.figure(figsize=(5.4, 5.4))
ax = fig.add_subplot(projection='3d')

im, = ax.plot([], [], [], 'o')

file = h5py.File(data_name, 'r')

with writer.saving(fig, video_name, dpi=100):
    for i in range(n + 1):
        data = file[str(i)]
        im.set_data(data[:,0], data[:,1])
        im.set_3d_properties(data[:,2])

        x_min = np.min(data[:,0])
        x_max = np.max(data[:,0])
        y_min = np.min(data[:,1])
        y_max = np.max(data[:,1])
        z_min = np.min(data[:,2])
        z_max = np.max(data[:,2])

        ax.set_xlim(1.2*x_min, 1.2*x_max)
        ax.set_ylim(1.2*y_min, 1.2*y_max)
        ax.set_zlim(1.2*z_min, 1.2*z_max)

        writer.grab_frame()
        print(f'{i}/{n}')

file.close()
