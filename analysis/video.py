import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as anim

n = 20000

ffmpeg_writer = anim.writers['ffmpeg']
writer = ffmpeg_writer(fps=60)

ax = (fig := plt.figure(figsize=(5.4, 5.4))).add_subplot(111)
ax.set_aspect('equal')

im, = ax.plot([], [], 'o')

file = h5py.File('gravity0.h5', 'r')

with writer.saving(fig, "gravity0.mp4", dpi=100):
    for i in range(n + 1):
        data = file[str(i)]
        im.set_data(data[:,0], data[:,1])

        x_min = np.min(data[:,0])
        x_max = np.max(data[:,0])
        y_min = np.min(data[:,1])
        y_max = np.max(data[:,1])

        if (x_max-x_min > y_max-y_min):
            ax.set_xlim(1.2*x_min, 1.2*x_max), ax.set_ylim(1.2*((y_min+y_max)/2-(x_max-x_min)/2), 1.2*((y_min+y_max)/2+(x_max-x_min)/2))
        else:
            ax.set_xlim(1.2*((x_min+x_max)/2-(y_max-y_min)/2), 1.2*((x_min+x_max)/2+(y_max-y_min)/2)), ax.set_ylim(1.2*y_min, 1.2*y_max)

        writer.grab_frame()
        print(i, '/', n)

file.close()
