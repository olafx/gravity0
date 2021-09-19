import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import erf
import h5py

k = 1
j = 1
V0 = -1
G = 1
r_max = 1.5
n = 100
N = 20
filename = '0.h5'

def rho(V):
    return 0 if V > 0 else np.pi**1.5*k/j**3*np.exp(2*j**2*(V0-V))*erf(j*np.sqrt(-2*V))\
                         -2*np.pi*k*np.sqrt(-2*V)*np.exp(2*j**2*V0)*(1/(j**2)-4*V/3)

def rhs(r, q):
    return [q[1], 4*np.pi*G*rho(q[0]) if r == 0 else 4*np.pi*G*rho(q[0])-2*q[1]/r]



r_eval = np.linspace(0, r_max, n)
sol = solve_ivp(rhs, [0, r_max], [V0, 0], t_eval=r_eval)

print('integration complete')



V_interp = interp1d(sol.t, sol.y[0], assume_sorted=True)

i = 0
while sol.y[0,i] < 0:
    i += 1
    if i == n:
        exit('could not continue since r0 > r_max')
r0 = sol.t[i-1]-sol.y[0,i-1]*(r_eval[1]-r_eval[0])/(sol.y[0,i]-sol.y[0,i-1])

print(f'cluster boundary is at {r0:.4f}, which is {n*r0/r_max:.2f} steps from the origin')

samples_r = np.empty(N)
samples_v = np.empty(N)
i = 0
attempts = 0
while i < N:
    r = np.random.uniform(0, r0)
    v = np.random.uniform(0, np.sqrt(-2*V_interp(r)))
    p = np.random.uniform(0, 1/np.e**2)
    if p < r**2*v**2*np.exp(-2*j**2*(V_interp(r)-V0))*(np.exp(-j**2*v**2)-np.exp(j**2*2*V_interp(r))):
        samples_r[i] = r
        samples_v[i] = v
        i += 1
    attempts += 1

print(f'r and v sampling completed with efficiency {100*i/attempts:.2f}%')



samples_r_polar = np.arccos(np.random.uniform(-1, 1, N))
samples_v_polar = np.arccos(np.random.uniform(-1, 1, N))
samples_r_azimuthal = np.random.uniform(0, 2*np.pi, N)
samples_v_azimuthal = np.random.uniform(0, 2*np.pi, N)

print('r and v angular sampling complete')



data = np.empty(6*N)

for i in range(N):
    data[6*i]   = samples_r[i]*np.sin(samples_r_polar[i])*np.cos(samples_r_azimuthal[i])
    data[6*i+1] = samples_r[i]*np.sin(samples_r_polar[i])*np.sin(samples_r_azimuthal[i])
    data[6*i+2] = samples_r[i]*np.cos(samples_r_polar[i])

for i in range(N):
    data[6*i+3] = samples_v[i]*np.sin(samples_v_polar[i])*np.cos(samples_v_azimuthal[i])
    data[6*i+4] = samples_v[i]*np.sin(samples_v_polar[i])*np.sin(samples_v_azimuthal[i])
    data[6*i+5] = samples_v[i]*np.cos(samples_v_polar[i])

print('spherical to Cartesian transformation complete')



fp = h5py.File(filename, 'w')
fp.create_dataset('ic', data=data, dtype=np.float64, shape=(2*N, 3))
fp.close()

print('data stored')
