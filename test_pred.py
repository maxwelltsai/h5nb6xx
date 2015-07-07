# verify the prediction
import numpy as np
from amuse.lab import *
from amuse.community.h5nb6xx.interface import H5nb6xx
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Mtot=2100|units.MSun
Rvir = 1|units.parsec # ridiculously small cluster, so we see something happening
W0=7.0
n = 3000
converter = nbody_system.nbody_to_si(Mtot,Rvir)

cluster_gravity = H5nb6xx(converter, redirection="none")
cluster_gravity.set_h5_filename('/Users/maxwell/Works/GalevNB/run/data.h5part')
cluster_gravity.commit_particles()
cluster_gravity.initialize_code()
cluster_gravity.get_number_of_particles()
x = []
y = []
z = []
xp = []
yp = []
zp = []
tspace = np.linspace(0,2000000, num=1000)
for tt in tspace:
    tphy = units.yr.new_quantity(tt)
    cluster_gravity.evolve_model(tphy)
    cluster_gravity.synchronize_model()
    tm, tx, ty, tz, tvx, tvy, tvz, tr = cluster_gravity.get_state(2)
    xp.append(tx.value_in(units.AU))
    yp.append(ty.value_in(units.AU))
    zp.append(tz.value_in(units.AU))

cluster_gravity.enableInterpolation(0)
for tt in tspace:
    tphy = units.yr.new_quantity(tt)
    cluster_gravity.evolve_model(tphy)
    tm, tx, ty, tz, tvx, tvy, tvz, tr = cluster_gravity.get_state(2)
    x.append(tx.value_in(units.AU))
    y.append(ty.value_in(units.AU))
    z.append(tz.value_in(units.AU))


cluster_gravity.stop()

plt.plot(tspace, x,'r.', label='x')
plt.plot(tspace, xp, 'r-', label='x (pred)')
plt.plot(tspace, y, 'g.', label='y')
plt.plot(tspace, yp, 'g-', label='y (pred)')
plt.plot(tspace, z, 'b.', label='z')
plt.plot(tspace, zp, 'b-', label='z (pred)')
plt.legend()
plt.xlabel('Time [years]')
plt.ylabel('Length [AU]')
plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(xp,yp,zp)
ax.plot(x,y,z,'go')
plt.show()
