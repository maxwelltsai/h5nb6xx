from amuse.lab import *
from amuse.community.h5nb6xx.interface import H5nb6xx
from matplotlib import pyplot as plt

Mtot=6.411723357779E-01*4000|units.MSun
Rvir = 1|units.parsec # ridiculously small cluster, so we see something happening
W0=7.0
n = 3000
converter = nbody_system.nbody_to_si(Mtot,Rvir)

cluster_gravity = H5nb6xx(converter, redirection="none")

# cluster_gravity.set_h5_filename('/Users/maxwell/Works/GalevNB/run/data.h5part')
cluster_gravity.set_h5_filename('/Users/maxwell/Works/visnb6/run/data.h5part')
cluster_gravity.commit_particles()
cluster_gravity.initialize_code()
cluster_gravity.get_number_of_particles()

cluster_gravity.set_host_star_flag(2277,1)

#t = 0 | units.Myr
dt = 10.| units.yr
t_end = 10*dt
cluster_gravity.evolve_model(1*t_end)
#cluster_gravity.evolve_model(0)
print cluster_gravity.total_mass

#cluster_gravity.get_gravity_at_point(0|units.parsec, 0.5|units.parsec, 0.5|units.parsec, 0.5|units.parsec)

print cluster_gravity.get_acceleration(1)
m, x, y, z, vx, vy, vz, r = cluster_gravity.get_state(1)
print x, y, z
cluster_gravity.set_host_star_flag(1,1)
print cluster_gravity.get_gravity_at_point(0|units.parsec, x+(0.001|units.parsec), y+(0.001|units.parsec), z+(0.001|units.parsec))
cluster_gravity.set_host_star_flag(1,0)
print cluster_gravity.get_gravity_at_point(0|units.parsec, x+(0.001|units.parsec), y+(0.001|units.parsec), z+(0.001|units.parsec))
cluster_gravity.set_host_star_flag(1,1)
print cluster_gravity.get_gravity_at_point(0|units.parsec, x+(0.001|units.parsec), y+(0.001|units.parsec), z+(0.001|units.parsec))
cluster_gravity.set_host_star_flag(1,0)
print cluster_gravity.get_gravity_at_point(0|units.parsec, x+(0.001|units.parsec), y+(0.001|units.parsec), z+(0.001|units.parsec))

cluster_gravity.stop()
