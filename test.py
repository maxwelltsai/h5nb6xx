from interface import H5nb6xx
from amuse.units import nbody_system
from amuse.units import units

inst = H5nb6xx(redirection="none")
inst.set_h5_filename('/Users/maxwell/Works/GalevNB/run/data.h5part')
print inst.get_h5_filename()
inst.initialize_code()
print inst.get_time()
print inst.get_number_of_particles()
#inst.evolve_model(1|nbody_system.time)
inst.evolve_model(1|nbody_system.time)

#1
print inst.get_potential_at_point(0|nbody_system.length, -0.206|nbody_system.length, -0.05|nbody_system.length, 0.306|nbody_system.length)

#20
print inst.get_potential_at_point(0|nbody_system.length, 0.487|nbody_system.length, -0.119|nbody_system.length, 0.028|nbody_system.length)

#1348
print inst.get_potential_at_point(0|nbody_system.length, -0.129|nbody_system.length, 0.296|nbody_system.length, -0.191|nbody_system.length)

#1349
print inst.get_potential_at_point(0|nbody_system.length, -1.118|nbody_system.length, -0.285|nbody_system.length, -0.866|nbody_system.length)



#1
print inst.get_gravity_at_point(0|nbody_system.length, -0.206|nbody_system.length, -0.05|nbody_system.length, 0.306|nbody_system.length)

#20
print inst.get_gravity_at_point(0|nbody_system.length, 0.487|nbody_system.length, -0.119|nbody_system.length, 0.028|nbody_system.length)

#1348
print inst.get_gravity_at_point(0|nbody_system.length, -0.129|nbody_system.length, 0.296|nbody_system.length, -0.191|nbody_system.length)

#1349
print inst.get_gravity_at_point(0|nbody_system.length, -1.118|nbody_system.length, -0.285|nbody_system.length, -0.866|nbody_system.length)
