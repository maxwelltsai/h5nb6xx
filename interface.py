from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class H5nb6xxInterface(CodeInterface,
        GravitationalDynamicsInterface,
        GravityFieldInterface):

    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="h5nb6xx_worker", **keyword_arguments)


    @legacy_function
    def get_h5_filename():
        function = LegacyFunctionSpecification()
        function.addParameter('h5_filename', dtype='string', direction=function.OUT,
            description = "The fullpath of the h5/h5part file.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_h5_filename():
        function = LegacyFunctionSpecification()
        function.addParameter('h5_filename', dtype='string', direction=function.IN,
            description = "The fullpath of the h5/h5part file.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_total_number_of_steps():
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='int32', direction=function.OUT,
            description = "The total number of steps in the h5/h5part file.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_total_number_of_steps():
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='int32', direction=function.IN,
            description = "The total number of steps in the h5/h5part file.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_duration_per_step():
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='float64', direction=function.OUT,
            description = "The duration between this step and the next step in N-body units.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_duration_per_step():
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='float64', direction=function.IN,
            description = "The duration between this step and the next step in N-body units.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eta():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def enableInterpolation():
        """
        Enables prediction/interpolation of coordinates.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('val', dtype='int32',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_host_star_flag():
        """
        Set a flag to determine whether a given ID is a host star.
        A host star is usually excluded in the get_gravity_at_point() calculations.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('host_star_id', dtype='int32',
                              direction=function.IN)
        function.addParameter('flag', dtype='int32',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_neighbors():
        """
        Get the IDs of the closest stars.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('host_star_id', dtype='int32',
                              direction=function.IN)
        function.addParameter('neighbor_star_id', dtype='int32',
                              direction=function.OUT)
        function.addParameter('nth_neighbor', dtype='int32',
                              direction=function.IN)
        function.addParameter('number_of_stars_inquired', dtype='int32',
                              direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_states():
        """
        Get the states of the specified particles.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The mass the particle")
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.addParameter('y', dtype='float64', direction=function.OUT)
        function.addParameter('z', dtype='float64', direction=function.OUT)
        function.addParameter('vx', dtype='float64', direction=function.OUT)
        function.addParameter('vy', dtype='float64', direction=function.OUT)
        function.addParameter('vz', dtype='float64', direction=function.OUT)
        function.addParameter('radius', dtype='float64', direction=function.OUT)
        function.addParameter('number_of_stars_inquired', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_state_a_adot():
        """
        Get the a and a dot of the specified particles.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('ax', dtype='float64', direction=function.OUT)
        function.addParameter('ay', dtype='float64', direction=function.OUT)
        function.addParameter('az', dtype='float64', direction=function.OUT)
        function.addParameter('jx', dtype='float64', direction=function.OUT)
        function.addParameter('jy', dtype='float64', direction=function.OUT)
        function.addParameter('jz', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_states_a_adot():
        """
        Get the a and a dot of the specified particles.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('ax', dtype='float64', direction=function.OUT)
        function.addParameter('ay', dtype='float64', direction=function.OUT)
        function.addParameter('az', dtype='float64', direction=function.OUT)
        function.addParameter('jx', dtype='float64', direction=function.OUT)
        function.addParameter('jy', dtype='float64', direction=function.OUT)
        function.addParameter('jz', dtype='float64', direction=function.OUT)
        function.addParameter('number_of_stars_inquired', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
            particle could not be found
        """
        return function


class H5nb6xx(GravitationalDynamics, GravityFieldCode):

    #def __init__(self, **options):
    #    InCodeComponentImplementation.__init__(self,  H5nb6xxInterface(**options), **options)

    def __init__(self, convert_nbody = None, **kargs):
        GravitationalDynamics.__init__(self,  H5nb6xxInterface(**kargs), convert_nbody, **kargs)

    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        object.add_method('RUN', 'get_particle_timestep')
        object.add_method('RUN', 'get_state_a_adot')
        GravityFieldCode.define_state(self, object)

        object.add_method('EDIT', 'set_state')
        object.add_method('EDIT', 'set_velocity')
        object.add_method('EDIT', 'set_mass')
        object.add_method('EDIT', 'set_position')
        object.add_method('CHANGED','before_get_parameter')

        object.add_transition('RUN', 'CHANGED', 'set_state', False)
        object.add_transition('RUN', 'CHANGED', 'set_velocity', False)
        object.add_transition('RUN', 'CHANGED', 'set_mass', False)
        object.add_transition('RUN', 'CHANGED', 'set_position', False)
        object.add_transition('CHANGED', 'RUN', 'synchronize_model')
        object.add_method('CHANGED', 'get_state')
        object.add_method('CHANGED', 'get_mass')
        object.add_method('CHANGED', 'get_position')
        object.add_method('CHANGED', 'get_velocity')
        object.add_method('CHANGED', 'get_particle_timestep')

    def define_parameters(self, object):

        # Set/get parameters specific to the module, not part of the
        # standard interface.  Accessors used here must be defined
        # above and reflected in interface.cc.  Python access is
        # (e.g.)
        #
        #        ph4.parameters.timestep_parameter = xxx

        object.add_method_parameter(
            "get_eta",                   # getter name in interface.cc
            "set_eta",                   # setter name in interface.cc
            "timestep_parameter",        # python parameter name
            "timestep parameter",        # description
            default_value = 0.05
        )

        object.add_method_parameter(
            "get_eps2",                  # already defined in standard interface
            "set_eps2",                  # already defined in standard interface
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )


        object.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )

        object.add_method_parameter(
            "get_time",
            None,
            "time",
            "current simulation time",
            default_value = 0.0
        )
        #self.stopping_conditions.define_parameters(object)

    def define_properties(self, object):
        object.add_property("get_kinetic_energy")
        object.add_property("get_potential_energy")
        object.add_property("get_total_energy")
        object.add_property("get_center_of_mass_position")
        object.add_property("get_center_of_mass_velocity")
        object.add_property("get_total_mass")

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        # Turn interface functions into methods.

        object.add_method(
            "set_eps2",
            (
                nbody_system.length * nbody_system.length
            ),
            (
                object.ERROR_CODE
            )
        )

        object.add_method(
            "get_eps2",
            (),
            (
                nbody_system.length * nbody_system.length,
                object.ERROR_CODE
            )
        )

        object.add_method(
            "get_state_a_adot",
            (
                object.NO_UNIT,
            ),
            (
                nbody_system.speed/nbody_system.time,
                nbody_system.speed/nbody_system.time,
                nbody_system.speed/nbody_system.time,
                nbody_system.speed/nbody_system.time/nbody_system.time,
                nbody_system.speed/nbody_system.time/nbody_system.time,
                nbody_system.speed/nbody_system.time/nbody_system.time,
                object.ERROR_CODE
            )
        )

