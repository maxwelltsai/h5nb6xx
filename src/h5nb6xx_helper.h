#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <float.h>

#ifndef _H5NB6XX_HELPER_H
#define _H5NB6XX_HELPER_H

#include "cuda_util.h"
class CUDA_Util;

class H5nb6xx_Helper 
{

    public:
        struct Dynamics {
            int *id; // sorted ID, id[NAME[index_of_part]] is the internal id of the particle
            int *id_original; // orig. ID read from the HDF5, given internal id, get NAME
            float *x;
            float *y;
            float *z;
            float *vx;
            float *vy;
            float *vz;
            float *ax;
            float *ay;
            float *az;
            float *jx;
            float *jy;
            float *jz;
            float *mass;
            double time;
            int step_id;
            hid_t h5_group_id;
            int n_records; // total number of particles recorded in a step (group)

            ~Dynamics() {
                printf("destructor start\n");
                if(id!=NULL) delete [] id;
                if(id!=NULL) delete [] id_original;
                if(x!=NULL) delete [] x;
                if(y!=NULL) delete [] y;
                if(z!=NULL) delete [] z;
                if(vx!=NULL) delete [] vx;
                if(vy!=NULL) delete [] vy;
                if(vz!=NULL) delete [] vz;
                if(ax!=NULL) delete [] ax;
                if(ay!=NULL) delete [] ay;
                if(az!=NULL) delete [] az;
                if(jx!=NULL) delete [] jx;
                if(jy!=NULL) delete [] jy;
                if(jz!=NULL) delete [] jz;
                if(mass!=NULL) delete [] mass;
                id = NULL;
                x = NULL;
                y = NULL;
                z = NULL;
                vx = NULL;
                vy = NULL;
                vz = NULL;
                ax = NULL;
                ay = NULL;
                az = NULL;
                jx = NULL;
                jy = NULL;
                jz = NULL;
                mass = NULL;
                printf("destructor end\n");
            }
        }; // storage of the dynamical data for current step

        struct Status{
            char h5_filename[256];
            double current_time; // current time passed by evolve_model()
            double begin_time; // begin time of the HDF5 snapshot
            double end_time; // end time of the HDF5 snapshot
            int current_step_id; // ID of the current step
            int prev_step_id; // ID of the previous step
            int next_step_id; // ID of the next step
            hid_t h5_file_id; // file ID of the opened HDF5 snapshot
            int nsteps; // total number of steps
            double t_step; //during between to steps
            // number of particles for the current step. Since this interface requires all step have the same
            // vector lengths, if (n_particles != data->n_records), the code will stop, otherwise the result
            // would be unpreditable.
            int n_particles; // total number of particles in the star cluster (determined by the first step).
            Dynamics* data;
            Dynamics* next_data;
            Dynamics* interp_data;

        }; // storage of global properties

        static H5nb6xx_Helper* GetInstance();
        H5nb6xx_Helper();

        Status get_status();
        Dynamics* get_data();
        Dynamics* get_data_next();
        // AMUSE helper functions
        int helper_get_h5_filename(char **h5_filename);
        int helper_set_h5_filename(char *h5_filename);
        int helper_get_eps2(double * epsilon_squared);
        int helper_set_eps2(double epsilon_squared);
        int helper_initialize_code();
        int helper_new_particle(int *index_of_the_particle, double mass, 
                double x, double y, double z, double vx, double vy, double vz, double radius);
        int helper_delete_particle(int index_of_the_particle);
        int helper_get_number_of_particles(int * number_of_particles);
        int helper_get_index_of_first_particle(int * index_of_the_particle);
        int helper_get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle);
        int helper_get_state(int index_of_the_particle, double * mass, double * x, 
                double * y, double * z, double * vx, double * vy, double * vz, double * radius);
        int helper_get_state_a_adot(int index_of_the_particle, double * ax, 
                double * ay, double * az, double * jx, double * jy, double * jz);
        int helper_get_states(int* index_of_the_particle, double * mass, double * x, 
                double * y, double * z, double * vx, double * vy, double * vz, double * radius, int n);
        int helper_get_states_a_adot(int* index_of_the_particle, double * ax, 
                double * ay, double * az, double * jx, double * jy, double * jz, int n);
        int helper_set_state(int index_of_the_particle, double mass, double x, 
                double y, double z, double vx, double vy, double vz, double radius);
        int helper_get_mass(int index_of_the_particle, double * mass);
        int helper_set_mass(int index_of_the_particle, double mass);
        int helper_get_position(int index_of_the_particle, double * x, double * y, double * z);
        int helper_set_position(int index_of_the_particle, double x, double y, double z);
        int helper_set_acceleration(int index_of_the_particle, double ax, double ay, double az);
        int helper_get_acceleration(int index_of_the_particle, double * ax, double * ay, double * az);
        int helper_get_potential(int index_of_the_particle, double * potential);
        int helper_commit_particles();
        int helper_get_time(double * time);
        int helper_get_kinetic_energy(double * kinetic_energy);
        int helper_get_potential_energy(double * potential_energy);
        int helper_get_center_of_mass_velocity(double * vx, double * vy, double * vz);
        int helper_get_center_of_mass_position(double * x, double * y, double * z);
        int helper_get_total_mass(double * mass);
        int helper_get_total_radius(double * radius);
        int helper_get_radius(int index_of_the_particle, double * radius);
        int helper_set_radius(int index_of_the_particle, double radius);
        int helper_cleanup_code();
        int helper_evolve_model(double to_time);
        int helper_get_velocity(int index_of_the_particle, double * vx, double * vy, double * vz);
        int helper_set_velocity(int index_of_the_particle, double vx, double vy, double vz);
        int helper_get_time_step(double * time_step);
        int helper_get_begin_time(double * output);
        int helper_set_begin_time(double input);
        int helper_commit_parameters();
        int helper_synchronize_model();
        int helper_recommit_particles();
        int helper_recommit_parameters();
        int helper_get_potential_at_point(double * eps, double * x, double * y, double * z, double * phi, int n);
        int helper_get_gravity_at_point(double * eps, double * x, double * y, double * z,
                double * forcex, double * forcey, double * forcez, int n);
        int helper_get_total_number_of_steps(int *val);
        int helper_set_total_number_of_steps(int val);
        int helper_get_duration_per_step(double *val);
        int helper_set_duration_per_step(double val);
        int helper_set_eta(double timestep_parameter);
        int helper_get_eta(double * timestep_parameter);
        int helper_set_enable_interpolation(bool val);
        int helper_set_host_star_flag(int host_star_id, int flag);
        int helper_get_neighbors(int *host_star_id, int* neighbor_star_id, int *nth_neighbor, int n);




    protected:
        herr_t HDF5_error;
        int MAX_PARTICLE_NUMBER;
        double _TINY_;
        int next_particle_id;
        int* host_star_flag; // 0 if not a host star, 1 if otherwise



        // HDF5 utility functions
        herr_t h5_op_func (hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data);
        int h5_get_total_number_of_groups(hid_t parent_id);
        hid_t h5_open_group_by_name(hid_t parent_id, const char* group_name);
        double h5_read_attribute_double(hid_t parent_id, const char *attrib_name);
        int h5_read_attribute_integer(hid_t parent_id, const char *attrib_name);
        int *h5_read_dataset_as_integer_vector(hid_t parent_id, const char *dset_name);
        float *h5_read_dataset_as_float_vector(hid_t parent_id, const char *dset_name);
        double h5_get_step_duration();
        int h5_get_dataset_vector_length(hid_t parent_id, const char *dset_name);
        int h5_load_step_by_id(int step_id, Dynamics* data);
        bool enable_interpolation;

        static H5nb6xx_Helper* instance;
        Status status;

#ifdef GPU
        CUDA_Util* cuda_util;
#endif
};

#endif
