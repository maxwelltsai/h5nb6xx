#include "src/h5nb6xx_helper.h"
#include "src/cuda_util.h"
///#include <helper_cuda.h>

/*
 * Interface code
 */


static H5nb6xx_Helper* helper = new H5nb6xx_Helper;


int get_h5_filename(char **h5_filename){
    return helper->helper_get_h5_filename(h5_filename);
}

int set_h5_filename(char *h5_filename){
    return helper->helper_set_h5_filename(h5_filename);;
}

int get_eps2(double * epsilon_squared){
    return helper->helper_get_eps2(epsilon_squared);
}

int set_eps2(double epsilon_squared){
    return helper->helper_set_eps2(epsilon_squared);
}

int initialize_code(){
    return helper->helper_initialize_code();
}

int new_particle(int * index_of_the_particle, double mass, double x, double y, double z, double vx, double vy, double vz, double radius) {
    return helper->helper_new_particle(index_of_the_particle, mass, x, y, z, vx, vy, vz, radius);
}


int delete_particle(int index_of_the_particle) {
    return helper->helper_delete_particle(index_of_the_particle);
}

int get_number_of_particles(int * number_of_particles){
    return helper->helper_get_number_of_particles(number_of_particles);
}

int get_index_of_first_particle(int * index_of_the_particle){
    return helper->helper_get_index_of_first_particle(index_of_the_particle);
}


int get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle) {
    return helper->helper_get_index_of_next_particle(index_of_the_particle, index_of_the_next_particle);
}

int get_state(int index_of_the_particle, double * mass, double * x, double * y, double * z, double * vx, double * vy, double * vz, double * radius){
    return helper->helper_get_state(index_of_the_particle, mass, x, y, z, vx, vy, vz, radius);
}

int get_state_a_adot(int index_of_the_particle, double * ax, double * ay, double * az, double * jx, double * jy, double * jz){
    return helper->helper_get_state_a_adot(index_of_the_particle, ax, ay, az, jx, jy, jz);
}

int get_states(int* index_of_the_particle, double * mass, double * x, double * y, double * z, double * vx, double * vy, double * vz, double * radius, int n){
    return helper->helper_get_states(index_of_the_particle, mass, x, y, z, vx, vy, vz, radius, n);
}

int get_states_a_adot(int* index_of_the_particle, double * ax, double * ay, double * az, double * jx, double * jy, double * jz, int n){
    return helper->helper_get_states_a_adot(index_of_the_particle, ax, ay, az, jx, jy, jz, n);
}

int set_state(int index_of_the_particle, double mass, double x, double y, double z, double vx, double vy, double vz, double radius){
    return helper->helper_set_state(index_of_the_particle, mass, x, y, z, vx, vy, vz, radius);
}

int get_mass(int index_of_the_particle, double * mass) {
    return helper->helper_get_mass(index_of_the_particle, mass);
}

int set_mass(int index_of_the_particle, double mass) {
    return helper->helper_set_mass(index_of_the_particle, mass);
}

int get_position(int index_of_the_particle, double * x, double * y, double * z) {
    return helper->helper_get_position(index_of_the_particle, x, y, z);

}int set_position(int index_of_the_particle, double x, double y, double z) {
    return helper->helper_set_position(index_of_the_particle, x, y, z);
}

int set_acceleration(int index_of_the_particle, double ax, double ay, double az) {
    return helper->helper_set_acceleration(index_of_the_particle, ax, ay, az);
}


int get_acceleration(int index_of_the_particle, double * ax, double * ay, double * az) {
    return helper->helper_get_acceleration(index_of_the_particle, ax, ay, az);
}

int get_potential(int index_of_the_particle, double * potential) {
    return helper->helper_get_potential(index_of_the_particle, potential);
}

int commit_particles() {
    return helper->helper_commit_particles();
}

int get_time(double * time) {
    return helper->helper_get_time(time);
}


int get_kinetic_energy(double * kinetic_energy) {
    return helper->helper_get_kinetic_energy(kinetic_energy);
}

int get_potential_energy(double * potential_energy) {
    return helper->helper_get_potential_energy(potential_energy);
}


int get_center_of_mass_velocity(double * vx, double * vy,  double * vz) {
    return helper->helper_get_center_of_mass_velocity(vx, vy, vz);
}

int get_center_of_mass_position(double * x, double * y, double * z) {
    return helper->helper_get_center_of_mass_position(x, y, z);
}


int get_total_mass(double * mass) {
    return helper->helper_get_total_mass(mass);
}

int get_total_radius(double * radius) {
    return helper->helper_get_total_radius(radius);
}

int get_radius(int index_of_the_particle, double * radius) {
    return helper->helper_get_radius(index_of_the_particle, radius);
}


int set_radius(int index_of_the_particle, double radius) {
    return helper->helper_set_radius(index_of_the_particle, radius);
}

int cleanup_code() {
    return helper->helper_cleanup_code();
}

int evolve_model(double to_time) {
    return helper->helper_evolve_model(to_time);
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, double * vz) {
    return helper->helper_get_velocity(index_of_the_particle, vx, vy, vz);
}


int set_velocity(int index_of_the_particle, double vx, double vy, double vz) {
    return helper->helper_set_velocity(index_of_the_particle, vx, vy, vz);
}

int get_time_step(double * time_step) {
    return helper->helper_get_time_step(time_step);
}


int get_begin_time(double * output) { 
    return helper->helper_get_begin_time(output);
}

int set_begin_time(double input) {
    return helper->helper_set_begin_time(input);
}

int commit_parameters()
{
    // Perform any needed setup after initial code parameters have been set.

    // Consistency check:

    return helper->helper_commit_parameters();
}

int synchronize_model()
{
    // Synchronize all particles at the current system time.  The
    // default is not to reinitialize the scheduler, as this will be
    // handled later, in recommit_particles().
    return helper->helper_synchronize_model();
}

int recommit_particles()
{
    // Reinitialize/reset the system after particles have been added
    // or removed.  The system should be synchronized at some reasonable
    // system_time, so we just need to recompute forces and update the
    // GPU and scheduler.  Note that we don't resize the jdata or
    // idata arrays.  To resize idata, just delete and create a new
    // one.  Resizing jdata is more complicated -- defer for now.

    return helper->helper_recommit_particles();
}


int recommit_parameters()
{
    // Perform any needed changes after code parameters have been reset.

    return helper->helper_recommit_parameters();

}


int get_potential_at_point(double * eps, double * x, double * y, double * z, double * phi, int n) {
    return helper->helper_get_potential_at_point(eps, x, y, z, phi, n);
}

int get_gravity_at_point(double * eps, double * x, double * y, double * z, double * forcex, double * forcey, double * forcez, int n) {
    return helper->helper_get_gravity_at_point(eps, x, y, z, forcex, forcey, forcez, n);
}

int get_total_number_of_steps(int *val) {
    return helper->helper_get_total_number_of_steps(val);
}

int set_total_number_of_steps(int val){
    return helper->helper_set_total_number_of_steps(val);
}


int get_duration_per_step(double *val){
    return helper->helper_get_duration_per_step(val);
}

int set_duration_per_step(double val) {
    return helper->helper_set_duration_per_step(val);
}

int set_eta(double timestep_parameter){
    return helper->helper_set_eta(timestep_parameter);
}

int get_eta(double * timestep_parameter){
    return helper->helper_get_eta(timestep_parameter);
}

int enableInterpolation(int val){
    if(!val){
        return helper->helper_set_enable_interpolation(false);
    } else {

        return helper->helper_set_enable_interpolation(true);
    }
}

int set_host_star_flag(int host_star_id, int flag){
    return helper->helper_set_host_star_flag(host_star_id, flag);
}

int get_neighbors(int *host_star_id, int *neighbor_star_id, int* nth_neighbor, int n){
    return helper->helper_get_neighbors(host_star_id, neighbor_star_id, nth_neighbor, n);
}


