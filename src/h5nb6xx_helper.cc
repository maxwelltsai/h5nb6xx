#include "h5nb6xx_helper.h"

H5nb6xx_Helper::H5nb6xx_Helper(void){
    HDF5_error = -1;
    MAX_PARTICLE_NUMBER = 81920;
    _TINY_ = 0.0001;
    next_particle_id = 0;
    enable_interpolation = true;
#ifdef GPU
    this->cuda_util = new CUDA_Util(this);
#endif
}

H5nb6xx_Helper* H5nb6xx_Helper::GetInstance(){
    static H5nb6xx_Helper* instance;
    if(!instance){
        instance = new H5nb6xx_Helper();
    }
    return instance;
}

H5nb6xx_Helper::Status H5nb6xx_Helper::get_status() {
    return this->status;
}

H5nb6xx_Helper::Dynamics* H5nb6xx_Helper::get_data() {
    return this->status.data;
}

H5nb6xx_Helper::Dynamics* H5nb6xx_Helper::get_data_next() {
    return this->status.next_data;
}

/************************************************************
  Operator function for H5Ovisit.  This function prints the
  name and type of the object passed to it.
 ************************************************************/
herr_t H5nb6xx_Helper::h5_op_func (hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data) {
    if (name[0] == '.') {         /* Skip counting the root group */
        //printf ("  (Group)\n");
    } else {
        switch (info->type) {
            case H5O_TYPE_GROUP:
                //printf ("%s  (Group)\n", name);
                status.nsteps++;
                break;
            case H5O_TYPE_DATASET:
                //printf ("%s  (Dataset)\n", name);
                break;
            case H5O_TYPE_NAMED_DATATYPE:
                //printf ("%s  (Datatype)\n", name);
                break;
            default:
                printf ("%s  (Unknown)\n", name);
            }
        }
    return 0;
}

int H5nb6xx_Helper::h5_get_total_number_of_groups(hid_t parent_id) {
    hid_t err;
    hsize_t ngroups;
    //err = H5Ovisit (parent_id, H5_INDEX_NAME, H5_ITER_NATIVE, h5_op_func, NULL);
    err = H5Gget_num_objs(parent_id, &ngroups);
    if(err>=0){
        status.nsteps =(int) ngroups;
        printf("Total number of groups: %d\n", status.nsteps);
        return status.nsteps;
    }
    else {
        printf("Error! Unable to get number of HDF5 groups. Returning -1 from h5_get_total_number_of_groups()\n");
        return -1;
    }
}

hid_t H5nb6xx_Helper::h5_open_group_by_name(hid_t parent_id, const char* group_name) {
    hid_t group_id;
    group_id = H5Gopen2(parent_id, group_name, H5P_DEFAULT);
    if (group_id == HDF5_error) { 
        printf("ERROR opening group %s, returning -1 from h5_open_group_by_name()\n", group_name);  
        return -1;
    }
    return group_id;
}


double H5nb6xx_Helper::h5_read_attribute_double(hid_t parent_id, const char *attrib_name) {
    double attrib_val=-1;
    hid_t err;
    hid_t attrib_id;
    attrib_id = H5Aopen (parent_id, attrib_name, H5P_DEFAULT);
    err = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,&attrib_val);
    H5Aclose (attrib_id);
    if (err >= 0) {
        return attrib_val;
    } else {
        printf("Error! Unable to read attribute %s, returning -1 from h5_read_attribute_double()...\n", attrib_name);
        return -1;
    }
}

int H5nb6xx_Helper::h5_read_attribute_integer(hid_t parent_id, const char *attrib_name) {
    int attrib_val=-1;
    hid_t err;
    hid_t attrib_id;
    attrib_id = H5Aopen (parent_id, attrib_name, H5P_DEFAULT);
    err = H5Aread(attrib_id,H5T_NATIVE_INT,&attrib_val);
    H5Aclose (attrib_id);
    if (err >= 0) {
        return attrib_val;
    } else {
        printf("Error! Unable to read attribute %s, returning -1 from h5_read_attribute_double()...\n", attrib_name);
        return -1;
    }
}

int* H5nb6xx_Helper::h5_read_dataset_as_integer_vector(hid_t parent_id, const char *dset_name) {
    hid_t dset_id;
    int * data;
    dset_id = H5Dopen2(parent_id, dset_name, H5P_DEFAULT);
    hid_t dspace_id = H5Dget_space(dset_id);
    int rank = H5Sget_simple_extent_ndims(dspace_id);
    hsize_t dims[rank];
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);

    // allocate the memory
    //data = (int *) malloc((int)dims[0] * sizeof(int));
    data = new int[(int)dims[0]];
    
    // Begin reading
    herr_t err = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    if(err>=0) return data;
    else {
        printf("Error! dataset %s cannot be read! Return NULL from h5_read_dataset_as_integer_vector()...\n", dset_name);
        return NULL;
    }
}

float* H5nb6xx_Helper::h5_read_dataset_as_float_vector(hid_t parent_id, const char *dset_name) {
    hid_t dset_id;
    float * data;
    dset_id = H5Dopen2(parent_id, dset_name, H5P_DEFAULT);
    hid_t dspace_id = H5Dget_space(dset_id);
    int rank = H5Sget_simple_extent_ndims(dspace_id);
    hsize_t dims[rank];
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);

    // allocate the memory
    //data = (float *) malloc((int)dims[0] * sizeof(float));
    data = new float[(int)dims[0]];
    
    // Begin reading
    herr_t err = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    if(err>=0) return data;
    else {
        printf("Error! dataset %s cannot be read! Return NULL from h5_read_dataset_as_float_vector()...\n", dset_name);
        return NULL;
    }

}

double H5nb6xx_Helper::h5_get_step_duration() {
    if(status.h5_file_id<=0) {
        printf("Error! Invalid h5_file_id = %d, returning -1 from h5_get_step_duration()\n", status.h5_file_id);
        return -1;
    }
    int ngroups = h5_get_total_number_of_groups(status.h5_file_id);
    if(ngroups==1) return 0;
    else {
        // open the second group, compare its time attribute with the first one
        hid_t step0 = h5_open_group_by_name(status.h5_file_id, "Step#0");
        hid_t step1 = h5_open_group_by_name(status.h5_file_id, "Step#1");
        double t0 = h5_read_attribute_double(step0, "Time");
        double t1 = h5_read_attribute_double(step1, "Time");
        H5Gclose(step0);
        H5Gclose(step1);
        return t1 - t0;
    }
}


int H5nb6xx_Helper::h5_get_dataset_vector_length(hid_t parent_id, const char *dset_name) {
    hid_t dset_id, dspace_id;
    dset_id = H5Dopen(parent_id, dset_name, H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    const int ndims = H5Sget_simple_extent_ndims(dspace_id);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    return (int)dims[0];
}

int H5nb6xx_Helper::h5_load_step_by_id(int step_id, H5nb6xx_Helper::Dynamics* data) {
    if (step_id < 0 || step_id >= status.nsteps) {
        printf("Error! HDF5 Step#%d cannot be read! Return NULL from h5_load_step_by_id()...\n", step_id); 
        return -1;
    }
    
    // prepare the group name
    char h5_group_name[32];
    sprintf(h5_group_name, "Step#%d", step_id);
    printf("Loading new group %s\n", h5_group_name);
    hid_t gid = h5_open_group_by_name(status.h5_file_id, h5_group_name);
    if (gid < 0){
        printf("Error! HDF5 group %s cannot be opened! Return NULL from h5_load_step_by_id()...\n", h5_group_name);
        return -1;
    }
    
    // When this point has been reached, the group is opened successfully
    data->h5_group_id = gid;

    // Prepare the properties of the group
    data->n_records = h5_get_dataset_vector_length(data->h5_group_id,"X");
    if (data->n_records != status.n_particles) {
        printf("FATAL ERROR: inconsistency number of particles (%d != %d)! Terminating...\n", status.n_particles, data->n_records);
        exit(-1); // terminate the code if the vector of this step is different from the first step
    }
    data->time = h5_read_attribute_double(data->h5_group_id, "Time");
    data->step_id = step_id;
    
    // reading out
    // Since the array index of Fortran starts from 1 instead of 0, here
    // the ID read out from the HDF5 file should be substract by 1 before using!!!
    // Use vec_x, vec_y, vec_z to reorganize the ID 
    int *id_h5_step, *id_sorted; // variable IDs read from the HDF5 step
    id_h5_step = h5_read_dataset_as_integer_vector(data->h5_group_id, "ID");
    id_sorted = new int[data->n_records+1]; // because fortran starts from 1
    if (id_sorted != NULL and id_h5_step != NULL) {
        for (int i = 0; i < data->n_records+1; i++) {
            id_sorted[i] = 0;
        }
        for (int i = 0; i < data->n_records; i++) {
            id_sorted[id_h5_step[i]] = i;
        }
        if (data->id != NULL) {
            delete [] data->id;
            data->id = NULL;
        }
        data->id = id_sorted;
        data->id_original = id_h5_step;
    } else {
        printf("Error! Either id_sorted == NULL or id_h5_step == NULL (load_step_by_id)\n");
    }

    float *vec[13];
    vec[0] = h5_read_dataset_as_float_vector(data->h5_group_id, "Mass");
    vec[1] = h5_read_dataset_as_float_vector(data->h5_group_id, "X");
    vec[2] = h5_read_dataset_as_float_vector(data->h5_group_id, "Y");
    vec[3] = h5_read_dataset_as_float_vector(data->h5_group_id, "Z");
    vec[4] = h5_read_dataset_as_float_vector(data->h5_group_id, "VX");
    vec[5] = h5_read_dataset_as_float_vector(data->h5_group_id, "VY");
    vec[6] = h5_read_dataset_as_float_vector(data->h5_group_id, "VZ");
    vec[7] = h5_read_dataset_as_float_vector(data->h5_group_id, "AX");
    vec[8] = h5_read_dataset_as_float_vector(data->h5_group_id, "AY");
    vec[9] = h5_read_dataset_as_float_vector(data->h5_group_id, "AZ");
    vec[10] = h5_read_dataset_as_float_vector(data->h5_group_id, "JX");
    vec[11] = h5_read_dataset_as_float_vector(data->h5_group_id, "JY");
    vec[12] = h5_read_dataset_as_float_vector(data->h5_group_id, "JZ");

    // allocate memory if necessary 
    if (data->mass==NULL) data->mass = new float[data->n_records];
    if (data->x==NULL) data->x = new float[data->n_records];
    if (data->y==NULL) data->y = new float[data->n_records];
    if (data->z==NULL) data->z = new float[data->n_records];
    if (data->vx==NULL) data->vx = new float[data->n_records];
    if (data->vy==NULL) data->vy = new float[data->n_records];
    if (data->vz==NULL) data->vz = new float[data->n_records];
    if (data->ax==NULL) data->ax = new float[data->n_records];
    if (data->ay==NULL) data->ay = new float[data->n_records];
    if (data->az==NULL) data->az = new float[data->n_records];
    if (data->jx==NULL) data->jx = new float[data->n_records];
    if (data->jy==NULL) data->jy = new float[data->n_records];
    if (data->jz==NULL) data->jz = new float[data->n_records];

    for (int i = 0; i < data->n_records; i++) {
        data->mass[i] = vec[0][id_sorted[i+1]];
        data->x[i] = vec[1][id_sorted[i+1]];
        data->y[i] = vec[2][id_sorted[i+1]];
        data->z[i] = vec[3][id_sorted[i+1]];
        data->vx[i] = vec[4][id_sorted[i+1]];
        data->vy[i] = vec[5][id_sorted[i+1]];
        data->vz[i] = vec[6][id_sorted[i+1]];
        data->ax[i] = vec[7][id_sorted[i+1]];
        data->ay[i] = vec[8][id_sorted[i+1]];
        data->az[i] = vec[9][id_sorted[i+1]];
        data->jx[i] = vec[10][id_sorted[i+1]];
        data->jy[i] = vec[11][id_sorted[i+1]];
        data->jz[i] = vec[12][id_sorted[i+1]];
    }
    for (int i = 0; i < 13; i++) {
        if (vec[i] != NULL) {
            delete [] vec[i];
            vec[i] = NULL;
        }
    }

    return 0;
}



int H5nb6xx_Helper::helper_get_h5_filename(char **h5_filename){
    *h5_filename = status.h5_filename;
    return 0;
}

int H5nb6xx_Helper::helper_set_h5_filename(char *h5_filename){
    strcpy(status.h5_filename, h5_filename);
    return 0;
}

int H5nb6xx_Helper::helper_get_eps2(double * epsilon_squared){
    *epsilon_squared = 0;
    return 0;
}

int H5nb6xx_Helper::helper_set_eps2(double epsilon_squared){
    return 0;
}

int H5nb6xx_Helper::helper_initialize_code(){
    //H5nb6xx_Helper::instance = this;
    // Try to open the h5 file. If not specified, the default h5 file is 'data.h5part'.
    printf("Code initialization begins...\n");
    char h5_fn_current[256];
    hid_t       file_id, h5_group_id;
    if(strlen(status.h5_filename)>0) {
        strcpy(h5_fn_current, status.h5_filename);
    } else {
        strcpy(h5_fn_current, "data.h5part");
    }
    file_id = H5Fopen(h5_fn_current, H5F_ACC_RDWR, H5P_DEFAULT);
    status.current_step_id = 0;
    status.prev_step_id = 0;
    status.next_step_id = 1;
    status.current_time = 0.0;
    if (file_id < 0) {
        printf("Error! Invalid HDF5 file ID %d, returning -1 from helper_initialize_code()\n", file_id);
        return -1; // Error opening HDF5 file
    } 
    // Open the first h5 group to determine some global properties
    char h5_group_name[32];
    status.h5_file_id = file_id;
    sprintf(h5_group_name, "Step#%d", status.current_step_id);
    h5_group_id = h5_open_group_by_name(status.h5_file_id, h5_group_name);
    
    // Determine the time range covered by the HDF5 snapshot
    status.begin_time = h5_read_attribute_double(h5_group_id, "Time");
    //status.current_step_time = h5_read_attribute_double(status.h5_group_id, "Time");
    status.t_step = h5_get_step_duration();
    status.end_time = status.t_step * status.nsteps; 
    printf("Snapshot covered from %f to %f\n", status.begin_time, status.end_time);

    // determine the number of particle and records in Step#0
    status.n_particles = h5_read_attribute_integer(h5_group_id, "TotalN");
    int n_records = h5_get_dataset_vector_length(h5_group_id, "X");
    
    // upscale the MAX_PARTICLE_NUMBER to be one magnitude larger
    if(status.n_particles>0) {
        MAX_PARTICLE_NUMBER = status.n_particles+100;
    } else if (n_records>0) { // TotalN attribute not available
        MAX_PARTICLE_NUMBER = (int) pow(10,ceil(log10(n_records)));
    }
    printf("Total number of particles: %d, MAX_PARTICLE_NUMBER = %d\n", status.n_particles, MAX_PARTICLE_NUMBER);
    H5Gclose(h5_group_id);
    status.current_step_id = -1;
    status.data = new H5nb6xx_Helper::Dynamics();
    status.next_data = new H5nb6xx_Helper::Dynamics();
    status.interp_data = new H5nb6xx_Helper::Dynamics();
    host_star_flag = new int[MAX_PARTICLE_NUMBER];
    for (int i = 0; i < MAX_PARTICLE_NUMBER; i++) {
        host_star_flag[i] = 0;
    }
    helper_evolve_model(0);
    return 0;
}

int H5nb6xx_Helper::helper_new_particle(int * index_of_the_particle, double mass, 
  double x, double y, double z, double vx, double vy, double vz, 
  double radius) {
    int i = next_particle_id;
    status.data->x[i] = x;
    status.data->y[i] = y;
    status.data->z[i] = z;
    status.data->vx[i] = vx;
    status.data->vy[i] = vy;
    status.data->vz[i] = vz;
    status.data->ax[i] = 0;
    status.data->ay[i] = 0;
    status.data->az[i] = 0;
    status.data->jx[i] = 0;
    status.data->jy[i] = 0;
    status.data->jz[i] = 0;
    status.data->mass[i] =mass;
    *index_of_the_particle = next_particle_id;
    next_particle_id++;
    return 0;
}


int H5nb6xx_Helper::helper_delete_particle(int index_of_the_particle) {
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    return -1; // Particle cannot be removed
}

int H5nb6xx_Helper::helper_get_number_of_particles(int * number_of_particles){
    *number_of_particles = status.n_particles;
    return 0;
}

int H5nb6xx_Helper::helper_get_index_of_first_particle(int * index_of_the_particle){
    *index_of_the_particle = 0;
    return 0;
}


int H5nb6xx_Helper::helper_get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle) {
    *index_of_the_next_particle = index_of_the_particle + 1;
    printf("get_index_of_next_particle called.\n");
    return 0;
}

int H5nb6xx_Helper::helper_get_state(int index_of_the_particle, double * mass, double * x, 
  double * y, double * z, double * vx, double * vy, double * vz, double * radius){
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    *mass = status.data->mass[index_of_the_particle];
    *x = status.data->x[index_of_the_particle];
    *y = status.data->y[index_of_the_particle];
    *z = status.data->z[index_of_the_particle];
    *vx = status.data->vx[index_of_the_particle];
    *vy = status.data->vy[index_of_the_particle];
    *vz = status.data->vz[index_of_the_particle];
    *radius = 0; // Radius attribute not supported yet
    return 0;
}


int H5nb6xx_Helper::helper_set_state(int index_of_the_particle, double mass, double x, 
        double y, double z, double vx, double vy, double vz, double radius){
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    status.data->x[index_of_the_particle] = x;
    status.data->y[index_of_the_particle] = y;
    status.data->z[index_of_the_particle] = z;
    status.data->vx[index_of_the_particle] = vx;
    status.data->vy[index_of_the_particle] = vy;
    status.data->vz[index_of_the_particle] = vz;
    //data.radius[index_of_the_particle] = radius;
    return 0; //modification of the HDF5 file not supported
}

int H5nb6xx_Helper::helper_get_mass(int index_of_the_particle, double * mass) {
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    *mass = status.data->mass[index_of_the_particle];
    return 0;
}

int H5nb6xx_Helper::helper_set_mass(int index_of_the_particle, double mass) {
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    status.data->mass[index_of_the_particle] = mass;
    return 0; // not supported
}

int H5nb6xx_Helper::helper_get_position(int index_of_the_particle, double * x, double * y, double * z) {
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    *x = status.data->x[index_of_the_particle];
    *y = status.data->y[index_of_the_particle];
    *z = status.data->z[index_of_the_particle];
    return 0;

}

int H5nb6xx_Helper::helper_set_position(int index_of_the_particle, double x, double y, double z) {
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    status.data->x[index_of_the_particle] = x;
    status.data->y[index_of_the_particle] = y;
    status.data->z[index_of_the_particle] = z;
    return 0; // not supported
}

int H5nb6xx_Helper::helper_set_acceleration(int index_of_the_particle, double ax, double ay, double az) {
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    status.data->ax[index_of_the_particle] = ax;
    status.data->ay[index_of_the_particle] = ay;
    status.data->az[index_of_the_particle] = az;
    return 0; //not supported
}


int H5nb6xx_Helper::helper_get_acceleration(int index_of_the_particle, double * ax, double * ay, double * az) {
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    *ax = status.data->ax[index_of_the_particle];
    *ay = status.data->ay[index_of_the_particle];
    *az = status.data->az[index_of_the_particle];
    return 0;
}

int H5nb6xx_Helper::helper_get_potential(int index_of_the_particle, double * potential) {
    index_of_the_particle -= 1; // the ID is already sorted, starting from 0
    return 0;
}

int H5nb6xx_Helper::helper_commit_particles() {
    return 0;
}

int H5nb6xx_Helper::helper_get_time(double * time) {
    if(status.data->h5_group_id == 0) return -1;
    *time = h5_read_attribute_double(status.data->h5_group_id, "Time");
    return 0;
}


int H5nb6xx_Helper::helper_get_kinetic_energy(double * kinetic_energy) {

    return 0;
}

int H5nb6xx_Helper::helper_get_potential_energy(double * potential_energy) {

    return 0;
}


int H5nb6xx_Helper::helper_get_center_of_mass_velocity(double * vx, double * vy, double * vz) {
    double mtot = 0;
    double cm_vx=0, cm_vy=0, cm_vz=0;
    for(int i=0;i<status.n_particles;i++) {
        mtot += status.data->mass[i];
        cm_vx += status.data->mass[i] * status.data->vx[i];
        cm_vy += status.data->mass[i] * status.data->vy[i];
        cm_vz += status.data->mass[i] * status.data->vz[i];
    }
    *vx = cm_vx/mtot;
    *vy = cm_vy/mtot;
    *vz = cm_vz/mtot;
    return 0;
}

int H5nb6xx_Helper::helper_get_center_of_mass_position(double * x, double * y, double * z) {
    double mtot = 0;
    double cm_x=0, cm_y=0, cm_z=0;
    for(int i=0;i<status.n_particles;i++) {
        mtot += status.data->mass[i];
        cm_x += status.data->mass[i] * status.data->x[i];
        cm_y += status.data->mass[i] * status.data->y[i];
        cm_z += status.data->mass[i] * status.data->z[i];
    }
    *x = cm_x/mtot;
    *y = cm_y/mtot;
    *z = cm_z/mtot;
    return 0;
}


int H5nb6xx_Helper::helper_get_total_mass(double * mass) {
    double mtot = 0;
    for(int i=0;i<status.n_particles;i++) {
        mtot += status.data->mass[i];
    }
    *mass = mtot;
    return 0;
}

int H5nb6xx_Helper::helper_get_total_radius(double * radius) {
    *radius = 0;
    return 0; // not supported yet
}

int H5nb6xx_Helper::helper_get_radius(int index_of_the_particle, double * radius) {
    index_of_the_particle = status.data->id[index_of_the_particle];
    *radius = 0;
    return 0;
}


int H5nb6xx_Helper::helper_set_radius(int index_of_the_particle, double radius) {
    index_of_the_particle = status.data->id[index_of_the_particle];
    //data.radius[index_of_the_particle] = radius;
    return 0; // not supported
}

int H5nb6xx_Helper::helper_cleanup_code() {
    return 0;
}

int H5nb6xx_Helper::helper_evolve_model(double to_time) {
    // On return, system_time will be greater than or equal to the specified time. 
    bool skip_loading = false;

    if (to_time < 0 || to_time > status.end_time) {
        printf("Error! invalid time (%f), returning -1 from helper_evolve_model()\n", to_time);
        return -1;
    }

    // Advance the system time to the specified time
    status.current_time = to_time;
    
    // Locate the group corresponding to the time specified, 
    // so that time_of_the_group <= to_time < time_of_next_group
    int step = floor(to_time / status.t_step);
    //int step = round(to_time / status.t_step);
    printf("current: %d, step: %d, next: %d\n", status.current_step_id, step, status.next_step_id);

    // if the next_data loaded previously happens to be the data needed for this time,
    // then use the next_data to replace the current data, and then load new next_data
    // with respect to the current step. Otherwise, both current_data and next_data have
    // to be loaded (usually happens at t = 0 or random seeking).

    /*
    if (step != status.current_step_id || NULL == status.data || NULL == status.next_data) { // loading needed
        status.prev_step_id = status.current_step_id;
        status.current_step_id = step;
        status.next_step_id = status.current_step_id + 1;
        H5nb6xx_Helper::Dynamics * tmpdata;

        if (NULL != status.next_data) {
            if (step == status.next_data->step_id){
                tmpdata = status.data;
                status.data = status.next_data;
                status.next_data = h5_load_step_by_id(status.next_step_id);
                if (status.next_data != NULL) {
                    delete tmpdata;
                } else {
                    status.next_data = tmpdata; // restore the current data as the next step data
                }
            }
        } else {
            status.data = h5_load_step_by_id(status.current_step_id);
            status.next_data = h5_load_step_by_id(status.next_step_id);
        }
    }
    */
    

    //if (step != status.current_step_id) { 
        status.prev_step_id = status.current_step_id;
        status.current_step_id = step;
        status.next_step_id = status.current_step_id + 1;
        h5_load_step_by_id(status.current_step_id, status.data);
        h5_load_step_by_id(status.next_step_id, status.next_data);
    //}
    

#ifdef GPU
    if(enable_interpolation) {
        this->cuda_util->cuda_predict(to_time);
    }
#endif

    return 0;
}

int H5nb6xx_Helper::helper_get_velocity(int index_of_the_particle, double * vx, double * vy, double * vz) {
    index_of_the_particle = status.data->id[index_of_the_particle];
    *vx = status.data->vx[index_of_the_particle];
    *vy = status.data->vy[index_of_the_particle];
    *vz = status.data->vz[index_of_the_particle];
    return 0;
}


int H5nb6xx_Helper::helper_set_velocity(int index_of_the_particle, double vx, double vy, double vz) {
    index_of_the_particle = status.data->id[index_of_the_particle];
    status.data->vx[index_of_the_particle] = vx;
    status.data->vy[index_of_the_particle] = vy;
    status.data->vz[index_of_the_particle] = vz;
    return 0; // not supported yet
}

int H5nb6xx_Helper::helper_get_time_step(double * time_step) {
    *time_step = status.t_step;
    return 0;
}


int H5nb6xx_Helper::helper_get_begin_time(double * output) { 
    *output = status.begin_time;
    return 0;
}

int H5nb6xx_Helper::helper_set_begin_time(double input) {

    return 0;
}

int H5nb6xx_Helper::helper_commit_parameters()
{
    // Perform any needed setup after initial code parameters have been set.

    // Consistency check:

    return 0;
}

int H5nb6xx_Helper::helper_synchronize_model()
{
    // Synchronize all particles at the current system time.  The
    // default is not to reinitialize the scheduler, as this will be
    // handled later, in recommit_particles().
    return 0;
}

int H5nb6xx_Helper::helper_recommit_particles()
{
    // Reinitialize/reset the system after particles have been added
    // or removed.  The system should be synchronized at some reasonable
    // system_time, so we just need to recompute forces and update the
    // GPU and scheduler.  Note that we don't resize the jdata or
    // idata arrays.  To resize idata, just delete and create a new
    // one.  Resizing jdata is more complicated -- defer for now.

    return 0;
}


int H5nb6xx_Helper::helper_recommit_parameters()
{
    // Perform any needed changes after code parameters have been reset.

    return 0;

}


int H5nb6xx_Helper::helper_get_potential_at_point(double * eps, double * x, double * y, double * z, double * phi, int n) {
    // Inquirying the potentials of n points specified by (x[0], y[0], z[0]), 
    // (x[1], y[1], z[1]), ..., (x[n-1, y[n-1], z[n-1])
/*
    if(n<1) return -1; // at least one point 
    phi = (double *) malloc(n * sizeof(double));
    for(int i=0; i<n; i++) {
        phi[i] = 0;
        double r2 = 0; // distance squared
        double r2i = 0;
        double ri = 0;
        for(int j=0; j<status.n_particles; j++) {
            // sum all particles in the cluster
            r2 = pow(data.x[j]-x[i], 2) + pow(data.y[j]-y[i], 2) + pow(data.z[j]-z[i], 2);
            r2i = 1.0/(r2 + eps[i] + _TINY_);
            ri = sqrt(r2i);
            phi[i] -= data.mass[j] * ri;
        }
        printf("potential ph[%d] = %f\n", i, phi[i]);
    }
*/
        *phi = 0;
        double r2 = 0; // distance squared
        double r2i = 0;
        double ri = 0;
        for(int j=0; j<status.data->n_records; j++) {
            // sum all particles in the cluster
            r2 = pow(status.data->x[j]-(*x), 2) + pow(status.data->y[j]-(*y), 2) + pow(status.data->z[j]-(*z), 2);
            r2i = 1.0/(r2 + *eps + _TINY_);
            ri = sqrt(r2i);
            *phi -= status.data->mass[j] * ri;
        }
    return 0;
}


int H5nb6xx_Helper::helper_get_gravity_at_point(double * eps, double * x, double * y, double * z,
             double * forcex, double * forcey, double * forcez, int n) {
    // Inquirying the accelerations of n points specified by (x[0], y[0], z[0]), 
    // (x[1], y[1], z[1]), ..., (x[n-1, y[n-1], z[n-1])
    if (n > 1) {
        forcex = (double *) malloc(n * sizeof(double));
        forcey = (double *) malloc(n * sizeof(double));
        forcez = (double *) malloc(n * sizeof(double));
        for(int i=0; i<n; i++) {
            forcex[i] = 0;
            forcey[i] = 0;
            forcez[i] = 0;
            double r2 = 0; // distance squared
            double r2i = 0;
            double ri = 0;
            double mri = 0;
            double mr3i = 0;
            for(int j=0; j<status.n_particles; j++) {
                // sum all particles in the cluster, but skip the host stars
                if (host_star_flag[j]>0) {
                    //printf("Skipping particle #%d because it is the host star\n", j+1);
                    continue; // already sorted, use j instead of (j+1)
                }
                r2 = pow(status.data->x[j]-x[i], 2.0) + pow(status.data->y[j]-y[i], 2.0) + pow(status.data->z[j]-z[i], 2.0);
                if (r2 < _TINY_) {
                    //printf("WARNING!! SMALL DISTANCE SKIPPED! j=%d, x[j]=%f, *x=%f, r2=%f\n", j, status.data->x[j], (*x), r2);
                    continue;
                }
                //r2i = 1.0/(r2 + (*eps) + _TINY_);
                r2i = 1.0/(r2 + (*eps));
                ri = sqrt(r2i);
                mri = status.data->mass[j] * ri;
                mr3i = mri * r2i;
                forcex[i] += mr3i * (status.data->x[j]-x[i]);
                forcey[i] += mr3i * (status.data->y[j]-y[i]);
                forcez[i] += mr3i * (status.data->z[j]-z[i]);
            }
        }
    } else if (n == 1) {
        *forcex = 0;
        *forcey = 0;
        *forcez = 0;
        double r2 = 0; // distance squared
        double r2i = 0;
        double ri = 0;
        double mri = 0;
        double mr3i = 0;
        for(int j=0; j<status.n_particles; j++) {
            // sum all particles in the cluster, but skip the host stars
            if (host_star_flag[j]>0) {
                //printf("Skipping particle #%d because it is the host star\n", j+1);
                continue; // already sorted, use j instead of (j+1)
            }
            r2 = pow(status.data->x[j]-(*x), 2.0) + pow(status.data->y[j]-(*y), 2.0) + pow(status.data->z[j]-(*z), 2.0);
            if (r2 < _TINY_) {
                //printf("WARNING!! SMALL DISTANCE SKIPPED! j=%d, x[j]=%f, *x=%f, r2=%f\n", j, status.data->x[j], (*x), r2);
                continue;
            }
            //r2i = 1.0/(r2 + (*eps) + _TINY_);
            r2i = 1.0/(r2 + (*eps));
            ri = sqrt(r2i);
            mri = status.data->mass[j] * ri;
            mr3i = mri * r2i;
            *forcex += mr3i * (status.data->x[j]-(*x));
            *forcey += mr3i * (status.data->y[j]-(*y));
            *forcez += mr3i * (status.data->z[j]-(*z));
        }
    }
    return 0;

}

int H5nb6xx_Helper::helper_get_total_number_of_steps(int *val) {
    *val = status.nsteps;
    return 0;
}

int H5nb6xx_Helper::helper_set_total_number_of_steps(int val){
    status.nsteps = val;
    return 0;
}


int H5nb6xx_Helper::helper_get_duration_per_step(double *val){
    *val = status.t_step;
    return 0;
}

int H5nb6xx_Helper::helper_set_duration_per_step(double val) {
    status.t_step = val;
    return 0;
}

int H5nb6xx_Helper::helper_set_eta(double timestep_parameter)
{
    return 0;
}

int H5nb6xx_Helper::helper_get_eta(double * timestep_parameter)
{
    *timestep_parameter = 0.05;
    return 0;
}

int H5nb6xx_Helper::helper_set_enable_interpolation(bool val)
{
    enable_interpolation = val;
    return 0;
}

int H5nb6xx_Helper::helper_set_host_star_flag(int host_star_id, int flag)
{
    host_star_id -= 1; // the ID is already sorted, starting from 0
    printf("setting %d as %d\n", host_star_id, flag);
    host_star_flag[host_star_id] = flag;
    return 0;
}

int H5nb6xx_Helper::helper_get_neighbors(int *host_star_id, int* neighbor_star_id, int n) {
    double d, d_min;
    double dx, dy, dz;
    int hsid = 0;
    //neighbor_star_id = new int[n];
    for(int i = 0; i < n; i++){
        d_min = DBL_MAX;
        neighbor_star_id[i] = -1;
        for(int j = 0; j < status.n_particles; j++) {
            hsid = status.data->id[host_star_id[i]];
            if(hsid == j) continue;
            dx = status.data->x[hsid] - status.data->x[j];
            dy = status.data->y[hsid] - status.data->y[j];
            dz = status.data->z[hsid] - status.data->z[j];
            d = sqrt(dx * dx + dy * dy + dz * dz);
            if (d < d_min) {
                d_min = d;
                neighbor_star_id[i] = status.data->id_original[j];
            }
        }
        printf("neighbor for %d is %d, rmin=%f\n", host_star_id[i], neighbor_star_id[i], d_min);
    }
    return 0;
}
