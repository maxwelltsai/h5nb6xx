#ifndef _CUDA_UTIL_H
#define _CUDA_UTIL_H

#include "h5nb6xx_helper.h"
#include <thrust/device_vector.h>
class H5nb6xx_Helper;

class CUDA_Util
{
    protected:
        H5nb6xx_Helper* h5nb6xx_helper;
        thrust::device_vector<float3> pos_interp;
        thrust::device_vector<float> dist_squared;
        thrust::device_vector<float> mass;
        thrust::device_vector<float> acc;

    public:
        CUDA_Util(H5nb6xx_Helper*);
        int cuda_predict(float);
        int cuda_get_nth_neighbor(int particle_id, int nth_neighbor);
        int* cuda_get_neighbors(int particle_id, int n_neighbors);
        int cuda_get_acceleration(int*, int*, int*, int, float);
        int cuda_get_potential(int*, int*, int*, int, float);
};

#endif
