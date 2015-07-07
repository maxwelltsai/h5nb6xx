#ifndef _CUDA_UTIL_H
#define _CUDA_UTIL_H

#include "h5nb6xx_helper.h"
class H5nb6xx_Helper;

class CUDA_Util
{
    protected:
        H5nb6xx_Helper* h5nb6xx_helper;

    public:
        CUDA_Util(H5nb6xx_Helper*);
        int cuda_predict(float);
        int cuda_get_acceleration(int*, int*, int*, int, float);
        int cuda_get_potential(int*, int*, int*, int, float);

};

#endif
