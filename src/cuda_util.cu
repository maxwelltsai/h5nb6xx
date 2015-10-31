#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>

#include "cuda_util.h"
#include "h5nb6xx_helper.h"

struct interpolate
{
    __host__ __device__
    float3 operator()(const thrust::tuple<float3, float3, float3, float3, float3, float3, float3, float3, float3>& vec)
    {
        float3 X = thrust::get<0>(vec);
        float3 V = thrust::get<1>(vec);
        float3 A = thrust::get<2>(vec);
        float3 J = thrust::get<3>(vec);
        float3 X1 = thrust::get<4>(vec);
        float3 V1 = thrust::get<5>(vec);
        float3 A1 = thrust::get<6>(vec);
        float3 J1 = thrust::get<7>(vec);
        float3 T  = thrust::get<8>(vec);

        float t = T.x;
        float t0 = T.y;
        float tstep = T.z;
        float dt = tstep;
        if (t<t0) t = t0; // fix a bug of the misalignment of HDF5 steps
        float tau = (t-t0)/tstep;
        float p0x = X.x;
        float p0y = X.y;
        float p0z = X.z;
        float p1x = V.x*dt;
        float p1y = V.y*dt;
        float p1z = V.z*dt;
        float p2x = 1.0/2*A.x*dt*dt;
        float p2y = 1.0/2*A.y*dt*dt;
        float p2z = 1.0/2*A.z*dt*dt;
        float p3x = 1.0/6*J.x*dt*dt*dt;
        float p3y = 1.0/6*J.y*dt*dt*dt;
        float p3z = 1.0/6*J.z*dt*dt*dt;
        float p4x = -1.0/6*(4*J.x+J1.x)*dt*dt*dt - 2.5*(2*A.x-A1.x)*dt*dt - 5*(4*V.x+3*V1.x)*dt - 35*(X.x-X1.x);
        float p4y = -1.0/6*(4*J.y+J1.y)*dt*dt*dt - 2.5*(2*A.y-A1.y)*dt*dt - 5*(4*V.y+3*V1.y)*dt - 35*(X.y-X1.y);
        float p4z = -1.0/6*(4*J.z+J1.z)*dt*dt*dt - 2.5*(2*A.z-A1.z)*dt*dt - 5*(4*V.z+3*V1.z)*dt - 35*(X.z-X1.z);
        float p5x = 0.5*(2*J.x+J1.x)*dt*dt*dt + (10*A.x-7*A1.x)*dt*dt + 3*(15*V.x+13*V1.x)*dt + 84*(X.x-X1.x);
        float p5y = 0.5*(2*J.y+J1.y)*dt*dt*dt + (10*A.y-7*A1.y)*dt*dt + 3*(15*V.y+13*V1.y)*dt + 84*(X.y-X1.y);
        float p5z = 0.5*(2*J.z+J1.z)*dt*dt*dt + (10*A.z-7*A1.z)*dt*dt + 3*(15*V.z+13*V1.z)*dt + 84*(X.z-X1.z);
        float p6x = -1.0/6*(4*J.x+3*J1.x)*dt*dt*dt - 0.5*(15*A.x-13*A1.x)*dt*dt - 2*(18*V.x+17*V1.x)*dt - 70*(X.x-X1.x);
        float p6y = -1.0/6*(4*J.y+3*J1.y)*dt*dt*dt - 0.5*(15*A.y-13*A1.y)*dt*dt - 2*(18*V.y+17*V1.y)*dt - 70*(X.y-X1.y);
        float p6z = -1.0/6*(4*J.z+3*J1.z)*dt*dt*dt - 0.5*(15*A.z-13*A1.z)*dt*dt - 2*(18*V.z+17*V1.z)*dt - 70*(X.z-X1.z);
        float p7x = 1.0/6*(J.x+J1.x)*dt*dt*dt + 2*(A.x-A1.x)*dt*dt + 10*(V.x+V1.x)*dt + 20*(X.x-X1.x);
        float p7y = 1.0/6*(J.y+J1.y)*dt*dt*dt + 2*(A.y-A1.y)*dt*dt + 10*(V.y+V1.y)*dt + 20*(X.y-X1.y);
        float p7z = 1.0/6*(J.z+J1.z)*dt*dt*dt + 2*(A.z-A1.z)*dt*dt + 10*(V.z+V1.z)*dt + 20*(X.z-X1.z);

        float x_pred = p0x + p1x*tau + p2x*pow(tau,2.0f) + p3x*pow(tau,3.0f) + p4x*pow(tau,4.0f) + p5x*pow(tau,5.0f) + p6x*pow(tau,6.0f) + p7x*pow(tau,7.0f);
        float y_pred = p0y + p1y*tau + p2y*pow(tau,2.0f) + p3y*pow(tau,3.0f) + p4y*pow(tau,4.0f) + p5y*pow(tau,5.0f) + p6y*pow(tau,6.0f) + p7y*pow(tau,7.0f);
        float z_pred = p0z + p1z*tau + p2z*pow(tau,2.0f) + p3z*pow(tau,3.0f) + p4z*pow(tau,4.0f) + p5z*pow(tau,5.0f) + p6z*pow(tau,6.0f) + p7z*pow(tau,7.0f);

        // fail-safe linear interpolation based on positions in case the v, a, j vectors are wrong
        // assuming that the interpolated x should be x0 <= x_interp <= x1
        //if ((x_pred<X.x||x_pred>X1.x) || (y_pred<X.y||y_pred>X1.y) || (z_pred<X.z||z_pred>X1.z)) {
        //    printf("crazy!!, x0=%f, x1=%f, x_pred=%f, t0=%f, t=%f, tau=%f\n", X.x, X1.x, x_pred, t0, t, tau);
        //    x_pred = X.x + tau * (X1.x - X.x);
        //    y_pred = X.y + tau * (X1.y - X.y);
        //    z_pred = X.z + tau * (X1.z - X.z);
        //}
        return make_float3(x_pred, y_pred, z_pred);
    }
};

struct tuple_to_float3 {
    __host__ __device__
    float3 operator()(thrust::tuple<float, float, float> vec) {
        float x = thrust::get<0>(vec);
        float y = thrust::get<1>(vec);
        float z = thrust::get<2>(vec);
        return make_float3(x, y, z);
    }
};

CUDA_Util::CUDA_Util(H5nb6xx_Helper* helper){
    this->h5nb6xx_helper = helper;
}
/*
int CUDA_Util::cuda_load_data() {

    return 0;
}

int CUDA_Util::cuda_free_data() {
    return 0;
}
*/

int CUDA_Util::cuda_predict(float to_time){
    H5nb6xx_Helper::Status istatus = this->h5nb6xx_helper->get_status();
    H5nb6xx_Helper::Dynamics* idata = this->h5nb6xx_helper->get_data();
    H5nb6xx_Helper::Dynamics* idata1 = this->h5nb6xx_helper->get_data_next();

    int n_particles = istatus.n_particles;
    float current_time = idata->time;
    std::cout<<"interpolation, dt="<<(to_time - current_time)<<" t0="<<current_time<<" t="<<to_time<<" t1="<<idata1->time<<std::endl;
    thrust::device_vector<float> x(idata->x, idata->x + n_particles);
    thrust::device_vector<float> y(idata->y, idata->y + n_particles);
    thrust::device_vector<float> z(idata->z, idata->z + n_particles);
    thrust::device_vector<float> vx(idata->vx, idata->vx + n_particles);
    thrust::device_vector<float> vy(idata->vy, idata->vy + n_particles);
    thrust::device_vector<float> vz(idata->vz, idata->vz + n_particles);
    thrust::device_vector<float> ax(idata->ax, idata->ax + n_particles);
    thrust::device_vector<float> ay(idata->ay, idata->ay + n_particles);
    thrust::device_vector<float> az(idata->az, idata->az + n_particles);
    thrust::device_vector<float> jx(idata->jx, idata->jx + n_particles);
    thrust::device_vector<float> jy(idata->jy, idata->jy + n_particles);
    thrust::device_vector<float> jz(idata->jz, idata->jz + n_particles);

    thrust::device_vector<float> x1(idata1->x, idata1->x + n_particles);
    thrust::device_vector<float> y1(idata1->y, idata1->y + n_particles);
    thrust::device_vector<float> z1(idata1->z, idata1->z + n_particles);
    thrust::device_vector<float> vx1(idata1->vx, idata1->vx + n_particles);
    thrust::device_vector<float> vy1(idata1->vy, idata1->vy + n_particles);
    thrust::device_vector<float> vz1(idata1->vz, idata1->vz + n_particles);
    thrust::device_vector<float> ax1(idata1->ax, idata1->ax + n_particles);
    thrust::device_vector<float> ay1(idata1->ay, idata1->ay + n_particles);
    thrust::device_vector<float> az1(idata1->az, idata1->az + n_particles);
    thrust::device_vector<float> jx1(idata1->jx, idata1->jx + n_particles);
    thrust::device_vector<float> jy1(idata1->jy, idata1->jy + n_particles);
    thrust::device_vector<float> jz1(idata1->jz, idata1->jz + n_particles);

    thrust::device_vector<float3> X(n_particles);
    thrust::device_vector<float3> V(n_particles);
    thrust::device_vector<float3> A(n_particles);
    thrust::device_vector<float3> J(n_particles);
    thrust::device_vector<float3> X1(n_particles);
    thrust::device_vector<float3> V1(n_particles);
    thrust::device_vector<float3> A1(n_particles);
    thrust::device_vector<float3> J1(n_particles);
    thrust::device_vector<float3> T(n_particles);
    thrust::device_vector<float> t(n_particles);
    thrust::device_vector<float> t0(n_particles);
    thrust::device_vector<float> tstep(n_particles);

    thrust::host_vector<float3> X_h(n_particles);

    thrust::fill(t.begin(), t.end(), to_time);
    thrust::fill(t0.begin(), t0.end(), current_time);
    thrust::fill(tstep.begin(), tstep.end(), istatus.t_step);

    std::cout<<"t="<<to_time<<" t0="<<current_time<<" tstep="<<istatus.t_step<<" tau="<<((to_time-current_time)/istatus.t_step)<<std::endl;

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(x.begin(), y.begin(), z.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(x.end(), y.end(), z.end())),
            X.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(vx.begin(), vy.begin(), vz.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(vx.end(), vy.end(), vz.end())),
            V.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(ax.begin(), ay.begin(), az.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(ax.end(), ay.end(), az.end())),
            A.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(jx.begin(), jy.begin(), jz.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(jx.end(), jy.end(), jz.end())),
            J.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(t.begin(), t0.begin(), tstep.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(t.end(), t0.end(), tstep.end())),
            T.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(x1.begin(), y1.begin(), z1.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(x1.end(), y1.end(), z1.end())),
            X1.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(vx1.begin(), vy1.begin(), vz1.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(vx1.end(), vy1.end(), vz1.end())),
            V1.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(ax1.begin(), ay1.begin(), az1.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(ax1.end(), ay1.end(), az1.end())),
            A1.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(jx1.begin(), jy1.begin(), jz1.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(jx1.end(), jy1.end(), jz1.end())),
            J1.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(t.begin(), t0.begin(), tstep.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(t.end(), t0.end(), tstep.end())),
            T.begin(), 
            tuple_to_float3());

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(X.begin(), V.begin(), A.begin(), J.begin(), X1.begin(), V1.begin(), A1.begin(), J1.begin(), T.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(X.end(), V.end(), A.end(), J.end(), X1.end(), V1.end(), A1.end(), J1.end(), T.end())),
            X.begin(),
            interpolate());

    thrust::copy(X.begin(), X.end(), X_h.begin());


    for(int i=0; i<n_particles;i++) {
        idata->x[i] = X_h[i].x;
        idata->y[i] = X_h[i].y;
        idata->z[i] = X_h[i].z;
    }

    // clean up memory
    x.clear();
    y.clear();
    z.clear();
    vx.clear();
    vy.clear();
    vz.clear();
    ax.clear();
    ay.clear();
    az.clear();
    jx.clear();
    jy.clear();
    jz.clear();
    x1.clear();
    y1.clear();
    z1.clear();
    vx1.clear();
    vy1.clear();
    vz1.clear();
    ax1.clear();
    ay1.clear();
    az1.clear();
    jx1.clear();
    jy1.clear();
    jz1.clear();
    X.clear();
    V.clear();
    A.clear();
    J.clear();
    X1.clear();
    V1.clear();
    A1.clear();
    J1.clear();
    T.clear();
    t.clear();
    t0.clear();
    tstep.clear();
    X_h.clear();

    //delete istatus;
    //delete idata;
    //delete idata1;
    return 0;
}


int CUDA_Util::cuda_get_acceleration(int* x, int* y, int* z, int n_points, float time) {

    return 0;
}

int CUDA_Util::cuda_get_potential(int* x, int* y, int* z, int n_points, float time) {

    return 0;
}

