#pragma once

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

#include "functions.hpp"
#include "Kernels_CUDA_Cpp.hpp"

//? ----------------------------------------------------------
//?
//? Description:
//?     A collection of explicit integrators
//?
//? ----------------------------------------------------------

template <typename rhs>
void explicit_Euler(rhs RHS, double* u, double* u_sol, double* u_temp, double dt, size_t N, bool GPU)
{
    //? RHS(u) = u_sol = du/dt
    RHS(u, u_sol);

    //? u_sol = u + du/dt * dt
    LeXInt::axpby(1.0, u, dt, u_sol, u_sol, N, GPU);
}


template <typename rhs>
void RK2(rhs RHS, double* u, double* u_sol, double* u_temp, double dt, size_t N, bool GPU)
{
    //? Assign names and variables
    double* k1 = &u_temp[0]; double* k2 = &u_temp[N];

    //? Internal stage 1: k1 = dt * RHS(u)
    RHS(u, k1);                                     //* k1 = RHS(u) = du/dt
    LeXInt::axpby(dt, k1, k1, N, GPU);              //* k1 = k1 * dt

    //? Internal stage 2: k2 = dt * RHS(u + k1)
    LeXInt::axpby(1.0, u, 1.0, k1, u_sol, N, GPU);  //* u_sol = u + k1
    RHS(u_sol, k2);                                 //* k2 = RHS(u + k1) = RHS(u_sol) 
    LeXInt::axpby(dt, k2, k2, N, GPU);              //* k2 = k2 * dt

    //? u_rk2 = u + 1./2.*(k1 + k2)
    LeXInt::axpby(1.0, u, 0.5, k1, 0.5, k2, u_sol, N, GPU);
}


template <typename rhs>
void RK4(rhs RHS, double* u, double* u_sol, double* u_temp, double dt, size_t N, bool GPU)
{
    //? Assign names and variables
    double* k1 = &u_temp[0];   double* k2 = &u_temp[N];
    double* k3 = &u_temp[2*N]; double* k4 = &u_temp[3*N];

    //? Internal stage 1: k1 = dt * RHS(u)
    RHS(u, k1);                                     //* k1 = RHS(u) = du/dt
    LeXInt::axpby(dt, k1, k1, N, GPU);              //* k1 = k1 * dt

    //? Internal stage 2: k2 = dt * RHS(u + k1/2)
    LeXInt::axpby(1.0, u, 0.5, k1, u_sol, N, GPU);  //* u_sol = u + k1/2
    RHS(u_sol, k2);                                 //* k2 = RHS(u + k1/2) = RHS(u_sol) 
    LeXInt::axpby(dt, k2, k2, N, GPU);              //* k2 = k2 * dt

    //? Internal stage 3: k3 = dt * RHS(u + k2/2)
    LeXInt::axpby(1.0, u, 0.5, k2, u_sol, N, GPU);  //* u_sol = u + k2/2
    RHS(u_sol, k3);                                 //* k3 = RHS(u + k2/2) = RHS(u_sol) 
    LeXInt::axpby(dt, k3, k3, N, GPU);              //* k3 = k3 * dt

    //? Internal stage 4: k4 = dt * RHS(u + k3)
    LeXInt::axpby(1.0, u, 1.0, k3, u_sol, N, GPU);  //* u_sol = u + k3
    RHS(u_sol, k4);                                 //* k4 = RHS(u + k3) = RHS(u_sol) 
    LeXInt::axpby(dt, k4, k4, N, GPU);              //* k4 = k4 * dt

    //? u_rk4 = u + 1./6.*(k1 + 2*k2 + 2*k3 + k4)
    LeXInt::axpby(1.0, u, 1./6., k1, 2./6., k2, 2./6., k3, 1./6., k4, u_sol, N, GPU);
}