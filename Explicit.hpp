#pragma once

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

#include "functions.hpp"

//? ----------------------------------------------------------
//?
//? Description:
//?     A collection of explicit integrators
//?
//? ----------------------------------------------------------

template <typename rhs>
void explicit_Euler(rhs RHS, double* u, double* u_sol, double* u_temp,double dt, size_t N)
{
    //? RHS(u) = u_temp = du/dt
    RHS(u, u_temp);

    //? u_sol = u + du/dt * dt
    axpby_Cpp(1.0, u, dt, u_temp, u_sol, N);

}

template <typename rhs>
void Runge_Kutta2(rhs RHS, double* u, double* u_sol, double* u_temp, double dt, size_t N)
{
    double* k1 = new double[N];
    double* k2 = new double[N];
  
    RHS(u, k1);
    axpby_Cpp(1.0, u, 0.5*dt, k1, u_temp, N);
    
    RHS(u_temp, k2);
    axpby_Cpp(1.0, u, dt, k2, u_sol, N);

}

template <typename rhs>
void Runge_Kutta4(rhs RHS, double* u, double* u_sol, double* u_temp, double dt, size_t N)
{
    double* k1 = new double[N];
    double* k2 = new double[N];
    double* k3 = new double[N];
    double* k4 = new double[N];

    RHS(u, k1);
    axpby_Cpp(1.0, u, 0.5*dt, k1, u_temp, N);
    
    RHS(u_temp, k2);
    axpby_Cpp(1.0, u, 0.5*dt, k2, u_temp, N);
   
    RHS(u_temp, k3);
    axpby_Cpp(1.0, u, dt, k3, u_temp, N);
    
    RHS(u_temp, k4);
    // axpby_Cpp(1.0, u, 1./6.*dt, k1+2.0*k2+2.0*k3+k4, u_sol, N)
    for (int i=0; i<N; i++)
    {
        u_sol[i] = u[i] + 1./6.*dt*(k1[i]+2.*k2[i]+2.*k3[i]+k4[i]);
    }
}

// TODO: Implement RK2, RK4 (the function declarations have to be similar to the previous one)l
