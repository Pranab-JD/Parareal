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

// TODO: Implement RK2, RK4 (the function declarations have to be similar to the previous one)l