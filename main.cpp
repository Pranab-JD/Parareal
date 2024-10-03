#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

//? LeXInt Timer and functions
#include "Timer.hpp"
#include "functions.hpp"

//? Problems
#include "Dif_Adv_2D.hpp"

//? Solvers
#include "Explicit.hpp"

using namespace std;

//? ====================================================================================== ?//

int main(int argc, char** argv)
{

    int index = atoi(argv[1]);              // N = 2^index * 2^index
    double n_cfl = atof(argv[2]);           // dt = n_cfl * dt_cfl
    double tol = atof(argv[3]);             // User-specified tolerance
    int num_time_steps = atoi(argv[4]);     // Final simulation time
    
    string movie = "no";                    // Default param = "no"
    if (sizeof(argv) == 5)
        movie = argv[5];                    // Set to "yes" to write data for plots/movie
    
    int num_threads;                        // # of OpenMP threads
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        #pragma omp single
        cout << "Using " << num_threads << "OpenMP threads." << endl ;
    }
    
    //! Set GPU spport to false
    bool GPU_access = false;

    //* Initialise parameters
    int n = pow(2, index);                          // # grid points (1D)
    int N = n*n;                                    // # grid points (2D)
    double xmin = -1;                               // Left boundary (limit)
    double xmax =  1;                               // Right boundary (limit)
    double ymin = -1;                               // Left boundary (limit)
    double ymax =  1;                               // Right boundary (limit)
    vector<double> X(n);                            // Array of grid points
    vector<double> Y(n);                            // Array of grid points
    vector<double> u_init(N);                       // Initial condition

    //* Set up X, Y arrays and initial condition
    for (int ii = 0; ii < n; ii++)
    {
        X[ii] = xmin + ii*(xmax - xmin)/n;
        Y[ii] = ymin + ii*(ymax - ymin)/n;
    }

    //* Initialise additional parameters
    double dx = X[12] - X[11];                              // Grid spacing
    double dy = Y[12] - Y[11];                              // Grid spacing
    double velocity = 10;                                   // Advection speed

    //* Temporal parameters
    double dif_cfl = (dx*dx * dy*dy)/(2*dx*dx + 2*dy*dy);   // Diffusion CFL
    double adv_cfl = dx*dy/(velocity * (dx + dy));          // Advection CFL
    double dt = n_cfl*min(dif_cfl, adv_cfl);                // Step size

    double time = 0;                                        // Simulation time elapsed                          
    int time_steps = 0;                                     // # time steps

    cout << endl << "N = " << N << ", tol = " << tol << ", Time steps = " << num_time_steps << endl;
    cout << "N_cfl = " << n_cfl << ", CFL: " << min(dif_cfl, adv_cfl) << ", dt = " << dt << endl << endl;

    int iters = 0;                                          //* # of iterations per time step
    int iters_total = 0;                                    //* Total # of iterations during the simulation

    //? Choose problem and integrator
    string problem = "Diff_Adv_2D";
    string integrator = "Explicit_Euler";

    //! Diffusion-Advection (+ Sources)
    RHS_Dif_Adv_2D RHS(n, dx, dy, velocity); 

    if (problem == "Diff_Adv_2D")
    {
        //? Initial condition
        for (int ii = 0; ii < n; ii++)
        {
            for (int jj = 0; jj< n; jj++)
            {
                u_init[n*ii + jj] = 1 + exp(-((X[ii] + 0.5)*(X[ii] + 0.5) + (Y[jj] + 0.5)*(Y[jj] + 0.5))/0.01);
            }
        }
    }
    else
    {
        cout << "Undefined problem! Please check that you have entered the correct problem." << endl;
    } 

    //! Allocate memory on CPU
    size_t N_size = N * sizeof(double);
    double* u = (double*)malloc(N_size);
    copy(u_init.begin(), u_init.end(), u);
    double* u_sol = (double*)malloc(N_size);            //* Solution vector
    double* u_temp;                                     //* Temporary vector(s) for integrators

    if (integrator == "Explicit_Euler")
    {
        u_temp = (double*)malloc(N_size);
    }
    else if (integrator == "RK2")
    {
        u_temp = (double*)malloc(2*N_size);
    }
    else if (integrator == "RK4")
    {
        u_temp = (double*)malloc(4*N_size);
    }
    else if (integrator == "Implicit_Euler")
    {
        u_temp = (double*)malloc(2*N_size);
    }
    else if (integrator == "CN")
    {
        u_temp = (double*)malloc(2*N_size);
    }
    else
    {
        cout << "Incorrect integrator!";
    }


    //! Create directories (for movies)
    if (movie == "yes")
    {
        int sys_value = system(("mkdir -p ./movie/"));
        string directory = "./movie/";
    }

    //! Time Loop
    LeXInt::timer time_loop;
    time_loop.start();

    for (int nn = 0; nn < num_time_steps; nn++)
    {

        //? ------------- List of integrators ------------- ?//

        if (integrator == "Explicit_Euler")
        {
            explicit_Euler(RHS, u, u_sol, u_temp, dt, N);
        }
        else if (integrator == "Implicit_Euler")
        {

        }

        //? ----------------------------------------------- ?//

        //* Update variables
        time = time + dt;
        time_steps = time_steps + 1;
        iters_total = iters_total + iters;

        //? Update solution
        copy_Cpp(u_sol, u, N);

        if (time_steps % 5 == 0)
        {
            cout << "Time step      : " << time_steps << endl;
            cout << "Simulation time: " << time << endl << endl;
        }

        //! Write data to files (for movies)
        if (time_steps % 100 == 0 && movie == "yes")
        {
        
            string output_data = "./movie/" +  to_string(time_steps) + ".txt";
            ofstream data;
            data.open(output_data); 
            for(int ii = 0; ii < N; ii++)
            {
                data << setprecision(16) << u[ii] << endl;
            }
            data.close();
        }
    }

    time_loop.stop();

    cout << endl << "==================================================" << endl;
    cout << "Simulation time            : " << time << endl;
    cout << "Total number of time steps : " << time_steps << endl;
    cout << "Total number of iterations : " << iters_total << endl;
    cout << "Total time elapsed (s)     : " << time_loop.total() << endl;
    cout << "==================================================" << endl << endl;


    //? Create directory to write simulation results/parameters
    int sys_value = system(("mkdir -p ./" + integrator + "/cores_" + to_string(num_threads)).c_str());
    string directory = "./" + integrator + "/cores_" + to_string(num_threads);
    string results = directory + "/Parameters.txt";
    ofstream params;
    params.open(results);
    params << "Grid points: " << N << endl;
    params << "Step size: " << dt << endl;
    params << "Tolerance (for implicit methods): " << tol << endl;
    params << "Simulation time: " << time << endl;
    params << "Total number of time steps: " << time_steps << endl;
    params << "Number of OpenMP threads: " << num_threads << endl;
    params << endl;
    params << "Total iterations (for implicit methods): " << iters_total << endl;
    params << setprecision(16) << "Runtime (s): " << time_loop.total() << endl;
    params.close();

    cout << "Simulations complete!" << endl;

    return 0;
}