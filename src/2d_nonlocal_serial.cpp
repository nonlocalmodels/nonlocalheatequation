//  Copyright (c) 2014 Hartmut Kaiser
//  Copyright (c) 2014 Patricia Grubel
//  Copyright (c) 2020 Pranav Gadikar
//
//  SPDX-License-Identifier: BSL-1.0
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm> 
#include <fstream>

#include "../include/point.h"
#include "../include/print_time_results.hpp"
#include "../include/writer.h"

//mathematical constant PI
#define PI 3.14159265

bool header = 1;
double k = 1 ;     // heat transfer coefficient
double dt = 1.;     // time step
double dh = 1.;     // grid spacing

class solver
{
public:
    // 3d coordinate for representation of point
    typedef util::Point3 point_3d;

    // 3d plane using the 3d point representation
    // Note that the plane is actually stored in a 1d fashion
    typedef std::vector<point_3d> plane_3d;

    // Our partition type
    typedef double temperature;

    // 2d space data structure
    typedef std::vector<temperature> space_2d;

    //alternate for space between time steps,
    //result of S[t%2] goes to S[(t+1)%2] at every timestep 't'
    space_2d S[2];

    // vector containining 3d coordinates of the points in the plane
    plane_3d P;

    //nx = length of grid in x dimension
    //ny = length of grid in y dimension
    //nt = number of timesteps
    //eps = Epsilon for influence zone of a point
    //nlog = Number of time steps to log the results
    long nx, ny, nt, eps, c_2d, nlog;
    bool current, next, test;
    double error_l2, error_linf;

    // file to store the simulation results in csv format
    const std::string simulate_fname = "../out_csv/simulate_2d.csv";

    // file to store l2 and l infinity norms per timestep
    const std::string score_fname = "../out_csv/score_2d.csv";

    //constructor to initialize class variables and allocate
    solver(long nx, long ny, long nt, long eps, long nlog)
    {
        this->nx = nx;
        this->ny = ny;
        this->nt = nt;
        this->eps = eps;
        this->nlog = nlog;
        this->c_2d = (k * 8)/ pow(eps * dh, 4);
        this->current = 0;
        this->error_l2 = 0.0;
        this->error_linf = 0.0;
        this->test = 0;

        P.resize(nx * ny);
        for(long sx = 0; sx < nx; ++sx)
        {
            for(long sy = 0; sy < ny; ++sy)
            {
                P[get_loc(sx, sy)] = point_3d(sx, sy, 0);
            }
        }

        for(auto &s : S)
        {
            s.resize(nx * ny);
        }
    }

    // function to compute l2 norm of the solution
    void compute_l2(long time)
    {
        error_l2 = 0;
        
        for(long sx = 0; sx < nx; ++sx)
            for(long sy = 0; sy < ny; ++sy)
                error_l2 += (S[next][get_loc(sx, sy)] - w(sx, sy, time)) * (S[next][get_loc(sx, sy)] - w(sx, sy, time));
    }

    // function to compute l infinity norm of the solution
    void compute_linf(long time)
    {
        error_linf = 0;
        
        for(long sx = 0; sx < nx; ++sx)
            for(long sy = 0; sy < ny; ++sy)
                error_linf += std::abs(S[next][get_loc(sx, sy)] - w(sx, sy, time));
    }

    //print error for testing
    void print_error(bool cmp)
    {
        std::cout << "l2: " << error_l2 << " linfinity: " << error_linf << std::endl;
        if(cmp)
            for(long sx = 0; sx < nx; ++sx)
                for(long sy = 0; sy < ny; ++sy)
                    std::cout << "Expected: " << w(sx, sy, nt)
                    << " Actual: " << S[nt % 2][get_loc(sx, sy)] << std::endl;
    }

    //print the solution for the user
    void print_soln()
    {
        for (long sx = 0; sx < nx; ++sx)
        {
            for (long sy = 0; sy < ny; ++sy)
            {
                std::cout << "S[" << sx << "][" << sy << "] = " << S[nt%2][get_loc(sx, sy)] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Function to visualize in csv format and conduct various experiments
    void log_vtk(long log_num)
    {
        const std::string fname = "../out_vtk/simulate_" + std::to_string(log_num);
        rw::writer::VtkWriter vtk_logger(fname);

        vtk_logger.appendNodes(&P);
        vtk_logger.appendPointData("Temperature", &S[next]);
        vtk_logger.addTimeStep(std::time(0));
        vtk_logger.close();
    }

    // Function to visualize in csv format and conduct various experiments
    void log_csv(long time)
    {
        std::ofstream outfile;
        outfile.open(simulate_fname, std::ios_base::app);
        
        for(long sx = 0; sx < nx; ++sx)
        {
            for(long sy = 0; sy < ny; ++sy)
            {
                outfile << time << ","
                << sx << ","
                << sy << ","
                << S[next][get_loc(sx, sy)] << ","
                << w(sx, sy, time) << ","
                << (S[next][get_loc(sx, sy)] - w(sx, sy, time)) * (S[next][get_loc(sx, sy)] - w(sx, sy, time)) << ","
                << std::abs(S[next][get_loc(sx, sy)] - w(sx, sy, time)) << ",\n";
            }
        }

        outfile.close();

        // add to the score csv only when it's required
        if(test)
        {
            compute_l2(time);
            compute_linf(time);

            std::ofstream outfile;

            outfile.open(score_fname, std::ios_base::app);
            outfile << time << "," << error_l2
                << "," << error_linf << ",\n";
            outfile.close();
        }
    }

    //input the initialization for 2d nonlocal equation
    void input_init()
    {
        test = 0;
        for(long sx = 0; sx < nx; ++sx)
        {
            for(long sy = 0; sy < ny; ++sy)
            {
                std::cin >> S[0][get_loc(sx, sy)];
            }
        }
    }

    //init for testing the 2d nonlocal equation
    void test_init()
    {
        test = 1;
        for(long sx = 0; sx < nx; ++sx)
        {
            for(long sy = 0; sy < ny; ++sy)
            {
                S[0][get_loc(sx, sy)] = sin(2 * PI * (sx * dh))
                    *sin(2 * PI * (sy * dh));
            }
        }
    }

    //our influence function to weigh various points differently
    static double influence_function(double distance)
    {
        return 1.0;
    }

    // operator for getting the location of 2d point in 1d representation
    inline long get_loc(long x, long y)
    {
        return x + y * nx;
    }

    //testing operator to verify correctness
    inline double w(long pos_x, long pos_y, long time)
    {
        return cos(2 * PI * (time * dt)) 
                * sin(2 * PI * (pos_x * dh))
                * sin(2 * PI * (pos_y * dh));
    }

    //condition to enforce the boundary conditions
    inline double boundary(long pos_x, long pos_y, double val = 2.0)
    {
        if(pos_x >= 0 && pos_x < nx && pos_y >= 0 && pos_y < ny)
            if(val != 2.0)
                return val;
            else
                return S[current][get_loc(pos_x, pos_y)];
        else
            return 0;
    }

    // Function to find distance in 2d plane given 2 points coordinates
    double distance(long cen_x, long cen_y, long pos_x, long pos_y)
    {
        return sqrt(((cen_x - pos_x) * (cen_x - pos_x))
                +   ((cen_y - pos_y) * (cen_y - pos_y)));
    }

    // Function to compute length of 3rd side of a right angled triangle
    // given hypotenuse and lenght of one side
    double len_1d_line(long len_x)
    {
        return sqrt((eps * eps) - (len_x * len_x));
    }

    //Following function adds the external source to the 2d nonlocal heat equation
    double sum_local_test(long pos_x, long pos_y, long time)
    {
        double result_local = - (2 * PI * sin(2 * PI * (time * dt)) 
                                * sin(2 * PI * (pos_x * dh))
                                * sin(2 * PI * (pos_y * dh)));
        double w_position = w(pos_x, pos_y, time);
        long len_line = 0;

        for(long sx = pos_x-eps; sx <= pos_x+eps; ++sx)
        {
            len_line = len_1d_line(std::abs((long)pos_x - (long)sx));
            for(long sy = pos_y-len_line; sy <= pos_y+len_line; ++sy)
            {
                result_local -= influence_function(distance(pos_x, pos_y, sx, sy))
                                * c_2d
                                * (boundary(sx, sy, w(sx, sy, time)) - w_position)
                                * (dh * dh);
            }
        }
        
        return result_local;
    }

    // Our operator to find sum of 'eps' radius circle in vicinity of point P(x,y)
    // Represent circle as a series of horizaontal lines of thickness 'dh'
    double sum_local(long pos_x, long pos_y)
    {
        double result_local = 0.0;
        long len_line = 0;

        for(long sx = pos_x-eps; sx <= pos_x+eps; ++sx)
        {
            len_line = len_1d_line(std::abs((long)pos_x - (long)sx));
            for(long sy = pos_y-len_line; sy <= pos_y+len_line; ++sy)
            {
                result_local += influence_function(distance(pos_x, pos_y, sx, sy))
                                * c_2d
                                * (boundary(sx, sy) - S[current][get_loc(pos_x, pos_y)])
                                * (dh * dh);
            }
        }
        
        return result_local;
    }

    // do all the work on 'nx * ny' data points for 'nt' time steps
    space_2d do_work()
    {
        // Actual time step loop
        for (long t = 0; t < nt; ++t)
        {
            current = t % 2;
            next = (t + 1) % 2;

            for (long sx = 0; sx < nx; ++sx)
            {
                for (long sy = 0; sy < ny; ++sy)
                {
                    S[next][get_loc(sx, sy)] = S[current][get_loc(sx, sy)] + (sum_local(sx, sy) * dt);
                    if(test)
                        S[next][get_loc(sx, sy)] += sum_local_test(sx, sy, t) * dt;
                }
            }

            if(t%nlog == 0)
            {
                log_vtk(t/ nlog);
                log_csv(t);
            }

        }

        next = nt % 2;

        //testing the code for correctness
        if(test)
        {
            compute_l2(nt);
            compute_linf(nt);
        }

        // Return the solution at time-step 'nt'.
        return S[nt % 2];
    }
};

int batch_tester(long nlog)
{
    std::uint64_t nx, ny, nt, eps, num_tests;
    bool test_failed = 0;
    
    std::cin >> num_tests;
    for(long i = 0; i < num_tests; ++i)
    {
        std::cin >> nx >> ny >> nt >> eps >> k >> dt >> dh;
        
        // Create the solver object
        solver solve(nx, ny, nt, eps, nlog);
        solve.test_init();
        
        solve.do_work();
        
        if (solve.error_l2 / (double)(nx * ny) > 1e-6)
        {
            test_failed = 1;
            break;
        }
    }

    // output regular expression whether the test passed or failed
    if(test_failed)
        std::cout << "Tests Failed" << std::endl;
    else
        std::cout << "Tests Passed" << std::endl;

    return hpx::finalize();
}

int hpx_main(hpx::program_options::variables_map& vm)
{
    std::uint64_t nx = vm["nx"].as<std::uint64_t>();   // Number of grid points in x dimension
    std::uint64_t ny = vm["ny"].as<std::uint64_t>();   // Number of grid points in y dimension
    std::uint64_t nt = vm["nt"].as<std::uint64_t>();   // Number of steps.
    std::uint64_t eps = vm["eps"].as<std::uint64_t>();   // Epsilon for influence zone of a point
    std::uint64_t nlog = vm["nlog"].as<std::uint64_t>();   // Number of time steps to log the results

    if (vm.count("no-header"))
        header = false;

    //batch testing for ctesting
    if (vm.count("test_batch"))
        return batch_tester(nlog);

    // Create the solver object
    solver solve(nx, ny, nt, eps, nlog);

    //Take inputs from stdin for testing
    if(vm.count("test"))
        solve.test_init();
    else
        solve.input_init();

    // Measure execution time.
    std::uint64_t t = hpx::util::high_resolution_clock::now();

    // Execute nt time steps on nx grid points.
    solver::space_2d solution = solve.do_work();

    std::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;

    //print the calculated error
    //print comparison only when the 'cmp' flag is set
    if(vm.count("test"))
        solve.print_error(vm["cmp"].as<bool>());

    // Print the final solution
    if (vm.count("results"))
        solve.print_soln();

    std::uint64_t const os_thread_count = hpx::get_os_thread_count();
    print_time_results(os_thread_count, elapsed, nx, ny, nt, header);

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    namespace po = hpx::program_options;

    po::options_description desc_commandline;
    desc_commandline.add_options()
        ("test", "use arguments from numerical solution for testing (default: false)")
        ("test_batch", "test the solution against numerical solution against batch inputs (default: false)")
        ("results", "print generated results (default: false)")
        ("cmp", po::value<bool>()->default_value(true),
         "Print expected versus actual outputs")
        ("nx", po::value<std::uint64_t>()->default_value(50),
         "Local x dimension")
        ("ny", po::value<std::uint64_t>()->default_value(50),
        "Local y dimension")
        ("nt", po::value<std::uint64_t>()->default_value(45),
         "Number of time steps")
        ("nlog", po::value<std::uint64_t>()->default_value(5),
         "Number of time steps to log the results")
        ("eps", po::value<std::uint64_t>()->default_value(5),
         "Epsilon for nonlocal equation")
        ("k", po::value<double>(&k)->default_value(1),
         "Heat transfer coefficient (default: 0.5)")
        ("dt", po::value<double>(&dt)->default_value(0.0005),
         "Timestep unit (default: 1.0[s])")
        ("dh", po::value<double>(&dh)->default_value(0.02),
         "Quantization across space")
        ( "no-header", "do not print out the csv header row")
    ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}