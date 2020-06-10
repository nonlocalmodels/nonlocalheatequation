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

#include "print_time_results.hpp"

//mathematical constant PI
#define PI 3.14159265

bool header = 1;
double k = 1 ;     // heat transfer coefficient
double dt = 1.;     // time step
double dh = 1.;     // grid spacing

class solver
{
public:
    // Our partition type
    typedef double partition;

    // 1d space data structure
    typedef std::vector<partition> d1_space;

    // Our data for one time step in 2d space
    typedef std::vector<d1_space> space;

    //alternate for space between time steps,
    //result of S[t%2] goes to S[(t+1)%2] at every timestep 't'
    space S[2];

    //nx = length of grid in x dimension
    //ny = length of grid in y dimension
    //nt = number of timesteps
    //eps = Epsilon for influence zone of a point
    long nx, ny, nt, eps, c_2d;
    bool current, next, test;
    double error;

    //constructor to initialize class variables and allocate
    solver(long nx, long ny, long nt, long eps)
    {
        this->nx = nx;
        this->ny = ny;
        this->nt = nt;
        this->eps = eps;
        this->c_2d = (k * 8)/ pow(eps * dh, 4);
        this->current = 0;
        this->error = 0.0;
        this->test = 0;

        for(auto &sx : S)
        {
            sx.resize(nx);
            for(auto &sy : sx)
            {
                sy.resize(ny);
            }
        }
    }

    //print error for testing
    void print_error(bool cmp)
    {
        std::cout << error << std::endl;
        if(cmp)
            for(long i = 0; i < nx; ++i)
                for(long j = 0; j < ny; ++j)
                    std::cout << "Expected: " << w(i, j, nt)
                    << " Actual: " << S[nt % 2][i][j] << std::endl;
    }

    //input the initialization for 2d nonlocal equation
    void input_init()
    {
        test = 0;
        for(auto &sx : S[0])
        {
            for(auto &sy : sx)
            {
                std::cin >> sy;
            }
        }
    }

    //init for testing the 2d nonlocal equation
    void test_init()
    {
        test = 1;
        for(long i = 0; i < nx; ++i)
        {
            for(long j = 0; j < ny; ++j)
            {
                S[0][i][j] = sin(2 * PI * (i * dh))
                    *sin(2 * PI * (j * dh));
            }
        }
    }

    //our influence function to weigh various points differently
    static double influence_function(double distance)
    {
        return 1.0;
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
                return S[current][pos_x][pos_y];
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

        for(long i = pos_x-eps; i <= pos_x+eps; ++i)
        {
            len_line = len_1d_line(std::abs((long)pos_x - (long)i));
            for(long j = pos_y-len_line; j <= pos_y+len_line; ++j)
            {
                result_local -= influence_function(distance(pos_x, pos_y, i, j))
                                * c_2d
                                * (boundary(i, j, w(i, j, time)) - w_position)
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

        for(long i = pos_x-eps; i <= pos_x+eps; ++i)
        {
            len_line = len_1d_line(std::abs((long)pos_x - (long)i));
            for(long j = pos_y-len_line; j <= pos_y+len_line; ++j)
            {
                result_local += influence_function(distance(pos_x, pos_y, i, j))
                                * c_2d
                                * (boundary(i, j) - S[current][pos_x][pos_y])
                                * (dh * dh);
            }
        }
        
        return result_local;
    }

    // do all the work on 'nx * ny' data points for 'nt' time steps
    space do_work()
    {
        // Actual time step loop
        for (long t = 0; t < nt; ++t)
        {
            current = t % 2;
            next = (t + 1) % 2;

            for (long i = 0; i < nx; ++i)
            {
                for (long j = 0; j < ny; ++j)
                {
                    S[next][i][j] = S[current][i][j] + (sum_local(i, j) * dt);
                    if(test)
                        S[next][i][j] += sum_local_test(i, j, t) * dt;
                }
            }

        }

        next = nt % 2;

        //testing the code for correctness
        if(test)
            for(long i = 0; i < nx; ++i)
                for(long j = 0; j < ny; ++j)
                    error += (S[next][i][j] - w(i, j, nt)) * (S[next][i][j] - w(i, j, nt));

        // Return the solution at time-step 'nt'.
        return S[nt % 2];
    }
};

int batch_tester()
{
    std::uint64_t nx, ny, nt, eps, num_tests;
    bool test_failed = 0;
    
    std::cin >> num_tests;
    for(long i = 0; i < num_tests; ++i)
    {
        std::cin >> nx >> ny >> nt >> eps >> k >> dt >> dh;
        
        // Create the solver object
        solver solve(nx, ny, nt, eps);
        solve.test_init();
        
        solve.do_work();
        
        if (solve.error / (double)(nx * ny) > 1e-6)
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

    if (vm.count("no-header"))
        header = false;

    //batch testing for ctesting
    if (vm.count("test_batch"))
        return batch_tester();

    // Create the solver object
    solver solve(nx, ny, nt, eps);

    //Take inputs from stdin for testing
    if(vm.count("test"))
        solve.test_init();
    else
        solve.input_init();

    // Measure execution time.
    std::uint64_t t = hpx::util::high_resolution_clock::now();

    // Execute nt time steps on nx grid points.
    solver::space solution = solve.do_work();

    std::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;

    //print the calculated error
    //print comparison only when the 'cmp' flag is set
    if(vm.count("test"))
        solve.print_error(vm["cmp"].as<bool>());

    // Print the final solution
    if (vm.count("results"))
        for (long i = 0; i < nx; ++i)
        {
            for (long j = 0; j < ny; ++j)
            {
                std::cout << "S[" << i << "][" << j << "] = " << solution[i][j] << " ";
            }
            std::cout << std::endl;
        }

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