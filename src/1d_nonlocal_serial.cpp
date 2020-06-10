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

bool header = 0;
double k = 1 ;     // heat transfer coefficient
double dt = 1.;     // time step
double dx = 1.;     // grid spacing

class solver
{
public:
    // Our partition type
    typedef double partition;

    // Our data for one time step
    typedef std::vector<partition> space;

    //alternate for space between time steps,
    //result of S[t%2] goes to S[(t+1)%2] at every timestep 't'
    space S[2];

    //nx = number of data points
    //nt = number of timesteps
    //eps = Epsilon for influence zone of a point
    long nx, nt, eps, c_1d;
    bool current, next, test;
    double error;

    //constructor to initialize class variables
    solver(long nx, long nt, long eps)
    {
        this->nx = nx;
        this->nt = nt;
        this->eps = eps;
        this->c_1d = (k * 3)/ (pow(eps * dx, 3));
        this->current = 0;
        this->error = 0.0;
        this->test = 0;

        for(auto &s : S)
        {
            s.resize(nx);
        }
    }

    //print error for testing
    void print_error(bool cmp)
    {
        std::cout << error << std::endl;
        if(cmp)
            for(long i = 0; i < nx; ++i)
                std::cout << "Expected: " << w(i, nt)
                << " Actual: " << S[nt % 2][i] << std::endl;
    }

    //input the initialization 
    void input_init()
    {
        this->test = 0;
        for(auto &i : S[0])
        {
            std::cin>>i;
        }
    }

    //init for testing the 1d nonlocal equation
    void test_init()
    {
        this->test = 1;
        for(long i = 0; i < nx; ++i)
        {
            S[0][i] = sin(2 * PI * (i * dx));
        }
    }

    //our influence function to weigh various points differently
    static double influence_function(double distance)
    {
        return 1;
    }

    //testing operator to verify correctness
    inline double w(long position, long time)
    {
        return cos(2 * PI * (time * dt)) * sin(2 * PI * (position * dx));
    }

    //condition to enforce the boundary conditions
    inline double boundary(long i, double val)
    {
        if(i >= 0 && i < nx)
            return val;
        else
            return 0;
    }

    //our test by introducing an external source
    double sum_local_test(long position, long time)
    {
        double result_local = - (2 * PI * sin(2 * PI * (time * dt)) * sin(2 * PI * (position * dx)));
        double w_position = w(position, time);
        for(long i = position-eps; i <= position+eps; ++i)
        {
            result_local -= influence_function(std::abs((long)position - (long)i)) 
                        * c_1d
                        * (boundary(i, w(i, time)) - w_position) 
                        * (dx);
        }
        return result_local;
    }

    // Our operator for getting result from local points on the locality
    double sum_local(long position)
    {
        double result_local = 0.0;
        double s_i = 0.0;
        for(long i = position-eps; i <= position+eps; ++i)
        {
            result_local += influence_function(std::abs((long)position - (long)i)) 
                        * c_1d
                        * (boundary(i, S[current][i]) - S[current][position]) 
                        * (dx);
        }
        return result_local;
    }

    // do all the work on 'nx' data points for 'nt' time steps
    space do_work()
    {
        // Actual time step loop
        for (long t = 0; t < nt; ++t)
        {
            current = t % 2;
            next = (t + 1) % 2;

            for (long i = 0; i < nx; ++i)
            {
                S[next][i] = S[current][i] + (sum_local(i) * dt);
                if(test)
                    S[next][i] += sum_local_test(i, t) * dt;
            }

        }

        next = nt % 2;

        //testing the code for correctness
        if(test)
            for(long i = 0; i < nx; ++i)
                error += (S[next][i] - w(i, nt)) * (S[next][i] - w(i, nt));
                
        // Return the solution at time-step 'nt'.
        return S[next];
    }
};

int batch_tester()
{
    std::uint64_t nx, nt, eps, num_tests;
    bool test_failed = 0;
    
    std::cin >> num_tests;
    for(long i = 0; i < num_tests; ++i)
    {
        std::cin >> nx >> nt >> eps >> k >> dt >> dx;
        
        // Create the solver object
        solver solve(nx, nt, eps);
        solve.test_init();
        
        solve.do_work();
        
        if (solve.error / (double)nx > 1e-6)
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
    std::uint64_t nx = vm["nx"].as<std::uint64_t>();   // Number of grid points.
    std::uint64_t nt = vm["nt"].as<std::uint64_t>();   // Number of steps.
    std::uint64_t eps = vm["eps"].as<std::uint64_t>();   // Epsilon for influence zone of a point

    if (vm.count("no-header"))
        header = false;
    
    //batch testing for ctesting
    if (vm.count("test_batch"))
        return batch_tester();

    // Create the solver object
    solver solve(nx, nt, eps);

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

    //print the calculated error for testing
    if(vm.count("test"))
        solve.print_error(vm["cmp"].as<bool>());

    // Print the final solution
    if (vm.count("results"))
        for (long i = 0; i < nx; ++i)
            std::cout << "S[" << i << "] = " << solution[i] << std::endl;

    std::uint64_t const os_thread_count = hpx::get_os_thread_count();
    print_time_results(os_thread_count, elapsed, nx, nt, header);

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
        ("nt", po::value<std::uint64_t>()->default_value(45),
         "Number of time steps")
        ("eps", po::value<std::uint64_t>()->default_value(5),
         "Epsilon for nonlocal equation")
        ("k", po::value<double>(&k)->default_value(1),
         "Heat transfer coefficient (default: 0.5)")
        ("dt", po::value<double>(&dt)->default_value(0.001),
         "Timestep unit (default: 1.0[s])")
        ("dx", po::value<double>(&dx)->default_value(0.02),
         "Local x dimension")
        ( "no-header", "do not print out the csv header row")
    ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}
