//  Copyright (c) 2014 Hartmut Kaiser
//  Copyright (c) 2014 Patricia Grubel
//  Copyright (c) 2014 Pranav Gadikar
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

bool header = 0;
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
    std::size_t nx, nt, eps;
    bool current, next;

    //constructor to initialize class variables
    solver(std::size_t nx, std::size_t nt, std::size_t eps)
    {
        this->nx = nx;
        this->nt = nt;
        this->eps = eps;
        this->current = 0;

        for(auto &s : S)
        {
            s.resize(nx);
        }
    }

    //init for testing the 1d nonlocal equation
    void init_space()
    {
        for(auto &i : S[0])
        {
            std::cin>>i;
        }
    }

    //default init for space
    void default_init()
    {
        for(std::size_t i = 0; i < nx; ++i)
        {
            S[0][i] = i;
        }
    }

    //our influence function
    static double influence_function(double distance)
    {
        return 1.0;
    }

    // Our operator
    double sum_local(std::size_t position)
    {
        double temp = 0.0;
        for(std::size_t i = std::max((std::size_t)0,position-eps); i <= std::min(position+eps,nx-1); ++i)
        {
            temp += influence_function(std::abs((long)position - (long)i)) 
                        * (1 / (pow(eps*dx,4))) 
                        * (S[current][i] - S[current][position]) 
                        * (dx);
        }
        return temp;
    }

    // do all the work on 'nx' data points for 'nt' time steps
    space do_work(std::size_t nx, std::size_t nt)
    {
        // Actual time step loop
        for (std::size_t t = 0; t < nt; ++t)
        {
            current = t % 2;
            next = (t + 1) % 2;

            for (std::size_t i = 0; i < nx; ++i)
                S[next][i] = S[current][i] + (sum_local(i) * dt);

        }

        // Return the solution at time-step 'nt'.
        return S[nt % 2];
    }
};


int hpx_main(hpx::program_options::variables_map& vm)
{
    std::uint64_t nx = vm["nx"].as<std::uint64_t>();   // Number of grid points.
    std::uint64_t nt = vm["nt"].as<std::uint64_t>();   // Number of steps.
    std::uint64_t eps = vm["eps"].as<std::uint64_t>();   // Epsilon for influence zone of a point

    if (vm.count("no-header"))
        header = false;

    // Create the solver object
    solver solve(nx, nt, eps);

    //Take inputs from stdin for testing
    if(vm.count("test"))
    {
        solve.init_space();
    }
    else
    {
        solve.default_init();
    }

    // Measure execution time.
    std::uint64_t t = hpx::util::high_resolution_clock::now();

    // Execute nt time steps on nx grid points.
    solver::space solution = solve.do_work(nx, nt);

    // Print the final solution
    if (vm.count("results"))
    {
        for (std::size_t i = 0; i < nx; ++i)
            std::cout << "S[" << i << "] = " << solution[i] << std::endl;
    }

    std::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;

    std::uint64_t const os_thread_count = hpx::get_os_thread_count();
    print_time_results(os_thread_count, elapsed, nx, nt, header);

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    namespace po = hpx::program_options;

    po::options_description desc_commandline;
    desc_commandline.add_options()
        ("test", "use arguments from stdin for testing (default: false)")
        ("results", "print generated results (default: false)")
        ("nx", po::value<std::uint64_t>()->default_value(100),
         "Local x dimension")
        ("nt", po::value<std::uint64_t>()->default_value(45),
         "Number of time steps")
        ("eps", po::value<std::uint64_t>()->default_value(5),
         "Epsilon for nonlocal equation")
        ("dt", po::value<double>(&dt)->default_value(1.0),
         "Timestep unit (default: 1.0[s])")
        ("dx", po::value<double>(&dx)->default_value(1.0),
         "Local x dimension")
        ( "no-header", "do not print out the csv header row")
    ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}
