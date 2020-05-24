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

bool header = 1;
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
    long nx, ny, nt, eps;
    bool current, next;

    //constructor to initialize class variables and allocate
    solver(long nx, long ny, long nt, long eps)
    {
        this->nx = nx;
        this->ny = ny;
        this->nt = nt;
        this->eps = eps;
        this->current = 0;

        for(auto &sx : S)
        {
            sx.resize(nx);
            for(auto &sy : sx)
            {
                sy.resize(ny);
            }
        }
    }

    //init for testing the 1d nonlocal equation
    void init_space()
    {
        for(auto &sx : S[0])
        {
            for(auto &sy : sx)
            {
                std::cin >> sy;
            }
        }
    }

    //default init for space to check correctness
    void default_init()
    {
        for(long i = 0; i < nx; ++i)
        {
            for(long j = 0; j < ny; ++j)
            {
                S[0][i][j] = i;
            }
        }
    }

    //our influence function to weigh various points differently
    static double influence_function(double distance)
    {
        return 1.0;
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

    // Our operator to find sum of 'eps' radius circle in vicinity of point P(x,y)
    // Represent circle as a series of horizaontal lines of thickness 'dh'
    double sum_local(long pos_x, long pos_y)
    {
        double temp = 0.0;
        long len_line = 0;
        // std::cout<<std::max((long)0,pos_x-eps) << " "<< std::min(pos_x+eps,nx-1) <<"\n";
        for(long i = std::max((long)0,pos_x-eps); i <= std::min(pos_x+eps,nx-1); ++i)
        {
            len_line = len_1d_line(std::abs((long)pos_x - (long)i));
            for(long j = std::max((long)0,pos_y-len_line); j <= std::min(pos_y+len_line,ny-1); ++j)
            {
                temp += influence_function(distance(pos_x, pos_y, i, j))
                                * (S[current][i][j] - S[current][pos_x][pos_y]);
                // std::cout << temp <<" ";
            }
            // std::cout<<"\n";
        }
        temp *= ((1 / (pow(eps*dh,4))) * (dh * dh));
        return temp;
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
                }
            }

        }

        // std::cout<<"\n" << sum_local(0,0) << "\n";
        // Return the solution at time-step 'nt'.
        return S[nt % 2];
    }
};


int hpx_main(hpx::program_options::variables_map& vm)
{
    std::uint64_t nx = vm["nx"].as<std::uint64_t>();   // Number of grid points in x dimension
    std::uint64_t ny = vm["ny"].as<std::uint64_t>();   // Number of grid points in y dimension
    std::uint64_t nt = vm["nt"].as<std::uint64_t>();   // Number of steps.
    std::uint64_t eps = vm["eps"].as<std::uint64_t>();   // Epsilon for influence zone of a point

    if (vm.count("no-header"))
        header = false;

    // Create the solver object
    solver solve(nx, ny, nt, eps);

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
    solver::space solution = solve.do_work();

    // Print the final solution
    if (vm.count("results"))
    {
        for (long i = 0; i < nx; ++i)
        {
            for (long j = 0; j < ny; ++j)
            {
                std::cout << "S[" << i << "][" << j << "] = " << solution[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    std::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;

    std::uint64_t const os_thread_count = hpx::get_os_thread_count();
    print_time_results(os_thread_count, elapsed, nx, ny, nt, header);

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    namespace po = hpx::program_options;

    po::options_description desc_commandline;
    desc_commandline.add_options()
        ("test", "use arguments from stdin for testing (default: false)")
        ("results", "print generated results (default: false)")
        ("nx", po::value<std::uint64_t>()->default_value(10),
         "Local x dimension")
        ("ny", po::value<std::uint64_t>()->default_value(10),
        "Local y dimension")
        ("nt", po::value<std::uint64_t>()->default_value(45),
         "Number of time steps")
        ("eps", po::value<std::uint64_t>()->default_value(5),
         "Epsilon for nonlocal equation")
        ("dt", po::value<double>(&dt)->default_value(1.0),
         "Timestep unit (default: 1.0[s])")
        ("dh", po::value<double>(&dh)->default_value(1.0),
         "quantization across space")
        ( "no-header", "do not print out the csv header row")
    ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}