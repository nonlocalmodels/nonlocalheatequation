//  Copyright (c) 2014 Hartmut Kaiser
//  Copyright (c) 2014 Patricia Grubel
//  Copyright (c) 2020 Pranav Gadikar
//
//  SPDX-License-Identifier: BSL-1.0
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Note on terminology  :-
// grid square  : small squares/ rectangles into which we have divided the larger mesh 
// mesh points  : points within a single grid square which is smallest unit of space in our equation

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include <hpx/include/parallel_algorithm.hpp>
#include <boost/range/irange.hpp>

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <algorithm>
#include <fstream>

#include "../include/point.h"
#include "../include/print_time_results.hpp"
#include "../include/writer.h"

bool header = 1;
double k = 1 ;      // heat transfer coefficient
double dt = 1.;     // time step
double dh = 1.;     // grid spacing
long nx = 50;       // Number of grid points in x dimension
long ny = 50;       // Number of grid points in y dimension
long np = 5;        // Number of partitions we want to break the grid into
long eps = 5;       // Epsilon for influence zone of a point
double c_2d = 0.0;  // 'c' to be used in nonlocal equation for 2d case when J(influence function) = 1 

// operator for getting the location of 2d point in 1d representation
inline long get_loc(long x, long y, long nx = nx)
{
    return x + y * nx;
}

// operator for getting 2d coordinate of grid square in which the space has been divided
inline long grid_loc(long x, long y)
{
    return get_loc(x / nx, y / ny, np);
}

// operator for getting 2d coordinate of a point within the grid square
inline long mesh_loc(long x, long y)
{
    return get_loc(x % nx, y % ny);
}

// A small square part of the original space whose temperature values will be computed
// on a single thread
class partition_space
{
    public:

    // coordinates of the grid point corresponding the np * np squares in the partition
    long gx, gy;

    explicit partition_space()
      : data_(new double[0]), size_(0), gx(0), gy(0)
    {}

    partition_space(std::size_t size, long gx, long gy)
      : data_(new double[size]), size_(size), gx(gx), gy(gy)
    {}

    // add test initialization here for executing tests
    partition_space(long nx, long ny, long gx, long gy)
      : data_(new double[nx * ny]),
        size_(nx * ny),
        gx(gx), gy(gy)
    {
        for(long sx = gx*nx; sx < nx + gx*nx; ++sx)
        {
            for(long sy = gy*ny; sy < ny + gy*ny; ++sy)
            {
                data_[mesh_loc(sx, sy)] = sin(2 * M_PI * (sx * dh))
                    *sin(2 * M_PI * (sy * dh));
            }
        }
    }

    // operators to return the data at some index from the given partition
    double& operator[](std::size_t idx) { return data_[idx]; }
    double operator[](std::size_t idx) const { return data_[idx]; }

    std::size_t size() const { return size_; }

private:
    std::shared_ptr<double[]> data_;
    std::size_t size_;

};

// partitioned data available for threads for reading only for efficient accesses
// and remove the need for unordered_map in the indexing part
typedef std::vector<partition_space> space_2d;
space_2d temperature_grid[2];

class solver
{
public:
    // 3d coordinate for representation of point
    typedef util::Point3 point_3d;

    // 3d plane using the 3d point representation
    // Note that the plane is actually stored in a 1d fashion
    typedef std::vector<point_3d> plane_3d;

    // Our data type used for temperature values
    typedef double temperature;

    // partitioning the 2d space into np * np small squares which constitute a partition of space
    typedef hpx::shared_future<partition_space> partition_fut;

    // 2d space data structure which consists of future of square partitions
    typedef std::vector<partition_fut> space_2d_fut;

    //alternate for space between time steps,
    //result of S[t%2] goes to S[(t+1)%2] at every timestep 't'
    space_2d_fut S[2];

    // vector containining 3d coordinates of the points in the plane
    plane_3d P;

    // nt = number of timesteps
    // nd = Depth till which dependency tree is allowed to grow
    // nlog = Number of time steps to log the results
    long nt, nd, nlog;
    bool current, next, test;
    // l2 norm and l infinity norm
    double error_l2, error_linf;

    //constructor to initialize class variables and allocate
    solver(long nt, long nd, long nlog, bool test)
    {
        this->nt = nt;
        this->nd = nd;
        this->nlog = nlog;
        this->current = 0;
        this->error_l2 = 0.0;
        this->error_linf = 0.0;
        this->test = test;
        this->next = 1;
        // setting the global variable c_2d which is a constant for this test case
        // as per the inputs
        c_2d = (k * 8)/ pow(eps * dh, 4);

        P.resize(nx * ny * np * np);
        for(long sx = 0; sx < nx*np; ++sx)
        {
            for(long sy = 0; sy < ny*np; ++sy)
            {
                P[sx + sy * (nx * np)] = point_3d(sx, sy, 0);
            }
        }

        for(space_2d_fut &s : S)
            s.resize(np * np);
        
        // resizing the global temperature grid for quick access for threads
        for(space_2d &t : temperature_grid)
            t.resize(np * np);

        auto range = boost::irange((long)0, np * np);
        hpx::parallel::for_each(hpx::parallel::execution::par, std::begin(range), std::end(range),
            [this](std::size_t i)
            {
                this->S[0][i] = hpx::make_ready_future(partition_space(nx, ny, i % np, i / np));
                temperature_grid[0][i] = hpx::util::unwrap(S[0][i]);
            }
        );
    }

    // function to compute l2 norm of the solution
    void compute_l2(long time)
    {
        error_l2 = 0;

        for(long sx = 0; sx < nx * np; ++sx)
            for(long sy = 0; sy < ny * np; ++sy)
                error_l2 += std::pow(S[next][grid_loc(sx, sy)].get()[mesh_loc(sx, sy)] - w(sx, sy, time), 2);
    }

    // function to compute l infinity norm of the solution
    void compute_linf(long time)
    {
        error_linf = 0;
        
        for(long sx = 0; sx < nx * np; ++sx)
            for(long sy = 0; sy < ny * np; ++sy)
                error_linf = std::max(std::abs(S[next][grid_loc(sx, sy)].get()[mesh_loc(sx, sy)] - w(sx, sy, time)), error_linf);
    }

    //print error for testing
    void print_error(bool cmp)
    {
        std::cout << "l2: " << error_l2 << " linfinity: " << error_linf << std::endl;
        if(cmp)
            for(long sx = 0; sx < nx * np; ++sx)
                for(long sy = 0; sy < ny * np; ++sy)
                    std::cout << "Expected: " << w(sx, sy, nt)
                    << " Actual: " << S[next][grid_loc(sx, sy)].get()[mesh_loc(sx, sy)] 
                    << std::endl;
    }

    //print the solution for the user
    void print_soln()
    {
        for (long sx = 0; sx < nx*ny; ++sx)
        {
            for (long sy = 0; sy < ny*ny; ++sy)
            {
                std::cout << "S[" << sx << "][" << sy << "] = " 
                << S[next][grid_loc(sx, sy)].get()[mesh_loc(sx, sy)] << " ";
            }
            std::cout << std::endl;
        }
    }

    // function to perform the logging in both csv and vtk format
    static void log_csv_vtk(std::vector<partition_space> temp_data, plane_3d coord, long time, bool test)
    {
        // file to store the simulation results in csv format
        const std::string simulate_fname = "../out_csv/simulate_2d.csv";

        // file to store l2 and l infinity norms per timestep
        const std::string score_fname = "../out_csv/score_2d.csv";
        
        std::ofstream outfile;
        outfile.open(simulate_fname, std::ios_base::app);
        
        for(long sx = 0; sx < nx*np; ++sx)
        {
            for(long sy = 0; sy < ny*np; ++sy)
            {
                outfile << time << ","
                << sx << ","
                << sy << ","
                << temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] << ","
                << w(sx, sy, time) << ","
                << std::pow((temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] - w(sx, sy, time)), 2) << ","
                << std::abs(temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] - w(sx, sy, time)) << ",\n";
            }
        }

        outfile.close();

        // VTK writer begins here

        // vtu file name for writing the vtk data
        const std::string fname = "../out_vtk/simulate_" + std::to_string(time);
        rw::writer::VtkWriter vtk_logger(fname);
        std::vector<double> vector_data(nx*ny*np*np);

        vtk_logger.appendNodes(&coord);
        
        for(long sx = 0; sx < nx*np; ++sx)
            for(long sy = 0; sy < ny*np; ++sy)
                vector_data[sy*nx*np + sx] = temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)];
                
        vtk_logger.appendPointData("Temperature", &vector_data);
        vtk_logger.addTimeStep(std::time(0));
        vtk_logger.close();

        // VTK writer ends here

        // add to the score csv only when it's required
        if(test)
        {
            double error_l2 = 0, error_linf = 0;

            for(long sx = 0; sx < nx * np; ++sx)
            {
                for(long sy = 0; sy < ny * np; ++sy)
                {
                    error_linf = std::max(std::abs(temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] - w(sx, sy, time)), error_linf);
                    error_l2 += std::pow(temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] - w(sx, sy, time), 2);
                }
            }

            std::ofstream outfile;

            outfile.open(score_fname, std::ios_base::app);
            outfile << time << "," << error_l2
                << "," << error_linf << ",\n";
            outfile.close();
        }
    }

    //our influence function to weigh various points differently
    static inline double influence_function(double distance)
    {
        return 1.0;
    }

    // inserting the rectangles which are just above, below and on the sides of the given partition
    // into the unordered_map used to index the squares already inserted into the vector
    static void add_neighbour_rectangle(long lx, long ly, long rx, long ry, space_2d_fut &all_squares
                       , space_2d_fut &neighbour_squares)
    {
        for(long sx = lx; sx < rx+nx; sx+=nx)
        {
            for(long sy = ly; sy < ry+ny; sy+=ny)
            {
                if(sx < 0 || sx >= nx*np || sy < 0 || sy >= ny*np) continue;

                neighbour_squares.push_back(all_squares[grid_loc(sx, sy)]);
            }
        }
    }

    //testing operator to verify correctness
    static inline double w(long pos_x, long pos_y, long time)
    {
        return cos(2 * M_PI * (time * dt)) 
                * sin(2 * M_PI * (pos_x * dh))
                * sin(2 * M_PI * (pos_y * dh));
    }

    //condition to enforce the boundary conditions for the divided partitions
    static inline double boundary(long pos_x, long pos_y, long next)
    {
        if(pos_x >= 0 && pos_x < nx * np && pos_y >= 0 && pos_y < ny * np)
        {
            return temperature_grid[next][grid_loc(pos_x, pos_y)][mesh_loc(pos_x, pos_y)];
        }
        else
        {
            return 0;
        }
    }

    //condition to enforce the boundary conditions for the tests
    static inline double boundary(long pos_x, long pos_y, double val)
    {
        if(pos_x >= 0 && pos_x < nx * np && pos_y >= 0 && pos_y < ny * np)
            return val;
        else
            return 0;
    }

    // Function to find distance in 2d plane given 2 points coordinates
    static inline double distance(long cen_x, long cen_y, long pos_x, long pos_y)
    {
        return sqrt(((cen_x - pos_x) * (cen_x - pos_x))
                +   ((cen_y - pos_y) * (cen_y - pos_y)));
    }

    // Function to compute length of 3rd side of a right angled triangle
    // given hypotenuse and lenght of one side
    static inline double len_1d_line(long len_x)
    {
        return sqrt((eps * eps) - (len_x * len_x));
    }

    //Following function adds the external source to the 2d nonlocal heat equation
    static double sum_local_test(long pos_x, long pos_y, long time)
    {
        double result_local = - (2 * M_PI * sin(2 * M_PI * (time * dt)) 
                                * sin(2 * M_PI * (pos_x * dh))
                                * sin(2 * M_PI * (pos_y * dh)));
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
    static double sum_local(long pos_x, long pos_y, long current)
    {
        double result_local = 0.0;
        long len_line = 0;
        double pos_val = boundary(pos_x, pos_y, current);

        for(long sx = pos_x-eps; sx <= pos_x+eps; ++sx)
        {
            len_line = len_1d_line(std::abs((long)pos_x - (long)sx));
            for(long sy = pos_y-len_line; sy <= pos_y+len_line; ++sy)
            {
                result_local += influence_function(distance(pos_x, pos_y, sx, sy))
                                * c_2d
                                * (boundary(sx, sy, current) - pos_val)
                                * (dh * dh);
            }
        }
        
        return result_local;
    }

    // partition of nx * ny cells which will be processed by a single thread
    static partition_space sum_local_partition(partition_space middle, const std::vector<partition_space> &neighbour_squares, 
                                               long time, bool test)
    {
        long size = middle.size();
        long gx = middle.gx, gy = middle.gy;
        partition_space &next = temperature_grid[(time+1)%2][get_loc(gx, gy, np)];

        next = partition_space(nx, ny, gx, gy);

        for(long sx = 0; sx < nx; ++sx)
        {
            for(long sy = 0; sy < ny; ++sy)
            {
                next[get_loc(sx, sy)] = middle[get_loc(sx, sy)] + (sum_local(sx + gx*nx, sy + gy*ny, time%2) * dt);
                if(test)
                    next[get_loc(sx, sy)] += sum_local_test(sx + gx*nx, sy + gy*ny, time) * dt;
            }
        }

        return next;
    }

    // do all the work on 'nx * ny' data points for each of 'np * np' partitions for 'nt' time steps
    void do_work()
    {
        // limit depth of dependency tree
        hpx::lcos::local::sliding_semaphore sem(nd);

        auto Op = hpx::util::unwrapping(&solver::sum_local_partition);

        // Actual time step loop
        for (long t = 0; t < nt; ++t)
        {
            current = t % 2;
            next = (t + 1) % 2;

            for(long gx = 0; gx < np; ++gx)
            {
                for(long gy = 0; gy < np; ++gy)
                {
                    // vector of dependent grid squares for a particular grid squares
                    space_2d_fut vector_deps;
                    
                    // add dependent neighbouring squares
                    // for now they are added assuming it to be a big rectangle
                    // in the future these additions can be in form of 4 sectors for 4 corners 
                    // and 4 rectangles for top, bottom, left and right of the rectangle
                    add_neighbour_rectangle(gx*nx-eps, gy*ny-eps, gx*nx+nx+eps, gy*ny+ny+eps, S[current], vector_deps);

                    // represent dependencies using hpx dataflow
                    S[next][get_loc(gx, gy, np)] = hpx::dataflow(hpx::launch::async, Op, S[current][get_loc(gx, gy, np)], vector_deps, t, test);
                }
            }

            // every nd time steps, attach additional continuation which will
            // trigger the semaphore once computation has reached this point
            if ((t % nd) == 0)
            {
                S[next][0].then(
                    [&sem, t](partition_fut &&)
                    {
                        // inform semaphore about new lower limit
                        sem.signal(t);
                    });
            }

            // suspend if the tree has become too deep, the continuation above
            // will resume this thread once the computation has caught up
            sem.wait(t);

            // launch a seperate thread for collecting the logging data from 
            // various threads and output in the required files
            if(t % nlog == 0)
                hpx::dataflow(hpx::launch::async, 
                              hpx::util::unwrapping(&solver::log_csv_vtk),
                              S[next],
                              P,
                              t,
                              test);
        }

        // required for correctness for the case when nt = 0
        // because in that case next won't be set in the above for loop
        next = nt % 2;

        // wait for all the grids to be computed before we computing the errors
        hpx::wait_all(S[next]);

        //testing the code for correctness
        if(test)
        {
            compute_l2(nt);
            compute_linf(nt);
        }
    }
};

// function to execute the ctests which are there in file '2d_async.txt' for this file
int batch_tester(long nd, long nlog)
{
    std::uint64_t nt, num_tests;
    bool test_failed = 0;
    
    std::cin >> num_tests;
    for(long i = 0; i < num_tests; ++i)
    {
        std::cin >> nx >> ny >> np >> nt >> eps >> k >> dt >> dh;
        
        // Create the solver object
        solver solve(nt, nd, nlog, 1);
        
        solve.do_work();
        
        if (solve.error_l2 / (double)(nx * ny * np * np) > 1e-6)
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
    std::uint64_t nt = vm["nt"].as<std::uint64_t>();   // Number of steps.
    std::uint64_t nd = vm["nd"].as<std::uint64_t>();   // Depth for which dependency tree is allowed to grow
    std::uint64_t nlog = vm["nlog"].as<std::uint64_t>();   // Number of time steps to log the results
    bool test = vm["test"].as<bool>();

    if (vm.count("no-header"))
        header = false;

    //batch testing for ctesting
    if (vm.count("test_batch"))
        return batch_tester(nd, nlog);

    // Create the solver object
    solver solve(nt, nd, nlog, test);

    // Measure execution time.
    std::uint64_t t = hpx::util::high_resolution_clock::now();

    // Execute nt time steps on nx grid points.
    solve.do_work();

    std::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;

    //print the calculated error
    //print comparison only when the 'cmp' flag is set
    if(test)
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
        ("test", po::value<bool>()->default_value(true),
         "use arguments from numerical solution for testing (default: false)")
        ("test_batch", "test the solution against numerical solution against batch inputs (default: false)")
        ("results", "print generated results (default: false)")
        ("cmp", po::value<bool>()->default_value(false),
         "Print expected versus actual outputs")
        ("nx", po::value<long>(&nx)->default_value(25),
         "Local x dimension")
        ("ny", po::value<long>(&ny)->default_value(25),
        "Local y dimension")
        ("nt", po::value<std::uint64_t>()->default_value(45),
         "Number of time steps")
        ("nd", po::value<std::uint64_t>()->default_value(5),
         "Number of time steps to allow the dependency tree to grow to")
        ("np", po::value<long>(&np)->default_value(2),
         "Number of partitions in x and y dimension")
        ("nlog", po::value<std::uint64_t>()->default_value(5),
         "Number of time steps to log the results")
        ("eps", po::value<long>(&eps)->default_value(5),
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