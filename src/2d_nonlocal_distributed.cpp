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
#include <hpx/serialization/serialize.hpp>
#include <hpx/type_support/unused.hpp>

#include <boost/shared_array.hpp>

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <valarray>
#include <algorithm>
#include <fstream>

#include "../include/point.h"
#include "../include/print_time_results.hpp"
#include "../include/writer.h"

bool header = 1;
double k = 1 ;      // heat transfer coefficient
double dt = 1.;     // time step
double dh = 1.;     // grid spacing
long nx = 25;       // Number of grid points in x dimension
long ny = 25;       // Number of grid points in y dimension
long np = 5;        // Number of partitions we want to break the grid into
long nl = 1;        // Number of localities available for distributing the data
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

// function to get the localoty id of a particular partition
inline std::size_t locidx(std::size_t i, std::size_t np, std::size_t nl)
{
    return (i * nl) / (np * np);
}

// A small square partition of the original space in which we have divided the large mesh into
class partition_space
{
private:
    typedef hpx::serialization::serialize_buffer<double> buffer_type;

public:

    // coordinates of the grid point corresponding the np * np squares in the partition
    long gx, gy;

    explicit partition_space()
      : data_(std::allocator<double>().allocate(0), 0, buffer_type::take), 
        size_(0), 
        gx(0), gy(0)
    {}

    partition_space(std::size_t size, long gx, long gy)
      : data_(std::allocator<double>().allocate(size), size, buffer_type::take), 
        size_(size), 
        gx(gx), gy(gy)
    {}

    // add test initialization here for executing tests
    partition_space(long nx, long ny, long gx, long gy)
      : data_(std::allocator<double>().allocate(nx*ny), nx*ny, buffer_type::take),
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
    // serialization support for transfering objects across localities
    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & data_ & size_ & gx & gy;
    }

    buffer_type data_;
    std::size_t size_;
};

// server side of the data to access the component using it's global address
struct partition_space_server
  : hpx::components::component_base<partition_space_server>
{
    // construct new instances
    partition_space_server() {}

    explicit partition_space_server(partition_space const& data)
      : data_(data)
    {}

    partition_space_server(std::size_t size, long gx, long gy)
      : data_(size, gx, gy)
    {}

    partition_space_server(long nx, long ny, long gx, long gy)
      : data_(nx, ny, gx, gy)
    {}

    // Access data. The parameter specifies what part of the data should be
    // accessed. As long as the result is used locally, no data is copied,
    // however as soon as the result is requested from another locality only
    // then the data will be copied
    partition_space get_data() const
    {
        return data_;
    }

    // Every member function which has to be invoked remotely needs to be
    // wrapped into a component action. The macro below defines a new type
    // 'get_data_action' which represents the (possibly remote) member function
    // partition::get_data().
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_space_server, get_data, get_data_action);

private:
    partition_space data_;
};

// HPX_REGISTER_COMPONENT() exposes the component creation
// through hpx::new_<>().
typedef hpx::components::component<partition_space_server> partition_space_server_type;
HPX_REGISTER_COMPONENT(partition_space_server_type, partition_space_server);

// HPX_REGISTER_ACTION() exposes the component member function for remote
// invocation.
typedef partition_space_server::get_data_action get_data_action;
HPX_REGISTER_ACTION(get_data_action);


// client side code for accessing the partition server from possibly a different locality
struct partition_space_client 
    : hpx::components::client_base<partition_space_client, partition_space_server>
{
    typedef hpx::components::client_base<partition_space_client, partition_space_server> base_type;

    partition_space_client() {}

    // Create new component on locality 'where' and initialize the held data
    partition_space_client(hpx::id_type where, std::size_t size, long gx, long gy)
      : base_type(hpx::new_<partition_space_server>(where, size, gx, gy))
    {}

    // Create new component on locality 'where' and initialize the held data
    partition_space_client(hpx::id_type where, long nx, long ny, long gx, long gy)
      : base_type(hpx::new_<partition_space_server>(where, nx, ny, gx, gy))
    {}

    // Create a new component on the locality co-located to the id 'where'. The
    // new instance will be initialized from the given partition_data.
    partition_space_client(hpx::id_type where, partition_space const& data)
      : base_type(hpx::new_<partition_space_server>(hpx::colocated(where), data))
    {}

    // Attach a future representing a (possibly remote) partition.
    explicit partition_space_client(hpx::future<hpx::id_type> && id)
      : base_type(std::move(id))
    {}

    // Unwrap a future<partition> (a partition already holds a future to the
    // id of the referenced object, thus unwrapping accesses this inner future).
    partition_space_client(hpx::future<partition_space_client> && c)
      : base_type(std::move(c))
    {}

    ///////////////////////////////////////////////////////////////////////////
    // Invoke the (remote) member function which gives us access to the data.
    // This is a pure helper function hiding the async.
    hpx::future<partition_space> get_data() const
    {
        partition_space_server::get_data_action act;
        return hpx::async(act, get_id());
    }
};

class solver
{
public:
    // 3d coordinate for representation of point
    typedef util::Point3 point_3d;

    // 3d plane using the 3d point representation
    // Note that the plane is actually stored in a 1d fashion
    typedef std::vector<point_3d> plane_3d;

    // partitioning the 2d space into np * np small squares which constitute a partition of space
    typedef hpx::shared_future<partition_space> partition_fut;

    // 2d space data structure which consists of future of square partitions
    typedef std::vector<partition_fut> space_2d_fut;

    // 2d space representation for temperature data
    typedef std::vector<partition_space_client> space_2d;

    //alternate for space between time steps,
    //result of S[t%2] goes to S[(t+1)%2] at every timestep 't'
    space_2d S[2];

    // vector containining 3d coordinates of the points in the plane
    plane_3d P;

    // nt = number of timesteps
    // nlog = Number of time steps to log the results
    long nt, nlog;
    bool current, next, test;
    // l2 norm and l infinity norm
    double error_l2, error_linf;

    // vector of id's of various available localitites for computation
    std::vector<hpx::id_type> localities;

    //constructor to initialize class variables and allocate
    solver(long nt, long nlog, bool test)
    {
        this->nt = nt;
        this->nlog = nlog;
        this->current = 0;
        this->error_l2 = 0.0;
        this->error_linf = 0.0;
        this->test = test;
        this->next = 1;

        // forming the vector of localities id's and assigning the 'nl'
        std::vector<hpx::id_type> all_localities = hpx::find_all_localities();

        // checking if the number of partitions is more than the number of localities
        for(long idx = 0; idx < std::min((long)all_localities.size(), np*np); ++idx)
        {
            this->localities.push_back(all_localities[idx]);
        }
        nl = this->localities.size();

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

        for(space_2d &s : S)
            s.resize(np * np);
        
        // initial conditions for testing the correctness
        for (long i = 0; i < np*np; ++i)
            S[0][i] = partition_space_client(localities[locidx(i, np, nl)], nx, ny, i % np, i / np);
    }

    // function to compute l2 norm of the solution
    void compute_l2(long time)
    {
        auto solution = hpx::util::unwrap(vector_get_data(S[next]));
        error_l2 = 0;

        for(long sx = 0; sx < nx * np; ++sx)
            for(long sy = 0; sy < ny * np; ++sy)
                error_l2 += std::pow(solution[grid_loc(sx, sy)][mesh_loc(sx, sy)] - w(sx, sy, time), 2);
    }

    // function to compute l infinity norm of the solution
    void compute_linf(long time)
    {
        auto solution = hpx::util::unwrap(vector_get_data(S[next]));
        error_linf = 0;
        
        for(long sx = 0; sx < nx * np; ++sx)
            for(long sy = 0; sy < ny * np; ++sy)
                error_linf = std::max(std::abs(solution[grid_loc(sx, sy)][mesh_loc(sx, sy)] - w(sx, sy, time)), error_linf);
    }

    //print error for testing
    void print_error(bool cmp)
    {        
        compute_l2(nt);
        compute_linf(nt);
        
        std::cout << "l2: " << error_l2 << " linfinity: " << error_linf << std::endl;

        if(cmp)
        {
            auto solution = hpx::util::unwrap(vector_get_data(S[next]));
            for(long sx = 0; sx < nx * np; ++sx)
            {
                for(long sy = 0; sy < ny * np; ++sy)
                {
                    std::cout << "sx: " << sx << " sy: " << sy << " Expected: " << w(sx, sy, nt)
                    << " Actual: " << solution[grid_loc(sx, sy)][mesh_loc(sx, sy)] 
                    << std::endl;
                }
                std::cout << std::endl;
            }
        }
    }

    //print the solution for the user
    void print_soln()
    {
        auto solution = hpx::util::unwrap(vector_get_data(S[next]));

        for (long sx = 0; sx < nx*np; ++sx)
        {
            for (long sy = 0; sy < ny*np; ++sy)
            {
                std::cout << "S[" << sx << "][" << sy << "] = " 
                << solution[grid_loc(sx, sy)][mesh_loc(sx, sy)] << " ";
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
    static void add_neighbour_rectangle(long lx, long ly, long rx, long ry, space_2d &all_squares
                       , space_2d &neighbour_squares)
    {
        for(long sx = lx; sx < rx + nx; sx += nx)
        {
            for(long sy = ly; sy < ry + ny; sy += ny)
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
    static inline double boundary(long pos_x, long pos_y, const std::vector<partition_space> &all_squares)
    {
        if(pos_x >= 0 && pos_x < nx * np && pos_y >= 0 && pos_y < ny * np)
            return all_squares[grid_loc(pos_x, pos_y)][mesh_loc(pos_x, pos_y)];
        else
            return 0;
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
        long len_line = 0;
        double w_position = w(pos_x, pos_y, time);

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
    static double sum_local(long pos_x, long pos_y, const std::vector<partition_space> &all_squares)
    {
        double result_local = 0.0;
        long len_line = 0;
        double pos_val = boundary(pos_x, pos_y, all_squares);

        for(long sx = pos_x-eps; sx <= pos_x+eps; ++sx)
        {
            len_line = len_1d_line(std::abs((long)pos_x - (long)sx));
            for(long sy = pos_y-len_line; sy <= pos_y+len_line; ++sy)
            {
                result_local += influence_function(distance(pos_x, pos_y, sx, sy))
                                * c_2d
                                * (boundary(sx, sy, all_squares) - pos_val)
                                * (dh * dh);
            }
        }
        
        return result_local;
    }

    static space_2d_fut vector_get_data(const std::vector<partition_space_client> &neighbour_squares)
    {
        space_2d_fut ret_get_data;
        ret_get_data.resize(neighbour_squares.size());

        for(long idx = 0; idx < neighbour_squares.size(); ++idx)
        {
            ret_get_data[idx] = neighbour_squares[idx].get_data();
        }

        return ret_get_data;
    }

    // partition of nx * ny cells which will be processed by a single thread
    static partition_space_client sum_local_partition(partition_space_client const& middle, 
                                               const std::vector<partition_space_client> &neighbour_squares, 
                                               long time, bool test)
    {                 
        // all_squares is a vector with partitions with pointers to valid partitions
        // to only those neighbours who affect the temperature at the middle partition
        // std::vector<partition_space> all_squares;
        // all_squares.resize(np * np);

        hpx::shared_future<partition_space> middle_data = middle.get_data();

        hpx::future<partition_space> next_middle = middle_data.then(
            hpx::util::unwrapping(
                [time, test](partition_space const& mid) -> partition_space
                {
                    std::vector<partition_space> all_squares;
                    all_squares.resize(np * np);
                    // All local operations are performed once the middle data of
                    // the previous time step becomes available.
                    long gx = mid.gx, gy = mid.gy;
                    partition_space next(mid.size(), gx, gy);
                    all_squares[get_loc(gx, gy, np)] = mid;
                    
                    for(long sx = eps+1; sx < nx-eps-1; ++sx)
                    {
                        for(long sy = eps+1; sy < ny-eps-1; ++sy)
                        {
                            next[get_loc(sx, sy)] = mid[get_loc(sx, sy)] + (sum_local(sx + gx*nx, sy + gy*ny, all_squares) * dt);
                            if(test)
                                next[get_loc(sx, sy)] += sum_local_test(sx + gx*nx, sy + gy*ny, time) * dt;
                        }
                    }
                    return next;
                }
            )
        );

        return hpx::dataflow(
            hpx::launch::async,
            hpx::util::unwrapping(
                [middle, time, test]( partition_space next, 
                                                    partition_space const& mid,
                                                    const std::vector<partition_space> &neighbour_squares
                                                    ) -> partition_space_client
                {
                    // Calculate the missing boundary elements once the
                    // corresponding data has become available.
                    std::vector<partition_space> all_squares;
                    all_squares.resize(np * np);
                    long gx = mid.gx, gy = mid.gy;

                    // assign only those partition pointers which are valid neighbours for a particular square in the grid
                    for(long idx = 0; idx < neighbour_squares.size(); ++idx)
                    {
                        all_squares[get_loc(neighbour_squares[idx].gx, neighbour_squares[idx].gy, np)] = neighbour_squares[idx];
                    }

                    if(nx <= eps)
                    {
                        for(long sx = 0; sx < nx; ++sx)
                        {
                            for(long sy = 0; sy < ny; ++sy)
                            {
                                next[get_loc(sx, sy)] = mid[get_loc(sx, sy)] + (sum_local(sx + gx*nx, sy + gy*ny, all_squares) * dt);
                                if(test)
                                    next[get_loc(sx, sy)] += sum_local_test(sx + gx*nx, sy + gy*ny, time) * dt;
                            }
                        }
                    }
                    else
                    {
                        for(long sx = 0; sx < eps+1; ++sx)
                        {
                            for(long sy = 0; sy < ny; ++sy)
                            {
                                next[get_loc(sx, sy)] = mid[get_loc(sx, sy)] + (sum_local(sx + gx*nx, sy + gy*ny, all_squares) * dt);
                                if(test)
                                    next[get_loc(sx, sy)] += sum_local_test(sx + gx*nx, sy + gy*ny, time) * dt;
                            }
                        }

                        for(long sx = eps + 1; sx < nx-eps-1; ++sx)
                        {
                            for(long sy = 0; sy < eps+1; ++sy)
                            {
                                next[get_loc(sx, sy)] = mid[get_loc(sx, sy)] + (sum_local(sx + gx*nx, sy + gy*ny, all_squares) * dt);
                                if(test)
                                    next[get_loc(sx, sy)] += sum_local_test(sx + gx*nx, sy + gy*ny, time) * dt;
                            }
                            for(long sy = ny-eps-1; sy < ny; ++sy)
                            {
                                next[get_loc(sx, sy)] = mid[get_loc(sx, sy)] + (sum_local(sx + gx*nx, sy + gy*ny, all_squares) * dt);
                                if(test)
                                    next[get_loc(sx, sy)] += sum_local_test(sx + gx*nx, sy + gy*ny, time) * dt;
                            }
                        }

                        for(long sx = nx-eps-1; sx < nx; ++sx)
                        {
                            for(long sy = 0; sy < ny; ++sy)
                            {
                                next[get_loc(sx, sy)] = mid[get_loc(sx, sy)] + (sum_local(sx + gx*nx, sy + gy*ny, all_squares) * dt);
                                if(test)
                                    next[get_loc(sx, sy)] += sum_local_test(sx + gx*nx, sy + gy*ny, time) * dt;
                            }
                        }
                    }

                    // The new partition_data will be allocated on the same locality
                    // as 'middle'.
                    return partition_space_client(middle.get_id(), next);
                }
            ),
            std::move(next_middle),
            middle_data,
            vector_get_data(neighbour_squares)
        );
    }

    void do_work();
};

HPX_PLAIN_ACTION(solver::sum_local_partition, sum_local_partition_action);

// do all the work on 'nx * ny' data points for each of 'np * np' partitions for 'nt' time steps
void solver::do_work()
{
    // Actual time step loop
    for (long t = 0; t < nt; ++t)
    {
        current = t % 2;
        next = (t + 1) % 2;

        sum_local_partition_action act;
        for(long gx = 0; gx < np; ++gx)
        {
            for(long gy = 0; gy < np; ++gy)
            {
                // we execute the action on the locality of the middle partition
                auto Op = hpx::util::bind_front(act, localities[locidx(get_loc(gx, gy, np), np, nl)]);

                // vector of dependent grid squares for a particular grid squares
                space_2d vector_deps;
                
                // add dependent neighbouring squares
                // for now they are added assuming it to be a big rectangle
                // in the future these additions can be in form of 4 sectors for 4 corners 
                // and 4 rectangles for top, bottom, left and right of the rectangle
                add_neighbour_rectangle(gx*nx - eps, gy*ny - eps, gx*nx+nx+eps, gy*ny+ny+eps, S[current], vector_deps);

                // represent dependencies using hpx dataflow
                S[next][get_loc(gx, gy, np)] = hpx::dataflow(hpx::launch::async, Op, S[current][get_loc(gx, gy, np)], vector_deps, t, test);
            }
        }

        // launch a seperate thread for collecting the logging data from 
        // various threads and output in the required files
        if(t % nlog == 0)
            hpx::dataflow(hpx::launch::async, 
                            hpx::util::unwrapping(&solver::log_csv_vtk),
                            vector_get_data(S[next]),
                            P,
                            t,
                            test);
    }

    // required for correctness for the case when nt = 0
    // because in that case next won't be set in the above for loop
    next = nt % 2;

    // wait for the solution to be ready for all the partitions
    for(int idx = 0; idx < np * np; ++idx)
        S[next][idx].get_data().wait();
}

// function to execute the ctests which are there in file '2d_async.txt' for this file
int batch_tester(long nlog)
{
    std::uint64_t nt, num_tests;
    bool test_failed = 0;
    
    std::cin >> num_tests;
    for(long i = 0; i < num_tests; ++i)
    {
        std::cin >> nx >> ny >> np >> nt >> eps >> k >> dt >> dh;
        
        // Create the solver object
        solver solve(nt, nlog, 1);
        
        // do the actual error
        solve.do_work();

        // compute and print the error for checking the difference between analytical and numerical solution
        solve.compute_l2(nt);
        
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
    std::uint64_t nt = vm["nt"].as<std::uint64_t>();        // Number of steps.
    std::uint64_t nlog = vm["nlog"].as<std::uint64_t>();    // Number of time steps to log the results
    bool test = vm["test"].as<bool>();                      // Boolean variable to indicate if numerical solution needs to be tested
    if(nx <= eps)
        std::cout << "[WARNING] Mesh size on a single node (nx * ny) is too small for given epsilon (eps)" << std::endl;

    if (vm.count("no-header"))
        header = false;

    //batch testing for ctesting
    if (vm.count("test_batch"))
        return batch_tester(nlog);

    // Create the solver object
    solver solve(nt, nlog, test);

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

    std::uint64_t const num_worker_threads = hpx::get_num_worker_threads();
    hpx::future<std::uint32_t> locs = hpx::get_num_localities();
    print_time_results(
        locs.get(), num_worker_threads, elapsed, nx, ny, np, nt, header);

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