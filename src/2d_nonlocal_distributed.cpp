//  Copyright (c) 2014 Hartmut Kaiser
//  Copyright (c) 2014 Patricia Grubel
//  Copyright (c) 2020 Pranav Gadikar
//
//  SPDX-License-Identifier: BSL-1.0
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Note on terminology  :-
// grid square  : small squares/ rectangles into which we have divided the
// larger mesh mesh points  : points within a single grid square which is
// smallest unit of space in our equation

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/serialization/serialize.hpp>
#include <hpx/type_support/unused.hpp>
#include <hpx/include/performance_counters.hpp>

#include <boost/shared_array.hpp>

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>
#include <sys/time.h>
#include <unordered_map>
#include <valarray>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <queue>
#include <utility>
#include <set>

#include "../include/point.h"
#include "../include/print_time_results.hpp"
#include "../include/writer.h"
#include "../include/Config.h"

#define ff first
#define ss second

bool header = 1;    //!< print csv heading
double k = 1;       //!< heat transfer coefficient
double dt = 1.;     //!< time step
double dh = 0.005;  //!< grid spacing
long nx = 20;       //!< Number of grid points in x dimension
long ny = 20;       //!< Number of grid points in y dimension
long npx =
    10;  //!< Number of partitions we want to break the grid in x-dimension
long npy =
    10;       //!< Number of partitions we want to break the grid in y-dimension
long nl = 1;  //!< Number of localities available for distributing the data
long eps = 5;       //!< Epsilon for influence zone of a point
double c_2d = 0.0;  //!< 'c' to be used in nonlocal equation for 2d case when
                    //!< J(influence function) = 1
std::vector<int> partition_domain;  //!< std::vector to store locality that
                                    //!< holds each partition
std::vector<hpx::performance_counters::performance_counter>
    perf_counters;  //!< std::vector of performance counters

//! operator for getting the location of 2d point in 1d representation
/**
 * @param x x coordinate in 2d cartesian coordinates
 * @param y y coordinate in 2d cartesian coordinates
 * @param nx length of the 2d cartesian plane in x dimension
 * @return location of 2d point in 1d array used to store the 2d mesh
 */
inline long get_loc(long x, long y, long nx = nx) { return x + y * nx; }

//! operator for getting 2d coordinate of grid square in which the space has
//! been divided
/**
 * Internally calls get_loc() function
 * @see get_loc()
 * @param x x coordinate in 2d cartesian coordinates
 * @param y y coordinate in 2d cartesian coordinates
 * @return location of 2d grid's square in 1d array used to store the 2d grid
 */
inline long grid_loc(long x, long y) { return get_loc(x / nx, y / ny, npx); }

//! operator for getting 2d coordinate of a point within the grid square
/**
 * Internally calls get_loc() function
 * @see get_loc()
 * @param x x coordinate in 2d cartesian coordinates
 * @param y y coordinate in 2d cartesian coordinates
 * @return location of 2d point within a grid square in 1d array used to store
 * the square subdomain
 */
inline long mesh_loc(long x, long y) { return get_loc(x % nx, y % ny); }

//! function to get the locality id of a particular partition
/**
 * @see partition_domain
 * @param i location of 2d point in 1d array used to store individual square
 * subdomains
 * @param nl number of localities avaible for distributing the computation
 * @return index of locality in the localities array
 */
inline std::size_t locidx(std::size_t i, std::size_t nl) {
  if (partition_domain.size() == 0)
    return (i * nl) / (npx * npy);
  else
    return partition_domain[i];
}

void setup_counters() {
  // We need to explicitly start all counters before we can use them. For
  // certain counters this could be a no-op, in which case start will return
  // 'false'.
  int num_loc = hpx::get_num_localities().get();

  // cannot test when the number of nodes is 1
  // because idle-rate performance counter is not available
  if (num_loc == 1) return;
  perf_counters.resize(num_loc);

  for (int idx = 0; idx < num_loc; ++idx) {
    perf_counters[idx] = hpx::performance_counters::performance_counter(
        "/threads{locality#" + std::to_string(idx) + "/total}/idle-rate");
    perf_counters[idx].start(hpx::launch::sync);
  }
}

/**
 * A small square partition of the original space in which we have divided the
 * large mesh into
 */
class partition_space {
 private:
  typedef hpx::serialization::serialize_buffer<double> buffer_type;
  //!< Data type of the serialized buffer

 public:
  long gx;  //!< x coordinate of the grid point corresponding the npx * npy
            //!< squares
  long gy;  //!< y coordinate of the grid point corresponding the npx * npy
            //!< squares

  /**
   * Empty constructor for initializing the partition_space class
   */
  explicit partition_space()
      : data_(std::allocator<double>().allocate(0), 0, buffer_type::take),
        size_(0),
        gx(0),
        gy(0) {}

  /**
   * Parameterized constructor for initializing the partition_space object
   * @param size size of the data buffer to be created in partition_space
   * @param gx x coordinate of the grid point corresponding the npx * npy
   * squares
   * @param gy y coordinate of the grid point corresponding the npx * npy
   * squares
   */
  partition_space(std::size_t size, long gx, long gy)
      : data_(std::allocator<double>().allocate(size), size, buffer_type::take),
        size_(size),
        gx(gx),
        gy(gy) {}

  //! function to add test initialization for executing tests
  /**
   * Parameterized constructor for initializing the partition_space object
   * @param nx size of square domain in x dimension
   * @param ny size of square domain in y dimension
   * @param gx x coordinate of the grid point corresponding the npx * npy
   * squares
   * @param gy y coordinate of the grid point corresponding the npx * npy
   * squares
   */
  partition_space(long nx, long ny, long gx, long gy)
      : data_(std::allocator<double>().allocate(nx * ny), nx * ny,
              buffer_type::take),
        size_(nx * ny),
        gx(gx),
        gy(gy) {
    for (long sx = gx * nx; sx < nx + gx * nx; ++sx) {
      for (long sy = gy * ny; sy < ny + gy * ny; ++sy) {
        data_[mesh_loc(sx, sy)] =
            sin(2 * M_PI * (sx * dh)) * sin(2 * M_PI * (sy * dh));
      }
    }
  }

  //! operator to return the data at some index from the given partition
  double& operator[](std::size_t idx) { return data_[idx]; }
  //! operator to return the data at some index from the given partition
  double operator[](std::size_t idx) const { return data_[idx]; }

  //! operator to return the size of the data partition in partition_space
  //! object
  std::size_t size() const { return size_; }

 private:
  friend class hpx::serialization::access;

  //! serialization support for transfering objects across localities
  template <typename Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& data_& size_& gx& gy;
  }

  buffer_type
      data_;  //!< buffer to store the temperature values for a square subdomain
  std::size_t size_;  //!< size of the buffer storing the square subdomain
};

/**
 * server side of the data to access the component using it's global address
 */
struct partition_space_server
    : hpx::components::component_base<partition_space_server> {
  partition_space_server() {}

  /**
   * Parameterized constructor to create a partition_space_server object
   * for given temperature values in a square subdomain
   * @param data partition_space object storing temperature values in a square
   * subdomain
   */
  explicit partition_space_server(partition_space const& data) : data_(data) {}

  /**
   * Parameterized constructor for initializing the partition_space_server
   * object
   * @param size size of the data buffer to be created in partition_space
   * @param gx x coordinate of the grid point corresponding the npx * npy
   * squares
   * @param gy y coordinate of the grid point corresponding the npx * npy
   * squares
   */
  partition_space_server(std::size_t size, long gx, long gy)
      : data_(size, gx, gy) {}

  /**
   * Parameterized constructor for initializing the partition_space_server
   * object
   * @param nx size of square domain in x dimension
   * @param ny size of square domain in y dimension
   * @param gx x coordinate of the grid point corresponding the npx * npy
   * squares
   * @param gy y coordinate of the grid point corresponding the npx * npy
   * squares
   */
  partition_space_server(long nx, long ny, long gx, long gy)
      : data_(nx, ny, gx, gy) {}

  //! Access data. The parameter specifies what part of the data should be
  //! accessed. As long as the result is used locally, no data is copied,
  //! however as soon as the result is requested from another locality only
  //! then the data will be copied
  partition_space get_data() const { return data_; }

  //! Every member function which has to be invoked remotely needs to be
  //! wrapped into a component action. The macro below defines a new type
  //! 'get_data_action' which represents the (possibly remote) member function
  //! partition::get_data().
  HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_space_server, get_data,
                                     get_data_action);

 private:
  partition_space data_;  //!< partition_space object storing temperature values
                          //!< which is exposed by the server
};

//! HPX_REGISTER_COMPONENT() exposes the component creation
//! through hpx::new_<>().
typedef hpx::components::component<partition_space_server>
    partition_space_server_type;
HPX_REGISTER_COMPONENT(partition_space_server_type, partition_space_server);

//! HPX_REGISTER_ACTION() exposes the component member function for remote
//! invocation.
typedef partition_space_server::get_data_action get_data_action;
HPX_REGISTER_ACTION(get_data_action);

/**
 * client side code for accessing the partition server from possibly a different
 * locality
 */
struct partition_space_client
    : hpx::components::client_base<partition_space_client,
                                   partition_space_server> {
  typedef hpx::components::client_base<partition_space_client,
                                       partition_space_server>
      base_type;  //!< client base type for the partition_space object

  partition_space_client() {}

  //! Create new component on locality 'where' and initialize the held data
  /**
   * Parameterized constructor for initializing the partition_space_server
   * object
   * @param where locality where the held data is initialized
   * @param size size of the data buffer to be created in partition_space
   * @param gx x coordinate of the grid point corresponding the npx * npy
   * squares
   * @param gy y coordinate of the grid point corresponding the npx * npy
   * squares
   */
  partition_space_client(hpx::id_type where, std::size_t size, long gx, long gy)
      : base_type(hpx::new_<partition_space_server>(where, size, gx, gy)) {}

  //! Create new component on locality 'where' and initialize the held data
  /**
   * Parameterized constructor for initializing the partition_space_server
   * object
   * @param where locality where the held data is initialized
   * @param nx size of square domain in x dimension
   * @param ny size of square domain in y dimension
   * @param gx x coordinate of the grid point corresponding the npx * npy
   * squares
   * @param gy y coordinate of the grid point corresponding the npx * npy
   * squares
   */
  partition_space_client(hpx::id_type where, long nx, long ny, long gx, long gy)
      : base_type(hpx::new_<partition_space_server>(where, nx, ny, gx, gy)) {}

  /**
   * Create a new component on the locality co-located to the id 'where'. The
   * new instance will be initialized from the given partition_data
   * @param where locality where the held data is initialized
   * @param data partition_space object storing temperature values in a square
   * subdomain
   */
  partition_space_client(hpx::id_type where, partition_space const& data)
      : base_type(
            hpx::new_<partition_space_server>(hpx::colocated(where), data)) {}

  //! Attach a future representing a (possibly remote) partition.
  explicit partition_space_client(hpx::future<hpx::id_type>&& id)
      : base_type(std::move(id)) {}

  //! Unwrap a future<partition> (a partition already holds a future to the
  //! id of the referenced object, thus unwrapping accesses this inner future).
  partition_space_client(hpx::future<partition_space_client>&& c)
      : base_type(std::move(c)) {}

  //! Invoke the (remote) member function which gives us access to the data.
  //! This is a pure helper function hiding the async.
  hpx::future<partition_space> get_data() const {
    partition_space_server::get_data_action act;
    return hpx::async(act, get_id());
  }
};

/**
 * class which contains the code to solve the 2d nonlocal equation
 * in a distributed setup. The class contains functions to run the
 * solver in a distributed setup with support for logging at regular
 * intervals.
 */
class solver {
 public:
  //! 3d coordinate for representation of point
  typedef util::Point3 point_3d;

  //! 3d plane using the 3d point representation
  //! Note that the plane is actually stored in a 1d fashion
  typedef std::vector<point_3d> plane_3d;

  //! partitioning the 2d space into npx * npy small squares which constitute a
  //! partition of space
  typedef hpx::shared_future<partition_space> partition_fut;

  //! 2d space data structure which consists of future of square partitions
  typedef std::vector<partition_fut> space_2d_fut;

  //! 2d space representation for temperature data
  typedef std::vector<partition_space_client> space_2d;

  space_2d S[2];  //!< alternate for space between time steps,
  //!< result of S[t%2] goes to S[(t+1)%2] at every timestep 't'

  plane_3d
      P;  //!< vector containining 3d coordinates of the points in the plane

  long nt;        //!< number of timesteps
  long nlog;      //!< Number of time steps to log the results
  long nbalance;  //!< Number of time steps to run load balancing algorithm
  bool current;   //!< boolean variable to store index of current mesh which is
                  //!< already computed
  bool next;  //!< boolean variable to store index of next mesh which is being
              //!< computed
  bool test;  //!< true if the analyitcal solution needs to be compared with the
              //!< numerical solution

  double error_l2;    //!< l2 norm
  double error_linf;  //!< l infinity norm

  std::vector<hpx::id_type>
      localities;  //!< std::vector of id's of various available localitites for
                   //!< computation

  //! constructor to allocate and initialize class variables
  /**
   * Internally calls get_loc() function
   * @see partition_space_client
   * @see partition_space_server
   * @see partition_space
   * @param nt number of timesteps to execute the sover
   * @param nlog number of iterations after which we want to log
   * @param nbalance number of time steps to run load balancing algorithm
   * @param test true if numerical solution needs to be tested against
   * analytical solution
   * @param filename name of file to input the mesh from
   */
  solver(long nt, long nlog, long nbalance, bool test,
         const std::string& filename = "None") {
    this->nt = nt;
    this->nlog = nlog;
    this->nbalance = nbalance;
    this->current = 0;
    this->error_l2 = 0.0;
    this->error_linf = 0.0;
    this->test = test;
    this->next = 1;

    // take parameters input from a file
    if (filename != "None") param_file_input(filename);

    // forming the vector of localities id's and assigning the 'nl'
    std::vector<hpx::id_type> all_localities = hpx::find_all_localities();

    this->localities.resize(all_localities.size());

    // checking if the number of partitions is less than the number of
    // localities
    for (long idx = 0; idx < std::min((long)all_localities.size(), npx * npy);
         ++idx) {
      this->localities[hpx::naming::get_locality_id_from_id(
          all_localities[idx])] = all_localities[idx];
    }
    nl = this->localities.size();

    // setting the global variable c_2d which is a constant for this test case
    // as per the inputs
    c_2d = (k * 8) / pow(eps * dh, 4);

    P.resize(nx * ny * npx * npy);
    for (long sx = 0; sx < nx * npx; ++sx) {
      for (long sy = 0; sy < ny * npy; ++sy) {
        P[sx + sy * (nx * npx)] = point_3d(sx, sy, 0);
      }
    }
    nl = this->localities.size();

    for (space_2d& s : S) s.resize(npx * npy);

    // initial conditions for testing the correctness
    for (long i = 0; i < npx * npy; ++i)
      S[0][i] = partition_space_client(this->localities[locidx(i, nl)], nx, ny,
                                       i % npx, i / npx);
  }

  /**
   * Internally calls get_loc() function
   * @param filename name of file to input the mesh from
   */
  void param_file_input(const std::string& filename) {
    char char_array[filename.length() + 1];
    strcpy(char_array, filename.c_str());

    std::ifstream infile;
    infile.open(char_array);

    if (infile) {
      std::size_t px, py, loc;
      infile >> nx >> ny >> npx >> npy >> dh;
      partition_domain.resize(npx * npy);

      for (long idx = 0; idx < npx; ++idx) {
        for (long idy = 0; idy < npy; ++idy) {
          infile >> px >> py >> loc;
          partition_domain[get_loc(px, py, npx)] = loc;
        }
      }
    }

    infile.close();
  }

  /**
   * function to compute l2 norm of the solution
   * @param time time at which we want to check the analytical solution against
   * numerical solution
   */
  void compute_l2(long time) {
    auto solution = hpx::util::unwrap(vector_get_data(S[next]));
    error_l2 = 0;

    for (long sx = 0; sx < nx * npx; ++sx)
      for (long sy = 0; sy < ny * npy; ++sy)
        error_l2 += std::pow(
            solution[grid_loc(sx, sy)][mesh_loc(sx, sy)] - w(sx, sy, time), 2);
  }

  /**
   * function to compute l infinity norm of the solution
   * @param time time at which we want to check the analytical solution against
   * numerical solution
   */
  void compute_linf(long time) {
    auto solution = hpx::util::unwrap(vector_get_data(S[next]));
    error_linf = 0;

    for (long sx = 0; sx < nx * npx; ++sx)
      for (long sy = 0; sy < ny * npy; ++sy)
        error_linf =
            std::max(std::abs(solution[grid_loc(sx, sy)][mesh_loc(sx, sy)] -
                              w(sx, sy, time)),
                     error_linf);
  }

  /**
   * print error for testing
   * @param cmp true if we want the print the analytical solution comparison
   * against the numerical solution
   */
  void print_error(bool cmp) {
    compute_l2(nt);
    compute_linf(nt);

    std::cout << "l2: " << error_l2 << " linfinity: " << error_linf
              << std::endl;

    if (cmp) {
      auto solution = hpx::util::unwrap(vector_get_data(S[next]));
      for (long sx = 0; sx < nx * npx; ++sx) {
        for (long sy = 0; sy < ny * npy; ++sy) {
          std::cout << "sx: " << sx << " sy: " << sy
                    << " Expected: " << w(sx, sy, nt) << " Actual: "
                    << solution[grid_loc(sx, sy)][mesh_loc(sx, sy)]
                    << std::endl;
        }
        std::cout << std::endl;
      }
    }
  }

  //! print the solution for the user
  void print_soln() {
    auto solution = hpx::util::unwrap(vector_get_data(S[next]));

    for (long sx = 0; sx < nx * npx; ++sx) {
      for (long sy = 0; sy < ny * npy; ++sy) {
        std::cout << "S[" << sx << "][" << sy
                  << "] = " << solution[grid_loc(sx, sy)][mesh_loc(sx, sy)]
                  << " ";
      }
      std::cout << std::endl;
    }
  }

  /**
   * function to perform the logging in both csv and vtk format
   * @param temp_data temperature data for 2d square subdomain
   * @param coord 3d coordinates for points in Point3 format
   * @param time time for which we are logging the results
   * @param test variable is true if we have to test the numerical solution
   * against the analytical solution
   */
  static void log_csv_vtk(std::vector<partition_space> temp_data,
                          plane_3d coord, long time, bool test) {
    // file to store the simulation results in csv format
    const std::string simulate_fname = "../out_csv/simulate_2d.csv";

    // file to store l2 and l infinity norms per timestep
    const std::string score_fname = "../out_csv/score_2d.csv";

    std::ofstream outfile;
    outfile.open(simulate_fname, std::ios_base::app);

    for (long sx = 0; sx < nx * npx; ++sx) {
      for (long sy = 0; sy < ny * npy; ++sy) {
        outfile << time << "," << sx << "," << sy << ","
                << temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] << ","
                << w(sx, sy, time) << ","
                << std::pow((temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] -
                             w(sx, sy, time)),
                            2)
                << ","
                << std::abs(temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] -
                            w(sx, sy, time))
                << ",\n";
      }
    }

    outfile.close();

    // VTK writer begins here

    // vtu file name for writing the vtk data
    const std::string fname = "../out_vtk/simulate_" + std::to_string(time);
    rw::writer::VtkWriter vtk_logger(fname);
    std::vector<double> vector_data(nx * ny * npx * npy);

    vtk_logger.appendNodes(&coord);

    for (long sx = 0; sx < nx * npx; ++sx)
      for (long sy = 0; sy < ny * npy; ++sy)
        vector_data[sy * nx * npx + sx] =
            temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)];

    vtk_logger.appendPointData("Temperature", &vector_data);
    vtk_logger.addTimeStep(std::time(0));
    vtk_logger.close();

    // VTK writer ends here

    // add to the score csv only when it's required
    if (test) {
      double error_l2 = 0, error_linf = 0;

      for (long sx = 0; sx < nx * npx; ++sx) {
        for (long sy = 0; sy < ny * npy; ++sy) {
          error_linf =
              std::max(std::abs(temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] -
                                w(sx, sy, time)),
                       error_linf);
          error_l2 += std::pow(
              temp_data[grid_loc(sx, sy)][mesh_loc(sx, sy)] - w(sx, sy, time),
              2);
        }
      }

      std::ofstream outfile;

      outfile.open(score_fname, std::ios_base::app);
      outfile << time << "," << error_l2 << "," << error_linf << ",\n";
      outfile.close();
    }
  }

  /**
   * function for checking the load balancing and visualizing the same
   * using hpx::performance counters
   */
  void test_load_balance() {
    double busy_rates[nl];  // idle rate values returned by performance counters
    double expected_busy_rate =
        0.0;  // hypotheticallly expected idle rate for equal distribution of
              // work among all nodes
    double max_diff = 0.0;  // variable to store the maximum difference between
                            // expected and actual idle rate
    std::cout << "Testing load balance:" << std::endl;

    // retrieve the performance counter values for idle-rates
    for (int idx = 0; idx < nl; ++idx) {
      hpx::performance_counters::counter_value value1 =
          perf_counters[idx].get_counter_value(hpx::launch::sync, true);
      busy_rates[idx] = 10000 - value1.get_value<double>();
      std::cout << "Test: counter value: " << busy_rates[idx] << std::endl;

      expected_busy_rate += busy_rates[idx];
    }
    expected_busy_rate /= nl;
    std::cout << "Expected busy rate " << expected_busy_rate << std::endl;

    for (int idx = 0; idx < nl; ++idx) {
      max_diff =
          std::max(std::abs(expected_busy_rate - busy_rates[idx]), max_diff);
    }

    std::cout << "Visualizing Load Balance across nodes" << std::endl;

    for (long idx = 0; idx < npx; ++idx) {
      for (long idy = 0; idy < npy; ++idy) {
        std::cout << partition_domain[get_loc(idx, idy, npx)] << " ";
      }
      std::cout << std::endl;
    }

    if (max_diff > 1500)
      std::cout << "Load not balanced correctly" << std::endl;
    else
      std::cout << "Load balanced correctly" << std::endl;
  }

  /**
   * function for bfs for extending the domain corresponding to each of the
   * nodes uniformly in all possible directions so that the change of shape in
   * the region corresponding to a node does not change much
   * @param node_id node id of the region whose boundary is being extended or
   * contracted
   * @param node_subdomain a square subdomain's id in the region belonging to
   * the current node
   * @param locality_ids locality id to which a particular square subdomain is
   * attached
   * @param visited_subdomain data structure to keep a track of sudomains
   * corresponding to current node
   * @param visited_node data structure to keep a track of visited nodes
   * @param work_realloc number of subdomains that need to be exchanged by
   * various nodes
   * @param total_subdomains counter for each nodes denoting number of square
   * subdomains for each of the nodes
   */
  void locality_subdomain_bfs(int node_id, int node_subdomain,
                              std::vector<std::vector<int> >& locality_ids,
                              int visited_subdomain[], int visited_node[],
                              int work_realloc[], int total_subdomains[]) {
    long ctr_2 = 1;
    long time = 0;

    // prioirity queue instead of normal queue so that we can reach boundary
    // points faster and expand all the boundaries evenly
    std::priority_queue<std::pair<bool, std::pair<long, long> > > bfs_q;
    std::pair<bool, std::pair<long, long> > q_element;
    bfs_q.push({1, {0, node_subdomain}});
    // variable to indicate whether we are currently expanding or contracting
    // the horizon
    bool invard_bfs = 0;
    time++;

    while (!bfs_q.empty() && work_realloc[node_id] != 0) {
      long ctr_1 = 0;
      for (long id = 0; id < ctr_2; ++id) {
        q_element = bfs_q.top();
        bfs_q.pop();
        long gx = q_element.ss.ss % npx, gy = q_element.ss.ss / npx;
        // we have reached the boundary and want to expand inward
        if (q_element.ff == 0 && work_realloc[node_id] < 0) invard_bfs = 1;
        visited_subdomain[q_element.ss.ss] = 1 + invard_bfs;

        for (long idx = std::max((long)0, gx - 1); idx < std::min(npx, gx + 2);
             ++idx) {
          for (long idy = std::max((long)0, gy - 1);
               idy < std::min(npy, gy + 2); ++idy) {
            if (std::abs(idx - gx) + std::abs(idy - gy) > 1 ||
                visited_subdomain[get_loc(idx, idy, npx)] == 1 + invard_bfs ||
                visited_node[locality_ids[idx][idy]] == 2)
              continue;

            bool same_locality = (locality_ids[idx][idy] == node_id);

            if (work_realloc[node_id] > 0) {
              bfs_q.push({same_locality, {time, get_loc(idx, idy, npx)}});
              ctr_1++;
            } else if (work_realloc[node_id] < 0 && invard_bfs == 0) {
              bfs_q.push({same_locality, {time, get_loc(idx, idy, npx)}});
              ctr_1++;
            } else if (work_realloc[node_id] < 0 && invard_bfs == 1 &&
                       same_locality && total_subdomains[node_id] > 1) {
              bfs_q.push({0, {time, get_loc(idx, idy, npx)}});
              ctr_1++;

              work_realloc[locality_ids[gx][gy]]--;
              work_realloc[node_id]++;

              total_subdomains[node_id]--;
              total_subdomains[locality_ids[gx][gy]]--;

              locality_ids[idx][idy] = locality_ids[gx][gy];
            }
          }
        }

        // we have reached the boundary and want to expand outward
        // at the same time, we want to expand only to those subdomains
        // whose localities have not been visited yet
        if (visited_node[locality_ids[gx][gy]] != 2 &&
            locality_ids[gx][gy] != node_id && work_realloc[node_id] > 0 &&
            total_subdomains[locality_ids[gx][gy]] > 1) {
          work_realloc[locality_ids[gx][gy]]++;
          work_realloc[node_id]--;

          total_subdomains[node_id]++;
          total_subdomains[locality_ids[gx][gy]]--;

          locality_ids[gx][gy] = node_id;
        }

        time++;
      }
      ctr_2 = ctr_1;
    }
  }

  /**
   * function for inserting the rectangles which are just above, below and on
   * the sides of the given partition into the unordered_map used to index the
   * squares already inserted into the vector
   * @param parent_id node id of the parent in the DFS
   * @param node_id node id of the region whose boundary is being extended or
   * contracted
   * @param locality_ids locality id to which a particular square subdomain is
   * attached
   * @param adj_set adjacency list for each of the nodes
   * @param visited_subdomain data structure to keep a track of sudomains
   * corresponding to current node
   * @param visited_node data structure to keep a track of visited nodes
   * @param work_realloc number of subdomains that need to be exchanged by
   * various nodes
   * @param node_subdomain std::vector containing a square subdomain's id in the
   * region belonging to the various nodes
   * @param total_subdomains counter for each nodes denoting number of square
   * subdomains for each of the nodes
   */
  void redistribution_dfs(int parent_id, int node_id,
                          std::vector<std::vector<int> >& locality_ids,
                          std::set<int> adj_set[], int visited_subdomain[],
                          int visited_node[], int work_realloc[],
                          int node_subdomain[], int total_subdomains[]) {
    visited_node[node_id] = 1;
    for (auto itr = adj_set[node_id].begin(); itr != adj_set[node_id].end();
         ++itr) {
      if (visited_node[*itr] == 0) {
        redistribution_dfs(node_id, *itr, locality_ids, adj_set,
                           visited_subdomain, visited_node, work_realloc,
                           node_subdomain, total_subdomains);
      }
    }

    // The first partition doesn't have any partition left to exchange with it
    if (parent_id == -1) return;
    // perform subdomain exchanges as required with other partitions
    locality_subdomain_bfs(node_id, node_subdomain[node_id], locality_ids,
                           visited_subdomain, visited_node, work_realloc,
                           total_subdomains);

    visited_node[node_id] = 2;
  }

  /**
   * Function to run the load balancing algorithm every nbalance time steps.
   * Function reassigns localities to various square subdomains(if any) to
   * have balanced computation across all the nodes in the distributed setup.
   * @param unbalanced_data std::vector of client-side data structures for
   * storing temperature data
   * @return std::vector of client-side data structures with some square
   * subdomains possibly reassigned to different nodes
   */
  std::vector<partition_space_client> load_balance(
      std::vector<partition_space_client> unbalanced_data) {
    // hpx::naming::get_locality_id_from_id
    double busy_rates[nl];  // idle rate values returned by performance counters
    int work_realloc[nl];   // amount of subdomains a particular node is capable
                            // of computing more
    int total_subdomains[nl];  // total number of subdomains present on the node
    double expected_busy_rate =
        0.0;  // hypotheticallly expected idle rate for equal distribution of
              // work among all nodes

    // retrieve the performance counter values for idle-rates
    for (int idx = 0; idx < nl; ++idx) {
      hpx::performance_counters::counter_value value1 =
          perf_counters[idx].get_counter_value(hpx::launch::sync, true);
      busy_rates[idx] = 10000 - value1.get_value<double>();

      expected_busy_rate += busy_rates[idx];
      total_subdomains[idx] = 0;
    }
    expected_busy_rate /= nl;

    std::set<int> adj_set[nl];  // adjacency list of the domain partitions graph
                                // as a set to prevent repetitions
    int visited_node[nl];  // true if a particular node has already been visited
                           // during DFS
    int node_subdomain[nl];  // contains subdomain id for each of the nodes to
                             // start the BFS on each node
    int visited_subdomain[npx *
                          npy];  // true if a particular subdomain has already
                                 // been visited during load balancing
    std::vector<std::vector<int> >
        locality_ids;  // locality id to which a particular element is attached
    locality_ids.resize(npx);
    for (long idx = 0; idx < npx; ++idx) locality_ids[idx].resize(npy);

    for (long gx = 0; gx < npx; ++gx) {
      for (long gy = 0; gy < npy; ++gy) {
        std::size_t locality_id = locidx(get_loc(gx, gy, npx), nl);
        bool is_boundary = 0;
        for (long idx = std::max((long)0, gx - 1); idx < std::min(npx, gx + 2);
             ++idx) {
          for (long idy = std::max((long)0, gy - 1);
               idy < std::min(npy, gy + 2); ++idy) {
            if (std::abs(idx - gx) + std::abs(idy - gy) > 1) continue;
            if (locidx(get_loc(idx, idy, npx), nl) != locality_id) {
              adj_set[locality_id].insert(locidx(get_loc(idx, idy, npx), nl));
            }
          }
        }

        node_subdomain[locality_id] = get_loc(gx, gy, npx);
        total_subdomains[locality_id]++;
        visited_subdomain[get_loc(gx, gy, npx)] = 0;
        locality_ids[gx][gy] = locality_id;
      }
    }

    int min_load_imbalance = INT32_MAX;
    int node_min_imbalance = 0;

    // calculate the excess and the deficiency of compute across all the nodes
    for (int idx = 0; idx < nl; ++idx) {
      double time_per_subdomain =
          busy_rates[idx] / (double)total_subdomains[idx];
      if (expected_busy_rate > busy_rates[idx] &&
          abs(expected_busy_rate - busy_rates[idx]) > 0.3 * time_per_subdomain)
        work_realloc[idx] =
            ceil((expected_busy_rate - busy_rates[idx]) / time_per_subdomain);
      else if (expected_busy_rate < busy_rates[idx] &&
               abs(expected_busy_rate - busy_rates[idx]) >
                   0.3 * time_per_subdomain)
        work_realloc[idx] =
            floor((expected_busy_rate - busy_rates[idx]) / time_per_subdomain);
      else
        work_realloc[idx] = 0;

      if (abs(work_realloc[idx]) < min_load_imbalance) {
        min_load_imbalance = abs(work_realloc[idx]);
        node_min_imbalance = idx;
      }

      visited_node[idx] = 0;
    }

    // positive work_realloc means requires more work for optimal distribution
    // perform the redistribution of domains and ressign the locality ids to
    // balance the load set the root node as the "compute node" with minimum
    // load imbalance
    redistribution_dfs(-1, node_min_imbalance, locality_ids, adj_set,
                       visited_subdomain, visited_node, work_realloc,
                       node_subdomain, total_subdomains);

    // perform the reassignment of localities for each of the subdomains and
    // return them
    for (long idx = 0; idx < npx; ++idx)
      for (long idy = 0; idy < npy; ++idy)
        unbalanced_data[get_loc(idx, idy, npx)] = partition_space_client(
            this->localities[locality_ids[idx][idy]],
            hpx::util::unwrap(
                unbalanced_data[get_loc(idx, idy, npx)].get_data()));

    if (partition_domain.size() == 0) partition_domain.resize(npx * npy);

    for (long idx = 0; idx < npx; ++idx)
      for (long idy = 0; idy < npy; ++idy)
        partition_domain[get_loc(idx, idy, npx)] = locality_ids[idx][idy];

    // reset all counters for more accurate time and remove the data transfer
    // times from the idle-rate
    for (int idx = 0; idx < nl; ++idx)
      hpx::performance_counters::counter_value value1 =
          perf_counters[idx].get_counter_value(hpx::launch::sync, true);

    return unbalanced_data;
  }

  /**
   * our influence function to weigh various points differently
   * @param distance distance of the point from point of computation
   * @return returns the influence of the current point for point under
   * computation
   */
  static inline double influence_function(double distance) { return 1.0; }

  /**
   * function for inserting the rectangles which are just above, below and on
   * the sides of the given partition into the unordered_map used to index the
   * squares already inserted into the vector
   * @param lx rectangle's lower left point's x coordinate
   * @param ly rectangle's lower left point's y coordinate
   * @param rx rectangle's upper right point's x coordinate
   * @param ry rectangle's upper right point's y coordinate
   * @param all_squares std::vector containing npx * npy squares from the
   * partitioned grid
   * @param neighbour_squares std::vector containing the neighbour square
   * subdomains from the current square
   */
  static void add_neighbour_rectangle(long lx, long ly, long rx, long ry,
                                      space_2d& all_squares,
                                      space_2d& neighbour_squares) {
    for (long sx = lx; sx < rx + nx; sx += nx) {
      for (long sy = ly; sy < ry + ny; sy += ny) {
        if (sx < 0 || sx >= nx * npx || sy < 0 || sy >= ny * npy) continue;

        neighbour_squares.push_back(all_squares[grid_loc(sx, sy)]);
      }
    }
  }

  /**
   * testing operator to verify correctness
   * @param pos_x x coordinate of the point of computation
   * @param pos_y y coordinate of the point of computation
   * @param time time at which we want the analytical solution
   * @return analytical solution's expected value at point (pos_x, pos_y)
   */
  static inline double w(long pos_x, long pos_y, long time) {
    return cos(2 * M_PI * (time * dt)) * sin(2 * M_PI * (pos_x * dh)) *
           sin(2 * M_PI * (pos_y * dh));
  }

  /**
   * condition to enforce the boundary conditions for the divided partitions
   * @param pos_x x coordinate of the point of computation
   * @param pos_y y coordinate of the point of computation
   * @param all_squares std::vector containing npx * npy squares from the
   * partitioned grid
   * @return 0 for points outside the mesh, else the value computed by the
   * numerical solution for point (pos_x, pos_y)
   */
  static inline double boundary(
      long pos_x, long pos_y, const std::vector<partition_space>& all_squares) {
    if (pos_x >= 0 && pos_x < nx * npx && pos_y >= 0 && pos_y < ny * npy)
      return all_squares[grid_loc(pos_x, pos_y)][mesh_loc(pos_x, pos_y)];
    else
      return 0;
  }

  /**
   * condition to enforce the boundary conditions for the tests of analytical
   * solution
   * @param pos_x x coordinate of the point of computation
   * @param pos_y y coordinate of the point of computation
   * @param val expected value for analytical solution at the point (pos_x,
   * pos_y)
   * @return 0 for points outside the mesh, else the value computed by the
   * analytical solution for point (pos_x, pos_y)
   */
  static inline double boundary(long pos_x, long pos_y, double val) {
    if (pos_x >= 0 && pos_x < nx * npx && pos_y >= 0 && pos_y < ny * npy)
      return val;
    else
      return 0;
  }

  /**
   * Function to find distance in 2d plane given 2 points coordinates
   * @param cen_x x coordinate of the point of computation
   * @param cen_y y coordinate of the point of computation
   * @param pos_x x coordinate of the point that influences the temperature at
   * (cen_x, cen_y)
   * @param pos_y y coordinate of the point that influences the temperature at
   * (cen_x, cen_y)
   * @return distance of point (cen_x, cen_y) from point (pos_x, pos_y)
   */
  static inline double distance(long cen_x, long cen_y, long pos_x,
                                long pos_y) {
    return sqrt(((cen_x - pos_x) * (cen_x - pos_x)) +
                ((cen_y - pos_y) * (cen_y - pos_y)));
  }

  //! Function to compute length of 3rd side of a right angled triangle
  //! given hypotenuse and lenght of one side
  static inline double len_1d_line(long len_x) {
    return sqrt((eps * eps) - (len_x * len_x));
  }

  /**
   * function adds the external source to the 2d nonlocal heat equation for
   * testing correctness of numerical solution against the analytical solution
   * @param pos_x x coordinate of the point whose numerical solution is to be
   * computed
   * @param pos_y y coordinate of the point whose numerical solution is to be
   * computed
   * @param time time at which we want the numerical solution
   * @return temperature value at point (pos_x, pos_y) for analytical solution
   */
  static double sum_local_test(long pos_x, long pos_y, long time) {
    double result_local =
        -(2 * M_PI * sin(2 * M_PI * (time * dt)) *
          sin(2 * M_PI * (pos_x * dh)) * sin(2 * M_PI * (pos_y * dh)));
    double w_position = w(pos_x, pos_y, time);

    for (long sx = pos_x - eps; sx <= pos_x + eps; ++sx) {
      long len_line = len_1d_line(std::abs((long)pos_x - (long)sx));
      for (long sy = pos_y - len_line; sy <= pos_y + len_line; ++sy) {
        result_local -=
            influence_function(distance(pos_x, pos_y, sx, sy)) * c_2d *
            (boundary(sx, sy, w(sx, sy, time)) - w_position) * (dh * dh);
      }
    }

    return result_local;
  }

  /**
   * Our operator to find sum of 'eps' radius circle in vicinity of point P(x,y)
   * Represent circle as a series of horizaontal lines of thickness 'dh'
   * @param pos_x x coordinate of the point whose numerical solution is to be
   * computed
   * @param pos_y y coordinate of the point whose numerical solution is to be
   * computed
   * @param all_squares std::vector containing npx * npy squares from the
   * partitioned grid
   * @return temperature value at point (pos_x, pos_y) for numerical solution
   */
  static double sum_local(long pos_x, long pos_y,
                          const std::vector<partition_space>& all_squares) {
    double result_local = 0.0;
    double pos_val = boundary(pos_x, pos_y, all_squares);

    for (long sx = pos_x - eps; sx <= pos_x + eps; ++sx) {
      long len_line = len_1d_line(std::abs((long)pos_x - (long)sx));
      for (long sy = pos_y - len_line; sy <= pos_y + len_line; ++sy) {
        result_local += influence_function(distance(pos_x, pos_y, sx, sy)) *
                        c_2d * (boundary(sx, sy, all_squares) - pos_val) *
                        (dh * dh);
      }
    }

    return result_local;
  }

  //! Function to retrieve the data from possibly remote partitions for
  //! vector of square sundomains
  static space_2d_fut vector_get_data(
      const std::vector<partition_space_client>& neighbour_squares) {
    space_2d_fut ret_get_data;
    ret_get_data.resize(neighbour_squares.size());

    for (long idx = 0; idx < neighbour_squares.size(); ++idx) {
      ret_get_data[idx] = neighbour_squares[idx].get_data();
    }

    return ret_get_data;
  }

  /**
   * partition of nx * ny cells which will be processed by a single thread to
   * enable full utilization of all cores across all nodes in the cluster
   * @param middle x coordinate of the point whose numerical solution is to be
   * computed
   * @param neighbour_squares std::vector containing the neighbour square
   * subdomains from the current square
   * @param time time for which we are logging the results
   * @param test variable is true if we have to test the numerical solution
   * against the analytical solution
   * @return client side data structure for temperature values for a square
   * subdomain
   */
  static partition_space_client sum_local_partition(
      partition_space_client const& middle,
      const std::vector<partition_space_client>& neighbour_squares, long time,
      bool test) {
    // all_squares is a vector with partitions with pointers to valid partitions
    // to only those neighbours who affect the temperature at the middle
    // partition

    hpx::shared_future<partition_space> middle_data = middle.get_data();

    hpx::future<partition_space> next_middle =
        middle_data.then(hpx::util::unwrapping(
            [time, test](partition_space const& mid) -> partition_space {
              std::vector<partition_space> all_squares;
              all_squares.resize(npx * npy);
              // All local operations are performed once the middle data of
              // the previous time step becomes available.
              long gx = mid.gx, gy = mid.gy;
              partition_space next(mid.size(), gx, gy);
              all_squares[get_loc(gx, gy, npx)] = mid;

              for (long sx = eps + 1; sx < nx - eps - 1; ++sx) {
                for (long sy = eps + 1; sy < ny - eps - 1; ++sy) {
                  next[get_loc(sx, sy)] =
                      mid[get_loc(sx, sy)] +
                      (sum_local(sx + gx * nx, sy + gy * ny, all_squares) * dt);
                  if (test)
                    next[get_loc(sx, sy)] +=
                        sum_local_test(sx + gx * nx, sy + gy * ny, time) * dt;
                }
              }
              return next;
            }));

    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapping([middle, time, test](
                                  partition_space next,
                                  partition_space const& mid,
                                  const std::vector<partition_space>&
                                      neighbour_squares)
                                  -> partition_space_client {
          // Calculate the missing boundary elements once the
          // corresponding data has become available.
          std::vector<partition_space> all_squares;
          all_squares.resize(npx * npy);
          long gx = mid.gx, gy = mid.gy;

          // assign only those partition pointers which are valid neighbours for
          // a particular square in the grid
          for (long idx = 0; idx < neighbour_squares.size(); ++idx) {
            all_squares[get_loc(neighbour_squares[idx].gx,
                                neighbour_squares[idx].gy, npx)] =
                neighbour_squares[idx];
          }

          if (nx <= eps) {
            for (long sx = 0; sx < nx; ++sx) {
              for (long sy = 0; sy < ny; ++sy) {
                next[get_loc(sx, sy)] =
                    mid[get_loc(sx, sy)] +
                    (sum_local(sx + gx * nx, sy + gy * ny, all_squares) * dt);
                if (test)
                  next[get_loc(sx, sy)] +=
                      sum_local_test(sx + gx * nx, sy + gy * ny, time) * dt;
              }
            }
          } else {
            for (long sx = 0; sx < eps + 1; ++sx) {
              for (long sy = 0; sy < ny; ++sy) {
                next[get_loc(sx, sy)] =
                    mid[get_loc(sx, sy)] +
                    (sum_local(sx + gx * nx, sy + gy * ny, all_squares) * dt);
                if (test)
                  next[get_loc(sx, sy)] +=
                      sum_local_test(sx + gx * nx, sy + gy * ny, time) * dt;
              }
            }

            for (long sx = eps + 1; sx < nx - eps - 1; ++sx) {
              for (long sy = 0; sy < eps + 1; ++sy) {
                next[get_loc(sx, sy)] =
                    mid[get_loc(sx, sy)] +
                    (sum_local(sx + gx * nx, sy + gy * ny, all_squares) * dt);
                if (test)
                  next[get_loc(sx, sy)] +=
                      sum_local_test(sx + gx * nx, sy + gy * ny, time) * dt;
              }
              for (long sy = ny - eps - 1; sy < ny; ++sy) {
                next[get_loc(sx, sy)] =
                    mid[get_loc(sx, sy)] +
                    (sum_local(sx + gx * nx, sy + gy * ny, all_squares) * dt);
                if (test)
                  next[get_loc(sx, sy)] +=
                      sum_local_test(sx + gx * nx, sy + gy * ny, time) * dt;
              }
            }

            for (long sx = nx - eps - 1; sx < nx; ++sx) {
              for (long sy = 0; sy < ny; ++sy) {
                next[get_loc(sx, sy)] =
                    mid[get_loc(sx, sy)] +
                    (sum_local(sx + gx * nx, sy + gy * ny, all_squares) * dt);
                if (test)
                  next[get_loc(sx, sy)] +=
                      sum_local_test(sx + gx * nx, sy + gy * ny, time) * dt;
              }
            }
          }

          // The new partition_data will be allocated on the same locality
          // as 'middle'.
          return partition_space_client(middle.get_id(), next);
        }),
        std::move(next_middle), middle_data,
        vector_get_data(neighbour_squares));
  }

  void do_work();
};

HPX_PLAIN_ACTION(solver::sum_local_partition, sum_local_partition_action);

//! Compute the numerical solution on 'nx * ny' data points for each of 'npx *
//! npy' partitions for 'nt' time steps
void solver::do_work() {
  // Actual time step loop
  for (long t = 0; t < nt; ++t) {
    current = t % 2;
    next = (t + 1) % 2;

    sum_local_partition_action act;
    for (long gx = 0; gx < npx; ++gx) {
      for (long gy = 0; gy < npy; ++gy) {
        // we execute the action on the locality of the middle partition
        auto Op = hpx::util::bind_front(
            act, localities[locidx(get_loc(gx, gy, npx), nl)]);

        // vector of dependent grid squares for a particular grid squares
        space_2d vector_deps;

        // add dependent neighbouring squares
        // for now they are added assuming it to be a big rectangle
        // in the future these additions can be in form of 4 sectors for 4
        // corners and 4 rectangles for top, bottom, left and right of the
        // rectangle
        add_neighbour_rectangle(gx * nx - eps, gy * ny - eps,
                                gx * nx + nx + eps, gy * ny + ny + eps,
                                S[current], vector_deps);

        // represent dependencies using hpx dataflow
        S[next][get_loc(gx, gy, npx)] = hpx::dataflow(
            hpx::launch::async, Op, S[current][get_loc(gx, gy, npx)],
            vector_deps, t, test);
      }
    }

    // launch a seperate thread for running the load balancing algorithm
    // and assign back the changed array of localities(if any) by the
    // load balancing algorithm
    if (t % nbalance == 0 && t != 0 && nl > 1) {
      for (int idx = 0; idx < npx * npy; ++idx) S[next][idx].get_data().wait();
      S[next] = load_balance(S[next]);
    }

    // launch a seperate thread for collecting the logging data from
    // various threads and output in the required files
    if (t % nlog == 0)
      hpx::dataflow(hpx::launch::async,
                    hpx::util::unwrapping(&solver::log_csv_vtk),
                    vector_get_data(S[next]), P, t, test);
  }

  // required for correctness for the case when nt = 0
  // because in that case next won't be set in the above for loop
  next = nt % 2;

  // wait for the solution to be ready for all the partitions
  for (int idx = 0; idx < npx * npy; ++idx) S[next][idx].get_data().wait();
}

//! function to execute the ctests which are there in file '2d_distributed.txt'
int batch_tester(long nlog, long nbalance) {
  std::uint64_t nt, num_tests;
  bool test_failed = 0;

  std::cin >> num_tests;
  for (long i = 0; i < num_tests; ++i) {
    std::cin >> nx >> ny >> npx >> npy >> nt >> eps >> k >> dt >> dh;

    // Create the solver object
    solver solve(nt, nlog, nbalance, 1);

    // do the actual error
    solve.do_work();

    // compute and print the error for checking the difference between
    // analytical and numerical solution
    solve.compute_l2(nt);

    if (solve.error_l2 / (double)(nx * ny * npx * npy) > 1e-6) {
      test_failed = 1;
      break;
    }
  }

  // output regular expression whether the test passed or failed
  if (test_failed)
    std::cout << "Tests Failed" << std::endl;
  else
    std::cout << "Tests Passed" << std::endl;

  return hpx::finalize();
}

//! main function in the HPX runtime
int hpx_main(hpx::program_options::variables_map& vm) {
  setup_counters();

  std::uint64_t nt = vm["nt"].as<std::uint64_t>();  // Number of steps.
  std::uint64_t nlog =
      vm["nlog"]
          .as<std::uint64_t>();  // Number of time steps to log the results
  std::uint64_t nbalance =
      vm["nbalance"].as<std::uint64_t>();  // Number of time steps to run the
                                           // load balancing algorithm
  std::string filename =
      vm["file"].as<std::string>();   // Filename to take input from
  bool test = vm["test"].as<bool>();  // Boolean variable to indicate if
                                      // numerical solution needs to be tested
  if (nx <= eps)
    std::cout << "[WARNING] Mesh size on a single node (nx * ny) is too small "
                 "for given epsilon (eps)"
              << std::endl;

  if (vm.count("no-header")) header = false;

  // batch testing for ctesting
  if (vm.count("test_batch")) return batch_tester(nlog, nbalance);

  // Create the solver object
  solver solve(nt, nlog, nbalance, test, filename);

  // Measure execution time.
  std::uint64_t t = hpx::util::high_resolution_clock::now();

  // Execute nt time steps on nx grid points.
  solve.do_work();

  std::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;

  if (vm.count("test_load_balance")) solve.test_load_balance();

  // print the calculated error
  // print comparison only when the 'cmp' flag is set
  if (test) solve.print_error(vm["cmp"].as<bool>());

  // Print the final solution
  if (vm.count("results")) solve.print_soln();

  std::uint64_t const num_worker_threads = hpx::get_num_worker_threads();
  hpx::future<std::uint32_t> locs = hpx::get_num_localities();
  print_time_results(locs.get(), num_worker_threads, elapsed, nx, ny, npx, npy,
                     nt, header);

  return hpx::finalize();
}

//! main function to create and initialize an HPX runtime
int main(int argc, char* argv[]) {
  std::cout << argv[0] << " (" << MAJOR_VERSION << "." << MINOR_VERSION << "."
            << UPDATE_VERSION << ")" << std::endl;
  namespace po = hpx::program_options;

  po::options_description desc_commandline;
  desc_commandline.add_options()(
      "test", po::value<bool>()->default_value(true),
      "use arguments from numerical solution for testing (default: false)")(
      "test_batch",
      "test the solution against numerical solution against batch inputs "
      "(default: false)")("test_load_balance",
                          "test the load balancing algorithm for approximately "
                          "equal load balancing "
                          "(default: false)")(
      "results", "print generated results (default: false)")(
      "cmp", po::value<bool>()->default_value(false),
      "Print expected versus actual outputs")(
      "file", po::value<std::string>()->default_value("None"),
      "name of file to take input from")(
      "nx", po::value<long>(&nx)->default_value(25), "Local x dimension")(
      "ny", po::value<long>(&ny)->default_value(25), "Local y dimension")(
      "nt", po::value<std::uint64_t>()->default_value(45),
      "Number of time steps")("npx", po::value<long>(&npx)->default_value(2),
                              "Number of partitions in x dimension")(
      "npy", po::value<long>(&npy)->default_value(2),
      "Number of partitions in y dimension")(
      "nlog", po::value<std::uint64_t>()->default_value(5),
      "Number of time steps to log the results")(
      "nbalance", po::value<std::uint64_t>()->default_value(INT64_MAX),
      "Number of time steps to run the load balancing algorithm")(
      "eps", po::value<long>(&eps)->default_value(5),
      "Epsilon for nonlocal equation")(
      "k", po::value<double>(&k)->default_value(1),
      "Heat transfer coefficient (default: 0.5)")(
      "dt", po::value<double>(&dt)->default_value(0.0005),
      "Timestep unit (default: 1.0[s])")(
      "dh", po::value<double>(&dh)->default_value(0.05),
      "Quantization across space")("no-header",
                                   "do not print out the csv header row");

  // Initialize and run HPX
  return hpx::init(desc_commandline, argc, argv);
}