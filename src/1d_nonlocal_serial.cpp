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
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <fstream>

#include "../include/point.h"
#include "../include/print_time_results.hpp"
#include "../include/writer.h"
#include "../include/Config.h"

bool header = 0;
double k = 1;    // heat transfer coefficient
double dt = 1.;  // time step
double dx = 1.;  // grid spacing

class solver {
 public:
  // 3d coordinate for representation of point
  typedef util::Point3 point_3d;

  // 3d line using the 3d point representation
  typedef std::vector<point_3d> line_3d;

  // Our temperature data type
  typedef double temperature;

  // Our data for one time step
  typedef std::vector<temperature> space_1d;

  // alternate for space between time steps,
  // result of S[t%2] goes to S[(t+1)%2] at every timestep 't'
  space_1d S[2];

  // vector containining 3d coordinates of the points on the line
  line_3d P;

  // nx = number of data points
  // nt = number of timesteps
  // eps = Epsilon for influence zone of a point
  // nlog = number of timesteps after which we want to log the states
  long nx, nt, eps, c_1d, nlog;
  bool current, next, test;
  // l2 norm and l infinity norm
  double error_l2, error_linf;

  // file to store the simulation results in csv format
  const std::string simulate_fname = "../out_csv/simulate_1d.csv";

  // file to store l2 and l infinity norms per timestep
  const std::string score_fname = "../out_csv/score_1d.csv";

  // constructor to initialize class variables
  solver(long nx, long nt, long eps, long nlog) {
    this->nx = nx;
    this->nt = nt;
    this->eps = eps;
    this->nlog = nlog;
    this->c_1d = (k * 3) / (pow(eps * dx, 3));
    this->current = 0;
    this->error_l2 = 0.0;
    this->error_linf = 0.0;
    this->test = 0;
    this->next = 1;

    P.resize(nx);
    for (long sx = 0; sx < nx; ++sx) {
      P[sx] = point_3d(sx, 0, 0);
    }

    for (auto& s : S) {
      s.resize(nx);
    }
  }

  void compute_l2(long time) {
    error_l2 = 0;

    for (long sx = 0; sx < nx; ++sx)
      error_l2 += (S[next][sx] - w(sx, time)) * (S[next][sx] - w(sx, time));
  }

  void compute_linf(long time) {
    error_linf = 0;

    for (long sx = 0; sx < nx; ++sx)
      error_linf = std::max(std::abs(S[next][sx] - w(sx, time)), error_linf);
  }

  // print error for testing
  void print_error(bool cmp) {
    std::cout << "l2: " << error_l2 << " linfinity: " << error_linf
              << std::endl;
    if (cmp)
      for (long sx = 0; sx < nx; ++sx)
        std::cout << "Expected: " << w(sx, nt) << " Actual: " << S[nt % 2][sx]
                  << std::endl;
  }

  // input the initialization
  void input_init() {
    this->test = 0;
    for (auto& sx : S[0]) {
      std::cin >> sx;
    }
  }

  // init for testing the 1d nonlocal equation
  void test_init() {
    this->test = 1;
    for (long sx = 0; sx < nx; ++sx) {
      S[0][sx] = sin(2 * M_PI * (sx * dx));
    }
  }

  // Function to visualize in csv format and conduct various experiments
  void log_vtk(long log_num) {
    const std::string fname = "../out_vtk/simulate_" + std::to_string(log_num);
    rw::writer::VtkWriter vtk_logger(fname);

    vtk_logger.appendNodes(&P);
    vtk_logger.appendPointData("Temperature", &S[next]);
    vtk_logger.addTimeStep(std::time(0));
    vtk_logger.close();
  }

  // Function to visualize in csv format and conduct various experiments
  void log_csv(long time) {
    std::ofstream outfile;
    outfile.open(simulate_fname, std::ios_base::app);

    for (long sx = 0; sx < nx; ++sx) {
      outfile << time << "," << sx << "," << S[next][sx] << "," << w(sx, time)
              << ","
              << (S[next][sx] - w(sx, time)) * (S[next][sx] - w(sx, time))
              << "," << std::abs(S[next][sx] - w(sx, time)) << ",\n";
    }

    outfile.close();

    // add to the score csv only when it's required
    if (test) {
      compute_l2(time);
      compute_linf(time);

      std::ofstream outfile;

      outfile.open(score_fname, std::ios_base::app);
      outfile << time << "," << error_l2 << "," << error_linf << ",\n";
      outfile.close();
    }
  }

  // our influence function to weigh various points differently
  static double influence_function(double distance) { return 1; }

  // testing operator to verify correctness
  inline double w(long position, long time) {
    return cos(2 * M_PI * (time * dt)) * sin(2 * M_PI * (position * dx));
  }

  // condition to enforce the boundary conditions
  inline double boundary(long pos_x, double val) {
    if (pos_x >= 0 && pos_x < nx)
      return val;
    else
      return 0;
  }

  // our test by introducing an external source
  double sum_local_test(long position, long time) {
    double result_local = -(2 * M_PI * sin(2 * M_PI * (time * dt)) *
                            sin(2 * M_PI * (position * dx)));
    double w_position = w(position, time);
    for (long sx = position - eps; sx <= position + eps; ++sx) {
      result_local -= influence_function(std::abs((long)position - (long)sx)) *
                      c_1d * (boundary(sx, w(sx, time)) - w_position) * (dx);
    }
    return result_local;
  }

  // Our operator for getting result from local points on the locality
  double sum_local(long position) {
    double result_local = 0.0;
    for (long sx = position - eps; sx <= position + eps; ++sx) {
      result_local +=
          influence_function(std::abs((long)position - (long)sx)) * c_1d *
          (boundary(sx, S[current][sx]) - S[current][position]) * (dx);
    }
    return result_local;
  }

  // do all the work on 'nx' data points for 'nt' time steps
  space_1d do_work() {
    // Actual time step loop
    for (long t = 0; t < nt; ++t) {
      current = t % 2;
      next = (t + 1) % 2;

      for (long sx = 0; sx < nx; ++sx) {
        S[next][sx] = S[current][sx] + (sum_local(sx) * dt);
        if (test) S[next][sx] += sum_local_test(sx, t) * dt;
      }

      if (t % nlog == 0) {
        log_vtk(t / nlog);
        log_csv(t);
      }
    }

    next = nt % 2;

    // testing the code for correctness
    if (test) {
      compute_l2(nt);
      compute_linf(nt);
    }

    // Return the solution at time-step 'nt'.
    return S[next];
  }
};

int batch_tester(long nlog) {
  std::uint64_t nx, nt, eps, num_tests;
  bool test_failed = 0;

  std::cin >> num_tests;
  for (long sx = 0; sx < num_tests; ++sx) {
    std::cin >> nx >> nt >> eps >> k >> dt >> dx;

    // Create the solver object
    solver solve(nx, nt, eps, nlog);
    solve.test_init();

    solve.do_work();

    if (solve.error_l2 / (double)nx > 1e-6) {
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

int hpx_main(hpx::program_options::variables_map& vm) {
  std::uint64_t nx = vm["nx"].as<std::uint64_t>();  // Number of grid points.
  std::uint64_t nt = vm["nt"].as<std::uint64_t>();  // Number of steps.
  std::uint64_t eps =
      vm["eps"].as<std::uint64_t>();  // Epsilon for influence zone of a point
  std::uint64_t nlog =
      vm["nlog"]
          .as<std::uint64_t>();  // Number of time steps to log the results

  if (vm.count("no-header")) header = false;

  // batch testing for ctesting
  if (vm.count("test_batch")) return batch_tester(nlog);

  // Create the solver object
  solver solve(nx, nt, eps, nlog);

  // Take inputs from stdin for testing
  if (vm.count("test"))
    solve.test_init();
  else
    solve.input_init();

  // Measure execution time.
  std::uint64_t t = hpx::util::high_resolution_clock::now();

  // Execute nt time steps on nx grid points.
  solver::space_1d solution = solve.do_work();

  std::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;

  // print the calculated error for testing
  if (vm.count("test")) solve.print_error(vm["cmp"].as<bool>());

  // Print the final solution
  if (vm.count("results"))
    for (long sx = 0; sx < nx; ++sx)
      std::cout << "S[" << sx << "] = " << solution[sx] << std::endl;

  std::uint64_t const os_thread_count = hpx::get_os_thread_count();
  print_time_results(os_thread_count, elapsed, nx, nt, header);

  return hpx::finalize();
}

int main(int argc, char* argv[]) {
  std::cout << argv[0] << " (" << MAJOR_VERSION << "." << MINOR_VERSION << "."
            << UPDATE_VERSION << ")" << std::endl;
  namespace po = hpx::program_options;

  po::options_description desc_commandline;
  desc_commandline.add_options()(
      "test",
      "use arguments from numerical solution for testing (default: false)")(
      "test_batch",
      "test the solution against numerical solution against batch inputs "
      "(default: false)")("results",
                          "print generated results (default: false)")(
      "cmp", po::value<bool>()->default_value(true),
      "Print expected versus actual outputs")(
      "nx", po::value<std::uint64_t>()->default_value(50), "Local x dimension")(
      "nt", po::value<std::uint64_t>()->default_value(45),
      "Number of time steps")("nlog",
                              po::value<std::uint64_t>()->default_value(5),
                              "Number of time steps to log the results")(
      "eps", po::value<std::uint64_t>()->default_value(5),
      "Epsilon for nonlocal equation")(
      "k", po::value<double>(&k)->default_value(1),
      "Heat transfer coefficient (default: 0.5)")(
      "dt", po::value<double>(&dt)->default_value(0.001),
      "Timestep unit (default: 1.0[s])")(
      "dx", po::value<double>(&dx)->default_value(0.02), "Local x dimension")(
      "no-header", "do not print out the csv header row");

  // Initialize and run HPX
  return hpx::init(desc_commandline, argc, argv);
}
