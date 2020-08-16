#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include <gmsh.h>
#include <metis.h>

//! operator for getting the location of 2d point in 1d representation
/**
 * @param x x coordinate in 2d cartesian coordinates
 * @param y y coordinate in 2d cartesian coordinates
 * @param nx length of the 2d cartesian plane in x dimension
 * @return location of 2d point in 1d array used to store the 2d mesh
 */
inline int get_loc(int x, int y, int nx) { return x + y * nx; }

/**
 * Function for writing out the domain decomposition of the coarse mesh
 * for efficient computation in a distributed setup
 * as opposed to the orignal fine mesh, along with the number if square
 * subdomains and the size of each square subdomain
 * @param eparts locality index to which the square subdomain belongs to
 * @param filename name of file to write the coarse mesh
 * @param mx size of each square subdomain in x dimension
 * @param my size of each square subdomain in y dimension
 * @param npx number of square subdomains along x dimension
 * @param npy number of square subdomains along y dimension
 * @param dh space discretization parameter for both x and y dimension
 */
void write_mesh(idx_t eparts[], char* filename, long mx, long my, long npx,
                long npy, double dh) {
  std::ofstream outfile;
  outfile.open(filename);

  if (outfile) {
    outfile << mx / npx << " " << my / npy << " " << npx << " " << npy << " "
            << dh << std::endl;

    for (long idx = 0; idx < npx; ++idx) {
      for (long idy = 0; idy < npy; ++idy) {
        // makes intermediate file more readable
        outfile << idx << " " << idy << " " << eparts[get_loc(idx, idy, npx)]
                << std::endl;
      }
    }
  }

  outfile.close();
}

int main(int argc, char* argv[]) {
  // Note: For the current implementation, the mesh is always assumed to be of
  // rectangle shape
  if (argc < 4) {
    std::cout << "Usage: " << argv[0]
              << " mesh_filename out_filename number_compute_nodes"
              << std::endl;
    return 0;
  }

  char* mesh_filename = argv[1];
  char* out_filename = argv[2];
  idx_t num_compute_nodes =
      atoi(argv[3]);  // number of computational nodes available for
                      // distributing the equation

  gmsh::initialize();
  gmsh::option::setNumber("General.Terminal", 1);
  gmsh::open(mesh_filename);

  // getting all nodes using GMSH API
  std::vector<std::size_t> nodeTags;
  std::vector<double> nodeCoords, nodeParams;
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1);

  // getting all elements using GMSH API
  std::vector<int> elemTypes;
  std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
  gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, -1, -1);

  idx_t idx_square_elements;  // index for square elements stored in the results
                              // from GMSH API
  idx_t mx, my;               // mesh sizes along x and y dimensions
  double dh;                  // discretization parameter for the mesh
  idx_t npx, npy;  // number of partitions along x and y dimension respectively
  double minx, maxx, miny, maxy;  // min and max values for x and y coordinates

  for (long idx = 0; idx < elemTypes.size(); ++idx) {
    // type for square elements is 3 in GMSH library
    if (elemTypes[idx] == 3) {
      idx_square_elements = idx;
      break;
    }
  }

  // calculating the difference between consecutive x and y coordinates to get
  // the discretization length (dh)
  dh = nodeCoords[(elemNodeTags[idx_square_elements][0] - 1) * 3] -
       nodeCoords[(elemNodeTags[idx_square_elements][1] - 1) * 3];
  dh = std::max(
      dh,
      std::abs(nodeCoords[(elemNodeTags[idx_square_elements][0] - 1) * 3 + 1] -
               nodeCoords[(elemNodeTags[idx_square_elements][1] - 1) * 3 + 1]));

  maxx = maxy = __DBL_MIN__;
  minx = miny = __DBL_MAX__;
  for (long idx = 0; idx < elemNodeTags[idx_square_elements].size(); idx++) {
    minx = std::min(
        minx, nodeCoords[(elemNodeTags[idx_square_elements][idx] - 1) * 3]);
    maxx = std::max(
        maxx, nodeCoords[(elemNodeTags[idx_square_elements][idx] - 1) * 3]);

    miny = std::min(
        miny, nodeCoords[(elemNodeTags[idx_square_elements][idx] - 1) * 3 + 1]);
    maxy = std::max(
        maxy, nodeCoords[(elemNodeTags[idx_square_elements][idx] - 1) * 3 + 1]);
  }

  mx = round((maxx - minx) / dh);
  my = round((maxy - miny) / dh);

  // getting the size of coarse mesh size
  std::cout << "\nSize of mesh is as follows:\n";
  std::cout << "x dimension : " << mx << "\ny dimension : " << my;

  std::cout << "\n\nNote:" << std::endl;
  std::cout << "Enter coarse grain size as a divisor for size of mesh along "
               "repective dimension"
            << std::endl;
  std::cout
      << "1 <= Coarse grain size on x-dimension <= Mesh size on x-dimension"
      << std::endl;
  std::cout
      << "1 <= Coarse grain size on y-dimension <= Mesh size on y-dimension"
      << std::endl;

  std::cout << "\n\nEnter coarse mesh size along x-dimension" << std::endl;
  std::cin >> npx;
  if (mx % npx != 0) {
    std::cout << "Mesh size along x direction not divisible by the coarse "
                 "grain size in x-direction"
              << std::endl;
    return 0;
  }
  npx = mx / npx;

  std::cout << "\n\nEnter coarse mesh size along y-dimension" << std::endl;
  std::cin >> npy;
  if (my % npy != 0) {
    std::cout << "Mesh size along y direction not divisible by the coarse "
                 "grain size in y-direction"
              << std::endl;
    return 0;
  }
  npy = my / npy;

  idx_t num_nodes_elements =
      npx * npy * 4;  // total number of nodes for all elements in the mesh
  idx_t num_elements =
      npx * npy;              // total number of elements in the coarse mesh
  idx_t objval, ncommon = 1;  // variables required for METIS API

  // data structures required for METIS API, refer to METIS manual
  idx_t eparts[num_elements], nparts[num_nodes_elements],
      eptr[num_elements + 1], eind[num_nodes_elements];

  // METIS API gives floating point exception with num_compute_nodes = 1
  for (long idx = 0; idx < num_elements; ++idx) eparts[idx] = 0;
  if (num_compute_nodes < 2) goto writemesh;

  for (long idx = 0; idx < npx; ++idx) {
    for (long idy = 0; idy < npy; ++idy) {
      eind[get_loc(idx, idy, npx) * 4 + 0] = get_loc(idx, idy, npx + 1);
      eind[get_loc(idx, idy, npx) * 4 + 1] = get_loc(idx + 1, idy, npx + 1);
      eind[get_loc(idx, idy, npx) * 4 + 2] = get_loc(idx + 1, idy + 1, npx + 1);
      eind[get_loc(idx, idy, npx) * 4 + 3] = get_loc(idx, idy + 1, npx + 1);
    }
  }

  for (long idx = 0; idx < num_elements + 1; idx++) {
    eptr[idx] = 4 * idx;
  }

  METIS_PartMeshDual(&num_elements, &num_nodes_elements, eptr, eind, NULL, NULL,
                     &ncommon, &num_compute_nodes, NULL, NULL, &objval, eparts,
                     nparts);

writemesh:
  gmsh::clear();
  gmsh::finalize();

  write_mesh(eparts, out_filename, mx, my, npx, npy, dh);
  return 0;
}