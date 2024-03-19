#include <fstream>
#include <iostream>

#include "aux.hpp"
#include "project.hpp"

int main(int argc, char* argv[]) {
  // Define variables to hold parameters
  std::string Gtau_file;
  std::string Gtau_psd_file;
  int max_it{};
  int flavors{};
  int ntau{};
  double tol{};
  bool verbose{false};

  // Parse parameters
  define_parameters(argc, argv, Gtau_file, Gtau_psd_file, max_it, flavors, ntau,
                    tol, verbose);

  // Print parameters
  std::cout << "##### Parameters #####" << std::endl;
  std::cout << "Gtau file: " << Gtau_file << std::endl;
  std::cout << "Gtau_psd file: " << Gtau_psd_file << std::endl;
  std::cout << "max_it: " << max_it << std::endl;
  std::cout << "flavors: " << flavors << std::endl;
  std::cout << "ntau: " << ntau << std::endl;
  std::cout << "tol: " << tol << std::endl;
  std::cout << "verbose: " << verbose << std::endl;

  // Read Gtau from file
  std::vector<std::vector<double>> Gtau(flavors, std::vector<double>(ntau, 0.));
  double beta = read_tau_from_file(Gtau, Gtau_file, verbose);

  // Do the projection
  std::cout << "##### Projecting #####" << std::endl;
  for (int f = 0; f < flavors; ++f) {
    double density = -Gtau[f][ntau - 1];
    std::cout << "---------------------" << std::endl;
    std::cout << "flavor: " << f << std::endl;
    std::cout << "density forced to: " << density << std::endl;
    hankel_dyproject(Gtau[f], max_it, density, tol, verbose);
  }
  dump_tau_to_file(Gtau, beta, Gtau_psd_file);
}
