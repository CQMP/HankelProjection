#include "aux.hpp"

void define_parameters(int argc, char *argv[], std::string &Gtau_file,
                       std::string &Gtau_psd_file, int &max_it, int &flavors,
                       int &ntau, double &tol, bool &verbose) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string arg = argv[i];
    try {
      if (arg == "--Gtau") {
        Gtau_file = argv[++i];
      } else if (arg == "--Gtau_psd") {
        Gtau_psd_file = argv[++i];
      } else if (arg == "--max_it") {
        max_it = std::stoi(argv[++i]);
      } else if (arg == "--flavors") {
        flavors = std::stoi(argv[++i]);
      } else if (arg == "--ntau") {
        ntau = std::stoi(argv[++i]);
      } else if (arg == "--tol") {
        tol = std::stod(argv[++i]);
      } else if (arg == "--verbose") {
        verbose = std::stoi(argv[++i]);
      } else {
        std::cerr << "Unknown option: " << arg << std::endl;
        exit(EXIT_FAILURE);
      }
    } catch (const std::invalid_argument &e) {
      std::cerr << "Invalid argument for option " << arg << ": " << e.what()
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  if (Gtau_file.empty() || Gtau_psd_file.empty() || max_it == 0 ||
      flavors == 0 || ntau == 0 || tol == 0) {
    std::cerr << "Missing required option(s)" << std::endl;
    exit(EXIT_FAILURE);
  }
}

double read_tau_from_file(std::vector<std::vector<double>> &G_tau,
                          const std::string &filename, bool verbose) {
  int flavors = G_tau.size();
  int n_tau = G_tau[0].size();

  std::ifstream in_file(filename);
  double tau{};
  std::cout << "##### Parsing input file #####" << std::endl;
  try {
    for (int i = 0; i < n_tau; ++i) {
      if (!(in_file >> tau)) {
        throw std::runtime_error(
            "Insufficient data in input file for specified flavors and n_tau.");
      }
      if (verbose) std::cout << tau << " ";
      for (int f = 0; f < flavors; ++f) {
        if (!(in_file >> G_tau[f][i])) {
          throw std::runtime_error(
              "Insufficient data in input file for specified flavors and "
              "n_tau.");
        }
        if (verbose) std::cout << G_tau[f][i] << " ";
      }
      if (verbose) std::cout << std::endl;
    }
    // Check if there are any additional elements in the file
    if (in_file >> tau) {
      throw std::runtime_error(
          "Excess data in input file for specified flavors and n_tau.");
    }
  } catch (const std::ifstream::failure &e) {
    // Exception handling for file read failures
    throw std::runtime_error("Problem reading from input file: " +
                             std::string(e.what()));
  }
  return tau;  // last tau corresponds to beta
}

void dump_tau_to_file(const std::vector<std::vector<double>> &G_tau,
                      double beta, const std::string &filename) {
  std::ofstream file(filename.c_str());
  int N_t = G_tau[0].size() - 1;
  int flavors = G_tau.size();
  for (std::size_t t = 0; t < N_t + 1; ++t) {
    double tau = t * beta / N_t;
    file << tau << " ";
    for (std::size_t f = 0; f < flavors; ++f) {
      file << G_tau[f][t] << " ";
    }
    file << std::endl;
  }
}
