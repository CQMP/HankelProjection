#pragma once
#include <Eigen/Dense>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

// auxiliary routines
void define_parameters(int argc, char *argv[], std::string &Gtau_file,
                       std::string &Gtau_psd_file, int &max_it, int &flavors,
                       int &ntau, double &tol, bool &verbose);

double read_tau_from_file(std::vector<std::vector<double>> &G_tau,
                          const std::string &filename, bool verbose);

void dump_tau_to_file(const std::vector<std::vector<double>> &G_tau,
                      double beta, const std::string &filename);

template <typename T>
int hankel_nhalf(const T &Gtau) {
  int N = Gtau.size();
  int Nhalf;
  if (Gtau.size() % 2 == 0) {
    Nhalf = N / 2;
  } else {
    Nhalf = N / 2 - 1;
  }
  return Nhalf;
}

template <typename T>
void vector_to_hankel_full(const T &f, Eigen::MatrixXd &H) {
  if (f.size() % 2 != 1) throw std::invalid_argument("ntau must be odd");
  int Nhalf = f.size() / 2 + 1;
  H = Eigen::MatrixXd::Zero(Nhalf, Nhalf);
  for (int p = 0; p < Nhalf; ++p) {
    for (int q = 0; q < Nhalf; ++q) {
      H(p, q) = f[p + q];
    }
  }
}

template <typename T>
void hankel_to_vector_full(const Eigen::MatrixXd &H, T &f) {
  if (f.size() % 2 != 1) throw std::invalid_argument("ntau must be odd");
  int Nhalf = f.size() / 2 + 1;
  for (int q = 0; q < Nhalf; ++q) {
    f[q] = H(0, q);
    f[q + Nhalf - 1] = H(Nhalf - 1, q);
  }
}
