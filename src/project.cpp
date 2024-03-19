#include "project.hpp"

#include "aux.hpp"

void project_to_PSD(Eigen::MatrixXd &H, bool verbose) {
  Eigen::MatrixXd Horig(H);
  H = 0.5 * (H + H.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ES(H);

  Eigen::VectorXd eigen_values = ES.eigenvalues();
  Eigen::MatrixXd eigen_vectors = ES.eigenvectors();

  double nevsum = 0.;
  for (int i = 0; i < H.rows(); ++i) {
    if (eigen_values[i] < 0) {
      nevsum += eigen_values[i];
      eigen_values[i] = 0;
    }
  }
  H = eigen_vectors * eigen_values.asDiagonal() * eigen_vectors.transpose();
  if (verbose) std::cout << "negative eigenvalues sum: " << nevsum << std::endl;
}

Eigen::MatrixXd project_general_to_PSD(const Eigen::MatrixXd &H0,
                                       bool verbose) {
  if (verbose) std::cout << "--- PSD projection ---" << std::endl;
  Eigen::MatrixXd H(H0);
  project_to_PSD(H, verbose);
  if (verbose) {
    std::cout << "diff norm (PSD): " << (H - H0).norm() << std::endl;
  }
  return H;
}

Eigen::MatrixXd project_top_to_PSD(const Eigen::MatrixXd &H0, bool verbose) {
  if (verbose) std::cout << "--- top submatrix PSD projection ---" << std::endl;
  Eigen::MatrixXd H(H0);
  Eigen::MatrixXd Htop = H.block(0, 1, H.rows() - 1, H.cols() - 1);
  Eigen::MatrixXd H0top = Htop;
  project_to_PSD(Htop, verbose);
  H.block(0, 1, H.rows() - 1, H.cols() - 1) = Htop;
  for (int i = 1; i < H.rows(); ++i) {
    H(i, 0) = H(0, i);
    H(H.rows() - 1, i - 1) = H(i - 1, H.rows() - 1);
  }
  if (verbose) {
    std::cout << "diff norm (PSD top submatrix): " << (H0top - Htop).norm()
              << std::endl;
    std::cout << "diff norm (PSD top): " << (H - H0).norm() << std::endl;
  }
  return H;
}

void project_to_hankel(Eigen::MatrixXd &H) {
  int s = H.rows();
  for (int i = 0; i < 2 * s; ++i) {
    double av = 0.;
    double ct = 0.;
    for (int j = 0; j < i + 1; ++j) {
      if ((i - j >= 0) && (j >= 0) && (i - j < s) && (j < s)) {
        av += H(i - j, j);
        ct++;
      }
    }
    for (int j = 0; j < i + 1; ++j) {
      if ((i - j >= 0) && (j >= 0) && (i - j < s) && (j < s)) {
        H(i - j, j) = av / ct;
      }
    }
  }
}

Eigen::MatrixXd project_general_to_hankel(const Eigen::MatrixXd &H0,
                                          bool verbose) {
  Eigen::MatrixXd H(H0);
  project_to_hankel(H);
  if (verbose) {
    std::cout << "--- hankel projection ---" << std::endl;
    std::cout << "diff norm (hankel): " << (H - H0).norm() << std::endl;
  }
  return H;
}

Eigen::MatrixXd project_to_density(const Eigen::MatrixXd &H0, double density,
                                   bool verbose) {
  Eigen::MatrixXd H(H0);
  H(0, 0) = 1 - density;
  H(H.rows() - 1, H.cols() - 1) = density;
  if (verbose) {
    std::cout << "--- density projection ---" << std::endl;
    std::cout << "diff norm (density): " << (H - H0).norm() << std::endl;
  }
  return H;
}

void hankel_dyproject(std::vector<double> &G, int max_iter, double density,
                      double tol, bool verbose) {
  if (G[0] > 0)
    throw std::runtime_error("expecting many-body convention, G(tau)<0");
  int Ndata = G.size();
  Eigen::VectorXd f(G.size());
  for (int i = 0; i < G.size(); ++i) f[i] = -1 * G[i];  // QMC convention

  Eigen::MatrixXd Hankel;
  vector_to_hankel_full(f, Hankel);
  int s = Hankel.rows();  // matrix size

  // copied from https://github.com/mjhough/Dykstra/blob/main/dykstra/Dykstra.py
  // x = x0.copy()
  Eigen::MatrixXd x = Hankel;
  Eigen::MatrixXd prev_x = Hankel;
  Eigen::MatrixXd x_again = Hankel;
  // p = len(P)
  int p = 4;  // 4 projections
  // y = np.zeros((p,x0.shape[0]))
  std::vector<Eigen::MatrixXd> y =
      std::vector<Eigen::MatrixXd>(p, Eigen::MatrixXd::Zero(s, s));

  // n = 0
  int n = 0;
  // cI = float('inf')
  double cI = std::numeric_limits<double>::infinity();
  // while n < max_iter and cI >= tol:
  while ((n < max_iter) && (cI >= tol)) {
    cI = 0;
    for (int i = 0; i < p; ++i) {
      // Update iterate
      prev_x = x;
      switch (i) {
        case 0:
          x = project_to_density(prev_x - y[i], density, verbose);
          break;
        case 1:
          x = project_general_to_hankel(prev_x - y[i], verbose);
          break;
        case 2:
          x = project_general_to_PSD(prev_x - y[i], verbose);
          break;
        case 3:
          x = project_top_to_PSD(prev_x - y[i], verbose);
          break;
        default:
          throw std::logic_error("invalid projection");
      }
      // # Update increment
      Eigen::MatrixXd prev_y = y[i];
      y[i] = x - (prev_x - prev_y);
      // # Stop condition
      double norm = (prev_y - y[i]).norm();
      cI += norm * norm;
    }
    n++;
    std::cout << "### iteration: " << n << " tolerance: " << cI
              << " correction: " << (prev_x - x).norm() << " ###" << std::endl;
  }

  x_again = project_to_density(x, density, false);
  std::cout << "### DONE ###" << std::endl
            << "diff norm (density): " << (x - x_again).norm() << std::endl;
  x_again = project_general_to_hankel(x, false);
  std::cout << "diff norm (hankel): " << (x - x_again).norm() << std::endl;
  x_again = project_general_to_PSD(x, false);
  std::cout << "diff norm (PSD): " << (x - x_again).norm() << std::endl;
  hankel_to_vector_full(x, f);

  for (int i = 0; i < Ndata; ++i) G[i] = -1 * f[i];  // many-body convention
}
