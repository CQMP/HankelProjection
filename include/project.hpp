#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <vector>

void project_to_PSD(Eigen::MatrixXd &H);
void project_to_hankel(Eigen::MatrixXd &H);
Eigen::MatrixXd project_to_density(const Eigen::MatrixXd &H0, double density,
                                   bool verbose);
Eigen::MatrixXd project_general_to_PSD(const Eigen::MatrixXd &H0, bool verbose);
Eigen::MatrixXd project_top_to_PSD(const Eigen::MatrixXd &H0, bool verbose);
Eigen::MatrixXd project_general_to_hankel(const Eigen::MatrixXd &H0,
                                          bool verbose);
void hankel_dyproject(std::vector<double> &G, int max_iter, double density,
                      double tol, bool verbose);
