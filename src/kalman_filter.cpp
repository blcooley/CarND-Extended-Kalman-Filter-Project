#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  I_ = MatrixXd::Identity(x_.size(), x_.size());
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  Eigen::VectorXd y(2);
  y = z - H_ * x_;
  Eigen::MatrixXd Ht = H_.transpose();
  Eigen::MatrixXd S = H_ * P_ * Ht + R_;
  Eigen::MatrixXd K = P_ * Ht * S.inverse();

  x_ = x_ + K * y;
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  Eigen::VectorXd zhat(3);
  double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double phi = atan2(x_(1), x_(0));
  double rhodot;
  if (abs(rho) > 1e-9) {
    rhodot = ( x_(0) * x_(2) + x_(1) * x_(3) ) / rho;
  } else {
    rhodot = 0;
  }

  zhat << rho, phi, rhodot;

  Eigen::VectorXd y(3);
  y = z - zhat;
  while (y(1) > M_PI) {
    y(1) = y(1) - 2 * M_PI;
  }

  while (y(1) < -M_PI) {
    y(1) = y(1) + 2 * M_PI;
  }


  Eigen::MatrixXd Ht = H_.transpose();
  Eigen::MatrixXd S = H_ * P_ * Ht + R_;
  Eigen::MatrixXd K = P_ * Ht * S.inverse();

  x_ = x_ + K * y;
  P_ = (I_ - K * H_) * P_;
}
