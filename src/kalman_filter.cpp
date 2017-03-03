//======================================================================================================================
#include "kalman_filter.h"
//======================================================================================================================

KalmanFilter::KalmanFilter()
{}

//----------------------------------------------------------------------------------------------------------------------

KalmanFilter::~KalmanFilter() {}

//----------------------------------------------------------------------------------------------------------------------

void KalmanFilter::Init(
  VectorXd& x_in,
  MatrixXd& P_in,
  MatrixXd& F_in,
  MatrixXd& H_in,
  MatrixXd& R_in,
  MatrixXd& Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

//----------------------------------------------------------------------------------------------------------------------

void KalmanFilter::Predict()
{
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

//----------------------------------------------------------------------------------------------------------------------

void KalmanFilter::Update(const VectorXd &z)
{
  const VectorXd z_pred = H_*x_;
  const VectorXd y = z - z_pred;
  const MatrixXd Ht = H_.transpose();
  const MatrixXd S = H_ * P_ * Ht + R_;
  const MatrixXd K = P_ * Ht * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  const MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

//----------------------------------------------------------------------------------------------------------------------

void KalmanFilter::UpdateWithAlreadyPredictedMeasurements(const VectorXd& z, const VectorXd& z_pred)
{
  /**
    TODO:
    * update after predicting
  */
}

//======================================================================================================================
