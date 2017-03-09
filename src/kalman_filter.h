//======================================================================================================================
#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"
//======================================================================================================================
using Eigen::MatrixXd;
using Eigen::VectorXd;
//======================================================================================================================

class KalmanFilter
{
public:
  KalmanFilter();
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   * @param Q_in Process covariance matrix
   */
  void Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &Q_in);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict();

  /**
   * Updates the state and
   * @param z The measurement at k+1
   * @param H measurement matrix
   * @param R measurement covariance matrix
   */
  void Update(const VectorXd& z, const MatrixXd& H, const MatrixXd& R);

  void UpdateWithAlreadyPredictedMeasurements(
    const VectorXd& z,
    const VectorXd& z_pred,
    const MatrixXd& H,
    const MatrixXd& R);

public:
  // state vector
  VectorXd x_;

  // state covariance matrix
  MatrixXd P_;

  // state transition matrix
  MatrixXd F_;

  // process covariance matrix
  MatrixXd Q_;

  // identity matrix
  MatrixXd I_;
};

//======================================================================================================================
#endif /* KALMAN_FILTER_H_ */
//======================================================================================================================
