//======================================================================================================================
#ifndef FusionEKF_H_
#define FusionEKF_H_
//======================================================================================================================
#include "measurement_package.h"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
//======================================================================================================================

class FusionEKF {
public:
  FusionEKF();
  virtual ~FusionEKF();

  // Run the whole flow of the Kalman Filter from here.
  void ProcessMeasurement(const MeasurementPackage& measurement_pack);

  // Kalman Filter update and prediction math lives in here.
  KalmanFilter ekf_;

private:
  void Initialize(const MeasurementPackage& measurement_package);
  void Predict(const MeasurementPackage& measurement_package);
  void Update(const MeasurementPackage& measurement_package);
  void UpdateProcessCovarianceMatrix(float dt, float noise_ax, float noise_ay);
  void UpdateStateTransitionMatrix(float dt);

private:
  bool is_initialized_;
  long previous_timestamp_;

  MatrixXd R_laser_;
  MatrixXd R_radar_;
  MatrixXd H_laser_;
  MatrixXd Hj_;
};

//======================================================================================================================
#endif /* FusionEKF_H_ */
//======================================================================================================================
