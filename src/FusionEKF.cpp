//======================================================================================================================
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
//======================================================================================================================
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
//======================================================================================================================

FusionEKF::FusionEKF()
: is_initialized_(false)
, previous_timestamp_(0)
, R_laser_(MatrixXd(2, 2))
, R_radar_(MatrixXd(3, 3))
, H_laser_(MatrixXd(2, 4))
, Hj_(MatrixXd(3, 4))
{
  H_laser_ << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0;

  // values taken from project suggestion
  const float var_px = 0.0225;
  const float var_py = 0.0225;
  R_laser_ << var_px,  0.0,  0.0, var_py;

  // values taken from project suggestion
  const float var_rho = 0.09;
  const float var_phi = 0.0009;
  const float var_rho_dot = 0.09;
  R_radar_ << var_rho, 0.0, 0.0,
    0.0, var_phi, 0.0,
    0.0, 0.0, var_rho_dot;
}

//----------------------------------------------------------------------------------------------------------------------

FusionEKF::~FusionEKF() {}

//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::ProcessMeasurement(const MeasurementPackage& m)
{
  if (!is_initialized_)
  {
    is_initialized_ = true;
    Initialize(m);
    return;
  }

  Predict(m);
  Update(m);
}

//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::Initialize(const MeasurementPackage& m)
{
  // first measurement
  auto x = VectorXd(4);

  auto P = MatrixXd(4,4);
  P << 1.0, 0.0, 0.0, 0.0,
       0.0, 1.0, 0.0, 0.0,
       0.0, 0.0, 1000.0, 0.0,
       0.0, 0.0, 0.0, 1000.0;

  auto Q = MatrixXd(4,4);

  auto F = MatrixXd(4,4);
  F << 1.0, 0.0, 1.0, 0.0,
       0.0, 1.0, 0.0, 1.0,
       0.0, 0.0, 1.0, 0.0,
       0.0, 0.0, 0.0, 1.0;

  previous_timestamp_ = m.timestamp_;

  if (m.sensor_type_ == MeasurementPackage::RADAR)
  {
    const double rho = m.raw_measurements_(0,0);
    const double phi = m.raw_measurements_(1,0);
    const double rho_dot = m.raw_measurements_(2,0);
    auto pos = Tools::PolarToCartesian(rho, phi);
    auto vel = Tools::PolarToCartesian(rho_dot, phi);
    x << pos(0,0), pos(1,0), vel(0,0), vel(1,0);
  }
  else if (m.sensor_type_ == MeasurementPackage::LASER)
  {
    x <<  m.raw_measurements_(0,0), m.raw_measurements_(1,0), 0.0, 0.0;
  }

  ekf_.Init(x, P, F, Q);
}

//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::Predict(const MeasurementPackage& m)
{
  const long current_time_stamp = m.timestamp_;
  const float dt = (current_time_stamp - previous_timestamp_) / 1.0e6;

  if (dt < 0.005f)
  {
    // Second of two (nearly) simultaneous measurements -> skip prediction.
    return;
  }

  previous_timestamp_ = current_time_stamp;

  UpdateProcessCovarianceMatrix(dt);
  UpdateStateTransitionMatrix(dt);

  ekf_.Predict();
}

//----------------------------------------------------------------------------------------------------------------------

VectorXd FusionEKF::PredictRadarMeasurement(const VectorXd& x) const
{
  const float px = x(0,0);
  const float py = x(1,0);
  const float vx = x(2,0);
  const float vy = x(3,0);
  const float eps = 1e-5;
  const float rho = sqrt(px * px + py * py);
  const float phi = atan2(py, px);
  const float rho_dot = (px * vx + py * vy) / (eps + rho);
  VectorXd result(3);
  result << rho, phi, rho_dot;
  return result;
}


//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::Update(const MeasurementPackage& m)
{
  if (m.sensor_type_ == MeasurementPackage::RADAR)
  {
    const auto& z = m.raw_measurements_;
    Hj_ = Tools::CalculateJacobian(ekf_.x_);
    const VectorXd z_pred = PredictRadarMeasurement(ekf_.x_);
    ekf_.UpdateWithAlreadyPredictedMeasurements(z, z_pred, Hj_, R_radar_);
  }
  else
  {
    ekf_.Update(m.raw_measurements_, H_laser_, R_laser_);
  }
}

//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::UpdateProcessCovarianceMatrix(float dt)
{
  const float dt_2 = dt * dt;
  const float dt_3 = dt_2 * dt;
  const float dt_4 = dt_3 * dt;
  const float var_a = 9; // as suggested by the project.

  ekf_.Q_ <<  dt_4/4*var_a, 0, dt_3/2*var_a, 0,
              0, dt_4/4*var_a, 0, dt_3/2*var_a,
              dt_3/2*var_a, 0, dt_2*var_a, 0,
              0, dt_3/2*var_a, 0, dt_2*var_a;
}

//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::UpdateStateTransitionMatrix(float dt)
{
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
}

//======================================================================================================================
