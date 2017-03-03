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

  const float s = 0.0225;
  const float sigma_px = s;
  const float sigma_py = s;
  R_laser_ << sigma_px,  0.0,  0.0, sigma_py;

  const float sigma_rho = s;
  const float sigma_phi = s;
  const float sigma_rho_dot = s;
  R_radar_ << sigma_rho, 0.0, 0.0,
    0.0, sigma_phi, 0.0,
    0.0, 0.0, sigma_rho_dot;
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

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::Initialize(const MeasurementPackage& m)
{
  // TODO:
  //   * Initialize the state ekf_.x_ with the first measurement.
  //   * Create the covariance matrix.
  //   * Remember: you'll need to convert radar from polar to cartesian coordinates.

  // first measurement
  cout << "EKF: " << endl;
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

  if (dt < 0.001f)
  {
    // Two simultaneous measurements -> skip prediction
    return;
  }

  previous_timestamp_ = current_time_stamp;

  UpdateProcessCovarianceMatrix(dt, 125, 125);
  UpdateStateTransitionMatrix(dt);

  ekf_.Predict();
}

//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::Update(const MeasurementPackage& m)
{
  if (m.sensor_type_ == MeasurementPackage::RADAR)
  {
    Hj_ = Tools::CalculateJacobian(ekf_.x_);
    ekf_.Update(m.raw_measurements_, Hj_, R_radar_);
  }
  else
  {
    ekf_.Update(m.raw_measurements_, H_laser_, R_laser_);
  }
}

//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::UpdateProcessCovarianceMatrix(float dt, float noise_ax, float noise_ay)
{
  const float dt_2 = dt * dt;
  const float dt_3 = dt_2 * dt;
  const float dt_4 = dt_3 * dt;

  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
    0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
    dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
    0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
}

//----------------------------------------------------------------------------------------------------------------------

void FusionEKF::UpdateStateTransitionMatrix(float dt)
{
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
}

//======================================================================================================================
