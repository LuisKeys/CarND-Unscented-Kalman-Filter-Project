#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  ///* State dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  ///* Augmented state dimension
  n_aug_ = 7;

  n_2_aug_plus_1_ = 2 * n_aug_ + 1;

  ///* predicted sigma points matrix

  Xsig_pred_ = MatrixXd(n_x_, n_2_aug_plus_1_);
  Xsig_pred_.setZero();

  ///* Weights of sigma points
  weights_ = VectorXd(n_2_aug_plus_1_);;
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_; 

  double std_radr_2 = std_radr_ * std_radr_;
  double std_radphi_2 = std_radphi_ * std_radphi_;
  double std_radrd_2 = std_radrd_ * std_radrd_;

  R_radar_ = MatrixXd(n_radar_, n_radar_);
  R_radar_ << std_radr_2, 0, 0,
              0, std_radphi_2, 0,
              0, 0, std_radrd_2;

  double std_laspx_2 = std_laspx_ * std_laspx_;
  double std_laspy_2 = std_laspy_ * std_laspy_;

  R_lidar_ = MatrixXd(n_lidar_, n_lidar_);
  R_lidar_ << std_laspx_2, 0,
              0, std_laspy_2;  

  sqrt_lamb_plus_n_aug_ = sqrt(lambda_ + n_aug_);

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * Angle normalization
 * @param {double*} angle The angle to be normalized to [-Pi, Pi] range
 * either radar or laser.
 */
void UKF::Anglenormalization(double *angle) {
    while (*angle > M_PI) * angle -= 2. * M_PI;
    while (*angle < -M_PI) * angle += 2. * M_PI;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::InitFilter(MeasurementPackage meas_package) {

    //Radar
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      //Distance
      double rho = meas_package.raw_measurements_[0];
      //Angle
      double phi = meas_package.raw_measurements_[1];
      //Velocity of angle
      double rho_dot = meas_package.raw_measurements_[2];
      // Coordinates convertion from polar to cartesian
      double px = rho * cos(phi); 
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);      
      double vx_2 = vx * vx;
      double vy_2 = vy * vy;
      double v  = sqrt(vx_2 + vy_2);

      x_ << px, py, v, 0, 0;
    }
    //Lidar
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 
            0, 0, 0;

      if (fabs(x_(0)) < small_number)
        x_(0) = 0.01;

      if (fabs(x_(1)) < small_number)
        x_(1) = 0.01;
    }

    //Initi state covariance with Identity matrix
    P_ = MatrixXd::Identity(n_x_, n_x_);
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  //Init if not yet initialized
  if(!is_initialized_) {
    //Initialize with first meassure
    //If first meassure is Radar
    InitFilter(meas_package);
    previous_timestamp_ = meas_package.timestamp_;    
    is_initialized_ = true;
    return;
  }

  //dt - expressed in seconds
  dt_ = ((float)(meas_package.timestamp_ - previous_timestamp_)) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;  

  Prediction(dt_);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  //If first meassure is Lidar
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //Sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_2_aug_plus_1_);
  VectorXd x_aug = VectorXd(n_aug_); 
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  x_aug.setZero();
  P_aug.setZero();  

  x_aug.head(n_x_) = x_;
  P_aug.topLeftCorner(n_x_, n_x_) = P_;

  MatrixXd Q = MatrixXd(2, 2);
  Q.setZero();
  const float std_a_2 = std_a_ * std_a_;
  const float std_yawdd_2 = std_yawdd_ * std_yawdd_; 
  Q(0, 0) = std_a_2;
  Q(1, 1) = std_yawdd_2;
  P_aug.bottomRightCorner(2, 2) = Q;

  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) << x_aug; 

  VectorXd sqrt_lambda_n_aug_A;

  for(int i = 0; i < n_aug_; i++) {
    sqrt_lambda_n_aug_A = sqrt_lamb_plus_n_aug_ * A.col(i);
    Xsig_aug.col(i + 1) = x_aug + sqrt_lambda_n_aug_A;
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt_lambda_n_aug_A;
  }
  
  // Predict sigma points
  for (int i = 0; i < n_2_aug_plus_1_; i++) {
    // Extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    double sin_yaw = sin(yaw);
    double cos_yaw = cos(yaw);
    double arg = yaw + yawd * delta_t;
   
    //Predicted state values
    double px_p, py_p;

    //Manage division by zero
    if (fabs(yawd) > small_number) { 
        double v_yawd = v / yawd;
        px_p = p_x + v_yawd * (sin(arg) - sin_yaw);
        py_p = p_y + v_yawd * (cos_yaw - cos(arg) );
    }
    else {
        double v_delta_t = v * delta_t;
        px_p = p_x + v_delta_t * cos_yaw;
        py_p = p_y + v_delta_t * sin_yaw;
    }

    double v_p = v;
    double yaw_p = arg;
    double yawd_p = yawd;

    //Add noise
    double coef_05_nu_a_delta_t = 0.5 * nu_a * delta_t;
    px_p += coef_05_nu_a_delta_t * cos_yaw;
    py_p += coef_05_nu_a_delta_t * sin_yaw;
    v_p += nu_a * delta_t;
    yaw_p += 0.5 * nu_yawdd * delta_t;
    yawd_p += nu_yawdd * delta_t;

    //Load predicted sigma points
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  //Predicted state mean
  x_ = Xsig_pred_ * weights_;

  //Predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_2_aug_plus_1_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Anglenormalization(&(x_diff(3)));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // Set measurement dimension
  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_lidar_, n_2_aug_plus_1_);
  UpdateUKF(meas_package, n_lidar_, Zsig);}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // Radar measure: r, phi, and r_dot
  //Matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_radar_, n_2_aug_plus_1_);
  //Transform sigma points into measurement space
  for (int i = 0; i < n_2_aug_plus_1_; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;
    //Measurement model
    double p_x_2 = p_x * p_x;
    double p_y_2 = p_y * p_y;

    double rho = sqrt(p_x_2 + p_y_2);
    double phi = atan2(p_y, p_x);
    double rho_dot = (p_x * v1 + p_y * v2 ) / rho;
    Zsig(0,i) = rho;
    Zsig(1,i) = phi;
    Zsig(2,i) = rho_dot;
 }

  UpdateUKF(meas_package, n_radar_, Zsig);  
}

  /**
   * Common update method
   * @param meas_package The measurement at k+1
   * @param n_z Device dimensions
   * @param Zsig Sigma points in meassure space
   */
void UKF::UpdateUKF(MeasurementPackage meas_package, int n_z, MatrixXd Zsig){

  //Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred  = Zsig * weights_;
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.setZero();

  for (int i = 0; i < n_2_aug_plus_1_; i++) { 
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // Angle normalization
    Anglenormalization(&(z_diff(1)));
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  //Measurement noise covariance
  MatrixXd R = MatrixXd(n_z, n_z);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    R = R_radar_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){ // Lidar
    R = R_lidar_;
  }

  S += R;
  
  //Create and calculation cross correlation Tc matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_2_aug_plus_1_; i++) { 
    //Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
     //Radar
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      Anglenormalization(&(z_diff(1)));
    }
    //State diff
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Anglenormalization(&(x_diff(3)));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  VectorXd z = meas_package.raw_measurements_;
  //Kalman gain;
  MatrixXd K = Tc * S.inverse();
  //Residual
  VectorXd z_diff = z - z_pred;
  //Radar
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    Anglenormalization(&(z_diff(1)));
  }
  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  //NIS
  //Radar
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    NIS_radar_ = z.transpose() * S.inverse() * z;
  }
  //Lidar
  else if (meas_package.sensor_type_  == MeasurementPackage::LASER) {
    NIS_lidar_ = z.transpose() * S.inverse() * z;
  }
}