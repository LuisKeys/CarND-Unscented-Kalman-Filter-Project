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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  ///* Augmented state dimension
  n_aug_ = 7;

  n_2_aug_plus_1_ = 2 * n_aug_ + 1;

  ///* predicted sigma points matrix

  MatrixXd Xsig_pred_ = MatrixXd(n_x_, n_2_aug_plus_1_);
  Xsig_pred_.setZero();

  ///* Weights of sigma points
  VectorXd weights_ = VectorXd(n_2_aug_plus_1_);;
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_; 

  is_initialized_ = false;

}

UKF::~UKF() {}

/**
 * Angle normalization
 * @param {double*} angle The angle to be normalized
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

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float px = ro * cos(theta);
      float py = ro * sin(theta);
      x_ << px, py, 
            0, 0, 0;
    }
    //If first meassure is Lidar
    else {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 
            0, 0, 0;
    }

    //Initi state covariance with Identity matrix
    P_ = MatrixXd::Identity(n_x_, n_x_);
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  //Init
  //Normalize values that may create a calculation issue
  if(meas_package.raw_measurements_.squaredNorm() < 0.000001) {
    meas_package.raw_measurements_.fill(0.001);
  }


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
  else {
    UpdateLidar(meas_package);
  }

}

/**
 * Create sigma points
 * @param {MatrixXd*} Xsig_out the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  VectorXd x_aug = VectorXd(n_aug_); 
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.setZero();  
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_2_aug_plus_1_);
  x_aug.setZero();
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

  const float sqrt_lamb_plus_n_aug = sqrt(lambda_ + n_aug_);

  for (int i = 1; i < n_aug_ + 1; ++i) {
    Xsig_aug.col(i) << x_aug + sqrt_lamb_plus_n_aug * A.col(i - 1);
  }

  const int upper_limit = n_aug_ + 1;
  for (int i = upper_limit; i < n_2_aug_plus_1_; ++i) {
    Xsig_aug.col(i) << x_aug - sqrt_lamb_plus_n_aug * A.col(i - upper_limit);
  }

  *Xsig_out = Xsig_aug;
}

/**
 * Predicts sigma points
 * @param Xsig_aug Augmented sigma points buffer
 * @param Xsig_pred_out Predicted sigma points buffer
 */
void UKF::PredictSigmaPoints(MatrixXd Xsig_aug, MatrixXd* Xsig_pred_out, double delta_t) {
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, n_2_aug_plus_1_);

  /*******************************************************************************
  * Student part begin
  ******************************************************************************/

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column

  //predict sigma points
  for (int i = 0; i < n_2_aug_plus_1_; i++) {
    //extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    //predicted state values
    double p_x_p, p_y_p;
    double yawd_delta_t = yawd * delta_t;      
    double sin_yaw = sin(yaw);
    double cos_yaw = cos(yaw);

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      p_x_p = p_x + v / yawd * (sin(yaw + yawd_delta_t) - sin_yaw);
      p_y_p = p_y + v / yawd * (cos_yaw - cos(yaw + yawd_delta_t));
    }
    else {
      p_x_p = p_x + v * delta_t * cos_yaw;
      p_y_p = p_y + v * delta_t * sin_yaw;
    }

    double v_p = v;
    double yaw_p = yaw + yawd_delta_t;
    double yawd_p = yawd;

    //add noise
    double half_nu_a_delt_2 = 0.5 * nu_a * delta_t * delta_t;
    p_x_p = p_x_p + half_nu_a_delt_2 * cos_yaw;
    p_y_p = p_y_p + half_nu_a_delt_2 * sin_yaw;
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0, i) = p_x_p;
    Xsig_pred(1, i) = p_y_p;
    Xsig_pred(2, i) = v_p;
    Xsig_pred(3, i) = yaw_p;
    Xsig_pred(4, i) = yawd_p;
  }

  //write result
  *Xsig_pred_out = Xsig_pred;  
}


/**
 * Predicts sigma points
 * @param Xsig_aug Augmented sigma points buffer
 * @param Xsig_pred_out Predicted sigma points buffer
 */
void UKF::PredictMeanAndCovariance() {

  ////predict state mean
  for (int i = 0; i < n_2_aug_plus_1_; ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  P_.setZero();

  //predict state covariance matrix
  for (int i = 0; i < n_2_aug_plus_1_; ++i) {
    P_ += weights_[i] * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //Sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_x_, n_2_aug_plus_1_);
  AugmentedSigmaPoints(&Xsig_aug);
  
  //Predict sigma points
  PredictSigmaPoints(Xsig_aug, &Xsig_pred_, delta_t);

  std::cout << delta_t << endl;
  std::cout << Xsig_pred_ << endl;

  //Predict State and State covariance
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // Set measurement dimension
  int n_z = 2;
  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_2_aug_plus_1_);
  UpdateUKF(meas_package, n_z, Zsig);}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // Radar measure: r, phi, and r_dot
  int n_z = 3;

  //Matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_2_aug_plus_1_);
  //Transform sigma points into measurement space
  for (int i = 0; i < n_2_aug_plus_1_; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
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

  UpdateUKF(meas_package, n_z, Zsig);  
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
  else {
    NIS_lidar_ = z.transpose() * S.inverse() * z;
  }
}