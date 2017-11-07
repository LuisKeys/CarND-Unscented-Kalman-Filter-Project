#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  const double small_number = 0.00001;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Previous time stamp to get delta time
  long previous_timestamp_;

  ///* Delta Time
  float dt_;

  ///* Augmented columns
  int n_2_aug_plus_1_;

  ///* Current NIS for radar
  double NIS_radar_;

  ///* Current NIS for laser
  double NIS_lidar_;
    /**
   * Constructor
   */

  ///* Radar dimensions
  int n_radar_ = 3;

  ///* Lidar dimensions
  int n_lidar_ = 2;

  ///* Radar measurement noise covariance matrix
  MatrixXd R_radar_;
  
  ///* Lidar measurement noise covariance matrix
  MatrixXd R_lidar_;

  ///* Square root of Lamda + n_aug
  double sqrt_lamb_plus_n_aug_;
    
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  void Anglenormalization(double *angle);

  /**
   * Init
   * @param meas_package The latest measurement data of either radar or laser
   */
  void InitFilter(MeasurementPackage meas_package);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * UpdateUKF
   * @param meas_package The measurement at k+1
   * @param n_z Device dimensions
   * @param Zsig Sigma points in meassure space
   */
  void UpdateUKF(MeasurementPackage meas_package, int n_z, MatrixXd Zsig);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
