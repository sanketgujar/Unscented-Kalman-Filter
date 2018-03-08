#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  ///* time when the state is true, in us
  time_us_ = 0.0;

  ///* State dimension
  n_x_ = 5; 

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_x_; 

  ///* Augmented state dimension
  n_aug_ = 7; //2*lambda +1

  ///* Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1 );
 
  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_ + 1 ); //5x11

  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;

  cout<<"Initialized perfectly"<<endl;


}

UKF::~UKF() {} //Deconstructor ---Empty

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
  // cout<<" Init : "<<is_initialized_<<endl;     
  if (!is_initialized_){
    // cout<<"Initialized to do";
    //first measurement and covaraince matrix
    x_ << 1 , 1, 1, 1, 0.1 ;
    P_  << 0.15, 0 , 0 , 0 ,0 ,
          0 ,0.15 , 0,0,0,
          0,0,1,0,0,
          0,0,0,1,0,
          0,0,0,0,1;

    time_us_ = meas_package.timestamp_; 
  

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float ro  = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);
      x_(0) = ro*cos(phi);
      x_(1) = ro*sin(phi);  
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    

    is_initialized_ = true;
    cout <<"First step initialisation";
    return ;
  }//initialisation
  float dt = (meas_package.timestamp_ - time_us_ ) / 1e6;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    UpdateLidar(meas_package);
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    UpdateRadar(meas_package);

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
  //Generate sigma points
  
  //creating sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_ , n_x_*2 + 1);

  //calulating square root of P
  MatrixXd A = P_.llt().matrixL();

  //fisrt col of sigma matrix
  Xsig.col(0) = x_;

  for(int i =0; i < n_x_ ; i ++ ){
    Xsig.col(i+1) = x_ + sqrt(lambda_ + n_x_ )*A.col(i);
    Xsig.col(i+n_x_ +1) =  x_ - sqrt(lambda_ + n_x_ )*A.col(i);
  }

  //Augment Sigma points now
  VectorXd X_aug_ = VectorXd(n_aug_);
  X_aug_.head(5) = x_ ;
  X_aug_(5) = 0;
  X_aug_(6) = 0;

  //Creating augmented covariance matrix
  MatrixXd P_aug_ = MatrixXd(n_aug_ , n_aug_);
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  //sqrt matrix
  MatrixXd L = P_aug_.llt().matrixL();

  // creating augmented sigma points
  MatrixXd Xsig_aug_  = MatrixXd(n_aug_ , n_aug_*2 +1);
  Xsig_aug_.col(0) = X_aug_;
  for (int i = 0 ; i < n_aug_ ; i++ ){
    Xsig_aug_.col(i+1) = X_aug_ + sqrt(lambda_ + n_aug_)*L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = X_aug_ - sqrt(lambda_ + n_aug_)*L.col(i);

  }

  //predicting sigma points..
  for (int i =0 ; i < 2*n_aug_ + 1 ; i ++){
    double p_x_ = Xsig_aug_(0,i);
    double p_y_ = Xsig_aug_(1,i);
    double v_   = Xsig_aug_(2,i);
    double yaw_ = Xsig_aug_(3,i);
    double yawd_= Xsig_aug_(4,i);
    double nu_a_= Xsig_aug_(5,i);
    double nu_yawdd_= Xsig_aug_(6,i);

    //predicting state values.
    double px_p_ , py_p_;

    //avoiding division by zero.
    if (fabs(yawd_) > 0.001){
      px_p_ = p_x_ + (v_/yawd_)*(sin(yaw_ + yawd_*delta_t) - sin(yaw_));
      py_p_ = p_y_ + (v_/yawd_)*( - cos(yaw_ + yawd_*delta_t) + cos(yaw_));
    }
    else{
      px_p_ = p_x_ + v_*delta_t*cos(yaw_);
      py_p_ = p_y_ + v_*delta_t*sin(yaw_);
    }

    double v_p_    = v_;
    double yaw_p_  = yaw_ + yawd_*delta_t;
    double yawd_p_ = yawd_;

    //adding noise....
    px_p_ += 0.5*nu_a_*delta_t*delta_t*cos(yaw_);
    py_p_ += 0.5*nu_a_*delta_t*delta_t*sin(yaw_);
    v_p_  += nu_a_*delta_t;
    yaw_p_+= 0.5*nu_yawdd_*delta_t*delta_t;
    yawd_p_ += nu_yawdd_*delta_t;


    //Now writing sigma point into right colums...
    Xsig_pred_(0,i) = px_p_;
    Xsig_pred_(1,i) = py_p_;
    Xsig_pred_(2,i) = v_p_;
    Xsig_pred_(3,i) = yaw_p_;  
    Xsig_pred_(4,i) = yawd_p_;
  }//for ending

  //Calculating Predicted Mean and Covariance....
  double weight_0 = lambda_ / (lambda_  + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1 ; i < 2*n_aug_ +1 ; i ++ ){
    weights_(i) = 0.5 /(n_aug_ + lambda_);
  }

  //predicting state mean
  x_.fill(0.0);
  for (int i =0 ; i < 2*n_aug_ + 1; i++)
    x_ += weights_(i)*Xsig_pred_.col(i);
  //x is the mean

  //predicting state covariance matrix.....
  P_.fill(0.0);
  for (int i = 0 ; i < 2*n_aug_ + 1; i++){
    VectorXd x_diff = Xsig_pred_.col(i)  - x_;

    //angle normalisation..
    while(x_diff(3) > M_PI) x_diff(3) -= 2*M_PI;
    while(x_diff(3) <-M_PI) x_diff(3) += 2*M_PI;
    P_ += weights_(i)*x_diff*x_diff.transpose();
  }

}//Prediction end......

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

  int n_z_  = 2; //measurment dimension
  MatrixXd Zsig_ = MatrixXd(n_z_ , 2*n_aug_ + 1);  //sigma points...

  VectorXd z_ = VectorXd(n_z_);
  z_ = meas_package.raw_measurements_;
  //transforming sigma points into measurment space....
  for (int i = 0 ; i < 2*n_aug_ + 1 ; i++){
    double p_x_ = Xsig_pred_(0,i);
    double p_y_ = Xsig_pred_(1,i);
    
    //measurment model//
    Zsig_(0,i) = p_x_; //x-position
    Zsig_(1,i) = p_y_; //y-position
  }

  //mean meaurment prediction
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_ + 1 ; i ++)
    z_pred_ += weights_(i)*Zsig_.col(i);


  //Covarince matrix S
  MatrixXd S_ = MatrixXd(n_z_ , n_z_);
  S_.fill(0.0);
  for (int i =0; i < 2*n_aug_ + 1 ; i ++){
    VectorXd z_diff  = Zsig_.col(i) - z_pred_;
    S_ += weights_(i)*z_diff*z_diff.transpose();
  }

  // adding measurment noise covariance matrix

  MatrixXd R_ = MatrixXd(n_z_ , n_z_);
  R_ << std_laspx_*std_laspx_ , 0 ,
        0,std_laspy_*std_laspy_;
  S_ += R_;

  //creating cross corelation matrix...
  MatrixXd Tc_ = MatrixXd(n_x_ , n_z_);
  Tc_.fill(0.0);
  for (int i = 0 ; i < 2*n_aug_ + 1 ; i++){
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    VectorXd x_diff = Xsig_pred_.col(i) - x_ ;
    Tc_ += weights_(i)*x_diff*z_diff.transpose();
  }//end of for

  //kalman gain
  MatrixXd k_ = Tc_*S_.inverse();

  VectorXd z_diff = z_ - z_pred_;

  //Update state and covariance matrix
  x_ +=  k_*z_diff;
  P_ -=  k_*S_*k_.transpose();

  NIS_laser_ = z_diff.transpose()*S_.inverse()*z_diff;
}//end for update lidar

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
  int n_z_  = 3; //measurment dimension
  MatrixXd Zsig_ = MatrixXd(n_z_ , 2*n_aug_ + 1);  //sigma points...

  VectorXd z_ = VectorXd(n_z_);
  z_ = meas_package.raw_measurements_;
  //transforming sigma points into measurment space....
  for (int i = 0 ; i < 2*n_aug_ + 1 ; i++){
    double p_x_ = Xsig_pred_(0,i);
    double p_y_ = Xsig_pred_(1,i);
    double v_   = Xsig_pred_(2,i);
    double yaw_ = Xsig_pred_(3,i);
    
    double vx_  = cos(yaw_)*v_;  
    double vy_  = sin(yaw_)*v_;

    //measurment model//
    Zsig_(0,i) = sqrt(p_x_*p_x_ + p_y_*p_y_); //r
    Zsig_(1,i) =  atan2(p_y_, p_x_); //phi
    Zsig_(2,i) = (p_x_*vx_ + p_y_*vy_) / sqrt(p_x_*p_x_ + p_y_*p_y_); //rho.  
  }

  //mean meaurment prediction
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_ + 1 ; i ++)
    z_pred_ += weights_(i)*Zsig_.col(i);


  //Covarince matrix S
  MatrixXd S_ = MatrixXd(n_z_ , n_z_);
  S_.fill(0.0);
  for (int i =0; i < 2*n_aug_ + 1 ; i ++){
    VectorXd z_diff  = Zsig_.col(i) - z_pred_;
    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) <-M_PI) z_diff(1) += 2.*M_PI; 
    S_ += weights_(i)*z_diff*z_diff.transpose();
  }

  // adding measurment noise covariance matrix

  MatrixXd R_ = MatrixXd(n_z_ , n_z_);
  R_ << std_radr_*std_radr_ , 0 ,0 ,
        0, std_radphi_*std_radphi_ , 0 ,
        0 ,0 ,std_radrd_*std_radrd_;
  S_ += R_;

  //creating cross corelation matrix...
  MatrixXd Tc_ = MatrixXd(n_x_ , n_z_);
  Tc_.fill(0.0);
  for (int i = 0 ; i < 2*n_aug_ + 1 ; i++){
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) <-M_PI) z_diff(1) += 2.*M_PI;
    VectorXd x_diff = Xsig_pred_.col(i) - x_ ;
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) <-M_PI) x_diff(3) += 2.*M_PI;
    Tc_ += weights_(i)*x_diff*z_diff.transpose();
  }//end of for

  //kalman gain
  MatrixXd k_ = Tc_*S_.inverse();

  VectorXd z_diff = z_ - z_pred_;
  while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while(z_diff(1) <-M_PI) z_diff(1) += 2.*M_PI;

  //Update state and covariance matrix
  x_ +=  k_*z_diff;
  P_ -=  k_*S_*k_.transpose();

  NIS_radar_ = z_diff.transpose()*S_.inverse()*z_diff;
}//Radar update end

