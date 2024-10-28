#include "trajectory_generator_waypoint.h"
#include <fstream>
#include <iostream>
#include <ros/console.h>
#include <ros/ros.h>
#include <stdio.h>//PolyQPGeneration//selfadd ???
#include <string>

using namespace std;
using namespace Eigen;

#define inf 1 >> 30

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint() {}

TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint() {}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
    const int d_order,           // the order of derivative
    const Eigen::MatrixXd &Path, // waypoints coordinates (3d) //selfadd:MatrixXd path(int(grid_path.size()), 3); path.row(k) = grid_path[k]; 
    const Eigen::MatrixXd &Vel,  // boundary velocity
    const Eigen::MatrixXd &Acc,  // boundary acceleration
    const Eigen::VectorXd &Time) // time allocation in each segment
{
  cout << "[Debug] inside PolyQPGeneration, just start" << endl; // for test
  // enforce initial and final velocity and accleration, for higher order
  // derivatives, just assume them be 0;
  int p_order = 2 * d_order - 1; // the order of polynomial
  int p_num1d = p_order + 1;     // the number of variables in each segment
  cout << "[Debug] inside PolyQPGeneration, 111" << endl; // for test
  int m = Time.size();
  cout << "[Debug] inside PolyQPGeneration, 222" << endl; // for test
  MatrixXd PolyCoeff(m, 3 * p_num1d); //每行是一段polynomial，一行内，三个三个一组(x，y，z)，一共p_num1d组

  /**
   *
   * STEP 3.2:  generate a minimum-jerk piecewise monomial polynomial-based
   * trajectory
   *
   * **/
  //get C,A,Q
  cout << "[Debug] inside PolyQPGeneration, 333" << endl; // for test
  Matrix3d I_3x3 = Matrix3d::Identity();
  MatrixXd I_12x12 = MatrixXd::Identity(12,12);
  cout << "[Debug] inside PolyQPGeneration, 444" << endl; // for test
  MatrixXd Q_Matrix = MatrixXd::Zero(3*p_num1d*m,3*p_num1d*m);// define size? //24mx24m
  MatrixXd A_Matrix = MatrixXd::Zero(3*p_num1d*m,3*p_num1d*m); //24mx24m
  MatrixXd C_Matrix_Trans = MatrixXd::Zero(24*m,12*m+12); //24mx12m+12
  cout << "[Debug] inside PolyQPGeneration, 555" << endl; // for test
  C_Matrix_Trans.block(0,0,3*4,3*4) = I_12x12;
  C_Matrix_Trans.block(24*m-24+12,3*(m+3),3*4,3*4) = I_12x12;
  cout << "[Debug] inside PolyQPGeneration, 666" << endl; // for test
  for(int i = 0;i<m;i++){
  double tm = Time(i);
  cout << "[Debug] inside polynomial trajectory, constructing Qi" << endl; // for test
  //Construct Qi:
  MatrixXd Q_i = MatrixXd::Zero(3 * p_num1d, 3 * p_num1d);
  Q_i.block(0,0,3,3) = 705600*pow(tm,6)*I_3x3;
  Q_i.block(0,3,3,3) = 302400*pow(tm,5)*I_3x3;
  Q_i.block(0,6,3,3) = 100800*pow(tm,4)*I_3x3;
  Q_i.block(0,9,3,3) = 20160*pow(tm,3)*I_3x3;

  Q_i.block(3,0,3,3) = 302400*pow(tm,5)*I_3x3;
  Q_i.block(3,3,3,3) = 129600*pow(tm,4)*I_3x3;
  Q_i.block(3,6,3,3) = 43200*pow(tm,3)*I_3x3;
  Q_i.block(3,9,3,3) = 8640*pow(tm,2)*I_3x3;

  Q_i.block(6,0,3,3) = 100800*pow(tm,4)*I_3x3;
  Q_i.block(6,3,3,3) = 43200*pow(tm,3)*I_3x3;
  Q_i.block(6,6,3,3) = 14400*pow(tm,2)*I_3x3;
  Q_i.block(6,9,3,3) = 2880*pow(tm,1)*I_3x3;

  Q_i.block(9,0,3,3) = 20160*pow(tm,3)*I_3x3;
  Q_i.block(9,3,3,3) = 8640*pow(tm,2)*I_3x3;
  Q_i.block(9,6,3,3) = 2880*pow(tm,1)*I_3x3;
  Q_i.block(9,9,3,3) = 576*I_3x3;
  // Q_i <<  705600*pow(Time(i),6)*I_3x3, 302400*pow(Time(i),5)*I_3x3, 100800*pow(Time(i),4)*I_3x3, 20160*pow(Time(i),3)*I_3x3,
  //         302400*pow(Time(i),5)*I_3x3, 129600*pow(Time(i),4)*I_3x3, 43200*pow(Time(i),3)*I_3x3, 8640*pow(Time(i),2)*I_3x3,
  //         100800*pow(Time(i),4)*I_3x3, 43200*pow(Time(i),3)*I_3x3, 14400*pow(Time(i),2)*I_3x3, 2880*pow(Time(i),1)*I_3x3,
  //         20160*pow(Time(i),3)*I_3x3,  8640*pow(Time(i),2)*I_3x3,  2880*pow(Time(i),1)*I_3x3,  576;
  cout << "[Debug] inside polynomial trajectory, constructing Q_Matrix" << endl; // for test
  Q_Matrix.block(3*p_num1d*i, 3*p_num1d*i, 3*p_num1d, 3*p_num1d) = Q_i;
  cout << "[Debug] inside polynomial trajectory, constructing Ai" << endl; // for test
  //Construct Ai:
  MatrixXd A_i = MatrixXd::Zero(3 * p_num1d, 3 * p_num1d);
  // A_i.block(0,3 * (p_num1d-1),3,3) = 1*I_3x3;
  // A_i.block(3,3 * (p_num1d-2),3,3) = 1*I_3x3;
  // A_i.block(6,3 * (p_num1d-3),3,3) = 2*I_3x3;
  // A_i.block(9,3 * (p_num1d-4),3,3) = 6*I_3x3;

  // A_i.block(12,3 * (p_num1d-1),3,3) = 1*I_3x3;

  A_i <<  1*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3,
          0*I_3x3, 1*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3,
          0*I_3x3, 0*I_3x3, 2*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3,
          0*I_3x3, 0*I_3x3, 0*I_3x3, 6*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3,
          1*I_3x3, pow(tm,1)*I_3x3, pow(tm,2)*I_3x3, pow(tm,3)*I_3x3, pow(tm,4)*I_3x3, pow(tm,5)*I_3x3, pow(tm,6)*I_3x3, pow(tm,7)*I_3x3,
          0*I_3x3, pow(tm,0)*I_3x3, 2*pow(tm,1)*I_3x3, 3*pow(tm,2)*I_3x3, 4*pow(tm,3)*I_3x3, 5*pow(tm,4)*I_3x3, 6*pow(tm,5)*I_3x3, 7*pow(tm,6)*I_3x3,
          0*I_3x3, 0*I_3x3, 2*pow(tm,0)*I_3x3, 6*pow(tm,1)*I_3x3, 12*pow(tm,2)*I_3x3, 20*pow(tm,3)*I_3x3, 30*pow(tm,4)*I_3x3, 42*pow(tm,5)*I_3x3,
          0*I_3x3, 0*I_3x3, 0*I_3x3, 6*pow(tm,0)*I_3x3, 24*pow(tm,1)*I_3x3, 60*pow(tm,2)*I_3x3, 120*pow(tm,3)*I_3x3, 210*pow(tm,4)*I_3x3;
          // pow(tm,7)*I_3x3, pow(tm,6)*I_3x3, pow(tm,5)*I_3x3, pow(tm,4)*I_3x3, pow(tm,3)*I_3x3, pow(tm,2)*I_3x3, pow(tm,1)*I_3x3, 1*I_3x3,
          // 7*pow(tm,6)*I_3x3, 6*pow(tm,5)*I_3x3, 5*pow(tm,4)*I_3x3, 4*pow(tm,3)*I_3x3, 3*pow(tm,2)*I_3x3, 2*pow(tm,1)*I_3x3, pow(tm,0)*I_3x3, 0*I_3x3,
          // 42*pow(tm,5)*I_3x3, 30*pow(tm,4)*I_3x3, 20*pow(tm,3)*I_3x3, 12*pow(tm,2)*I_3x3, 6*pow(tm,1)*I_3x3, 2*pow(tm,0)*I_3x3, 0*I_3x3, 0*I_3x3,
          // 210*pow(tm,4)*I_3x3, 120*pow(tm,3)*I_3x3, 60*pow(tm,2)*I_3x3, 24*pow(tm,1)*I_3x3, 6*pow(tm,0)*I_3x3, 0*I_3x3, 0*I_3x3, 0*I_3x3;
  cout << "[Debug] inside polynomial trajectory, constructing A_Matrix" << endl; // for test
  A_Matrix.block(3*p_num1d*i, 3*p_num1d*i, 3*p_num1d, 3*p_num1d) = A_i; 
  
  cout << "[Debug] inside polynomial trajectory, constructing C_Matrix"<<"m:"<<m<<"i:"<<i << endl; // for test
  //Construct C_Matrix: //24mx12m+12
  if(i>0 && i<m){
    C_Matrix_Trans.block(24*i+0+ 0*3 +0, 12+ 3*(i-1) +0, 3, 3) = I_3x3; // d_{i,0}^{0},x/y/z
    C_Matrix_Trans.block(24*i+0+ 1*3 +0, 3*(m+7)+9*(i-1)+ 0*3 +0, 3, 3) = I_3x3; // d_{i,0}^{1},x/y/z
    C_Matrix_Trans.block(24*i+0+ 2*3 +0, 3*(m+7)+9*(i-1)+ 1*3 +0, 3, 3) = I_3x3; // d_{i,0}^{2},x/y/z
    C_Matrix_Trans.block(24*i+0+ 3*3 +0, 3*(m+7)+9*(i-1)+ 2*3 +0, 3, 3) = I_3x3; // d_{i,0}^{3},x/y/z
    // C_Matrix.block(24*i+0+ 4*3 +0, 3*(m+7)+9*(i-1)+ 4*3 +0, 3, 3) = I_3x3; // d_{i,0}^{4},x/y/z
  }
  if(i<m-1){
    C_Matrix_Trans.block(24*i+12+ 0*3 +0, 12+ 3*i +0, 3, 3) = I_3x3; // d_{i,T}^{0},x/y/z
    C_Matrix_Trans.block(24*i+12+ 1*3 +0, 3*(m+7)+9*i+ 0*3 +0, 3, 3) = I_3x3; // d_{i,T}^{1},x/y/z
    C_Matrix_Trans.block(24*i+12+ 2*3 +0, 3*(m+7)+9*i+ 1*3 +0, 3, 3) = I_3x3; // d_{i,T}^{2},x/y/z
    C_Matrix_Trans.block(24*i+12+ 3*3 +0, 3*(m+7)+9*i+ 2*3 +0, 3, 3) = I_3x3; // d_{i,T}^{3},x/y/z
    // C_Matrix.block(24*i+12+ 4*3 +0, 3*(m+7)+9*i+ 4*3 +0, 3, 3) = I_3x3; // d_{i,T}^{4},x/y/z
  }
  }
  //C A^T Q A^-1 C^T = R
  cout << "[Debug] inside polynomial trajectory, constructing C_Trans" << endl; // for test
  MatrixXd C_Matrix = C_Matrix_Trans.transpose();
  cout << "[Debug] inside polynomial trajectory, constructing A_inv" << endl; // for test
  MatrixXd A_inv(A_Matrix.rows(), A_Matrix.cols());
  Eigen::FullPivLU<Eigen::MatrixXd> lu(A_Matrix);
  A_inv = lu.inverse();
  cout << "[Debug] inside polynomial trajectory, constructing A_invTrans" << endl; // for test
  MatrixXd A_invTrans = A_inv.transpose(); // selfadd:matrix inverse note!
  cout << "[Debug] inside polynomial trajectory, constructing R_Matrix" << endl; // for test
  MatrixXd R_Matrix = C_Matrix*A_invTrans*Q_Matrix*A_inv*C_Matrix_Trans;

  //dP∗=− RPP^−1 * RFP^T *dF
  cout << "[Debug] inside polynomial trajectory, constructing dF" << endl; // for test
  MatrixXd R_fp = R_Matrix.block(0,3*m+21,3*m+21,9*m-9);
  MatrixXd R_pp = R_Matrix.block(3*m+21,3*m+21,9*m-9,9*m-9);
  MatrixXd Rpp_inv = R_pp.inverse();
  MatrixXd Rfp_trans = R_fp.transpose();
  VectorXd dF = VectorXd::Zero(3*(m+7));
  // fill in dF 
  dF.block(0,0,3,1) = Path.row(0).transpose(); //Position0
  dF.block(3*(m+3),0,3,1) = Path.row(m-1).transpose(); //PositionM
  for(int j =1;j<m;j++){
    dF.block(12+3*(j-1),0,3,1) =Path.row(j).transpose(); //position condition(except p0 and pM)
  }
  cout << "[Debug] inside polynomial trajectory, constructing dP_star" << endl; // for test
  VectorXd dP_star = -1*Rpp_inv*Rfp_trans*dF;

  //P = A^-1*C^T*[dp,df]^T
  cout << "[Debug] inside polynomial trajectory, constructing dPF" << endl; // for test
  VectorXd dPF(12*m+12); 
  dPF << dF,dP_star;
  cout << "[Debug] inside polynomial trajectory, constructing PolyCoeff" << endl; // for test
  PolyCoeff = A_inv*C_Matrix_Trans*dPF;
  PolyCoeff.resize(m, 3 * p_num1d);
  return PolyCoeff;
}

double TrajectoryGeneratorWaypoint::getObjective() {
  _qp_cost = (_Px.transpose() * _Q * _Px + _Py.transpose() * _Q * _Py +
              _Pz.transpose() * _Q * _Pz)(0);
  return _qp_cost;
}

Vector3d TrajectoryGeneratorWaypoint::getPosPoly(MatrixXd polyCoeff, int k,
                                                 double t){
  Vector3d ret;
  int _poly_num1D = (int)polyCoeff.cols() / 3;
  for (int dim = 0; dim < 3; dim++) {
    VectorXd coeff = (polyCoeff.row(k)).segment(dim * _poly_num1D, _poly_num1D);
    VectorXd time = VectorXd::Zero(_poly_num1D);

    for (int j = 0; j < _poly_num1D; j++){
        if (j == 0){
          time(j) = 1.0;
        }
        else{
          time(j) = pow(t, j); //selfadd:[1, t, t^2, t^3, ...]
        }
      ret(dim) = coeff.dot(time); //selfadd: 使用点积计算多项式在时间 t 的值。这等价于计算 a0 + a1*t + a2*t^2 + ...
      // cout << "dim:" << dim << " coeff:" << coeff << endl;
    }
  }
  return ret;
}

Vector3d TrajectoryGeneratorWaypoint::getVelPoly(MatrixXd polyCoeff, int k,
                                                 double t) {
  Vector3d ret;
  int _poly_num1D = (int)polyCoeff.cols() / 3;
  for (int dim = 0; dim < 3; dim++) {
    VectorXd coeff = (polyCoeff.row(k)).segment(dim * _poly_num1D, _poly_num1D);
    VectorXd time = VectorXd::Zero(_poly_num1D);

    for (int j = 0; j < _poly_num1D; j++)
      if (j == 0)
        time(j) = 0.0;
      else
        time(j) = j * pow(t, j - 1);

    ret(dim) = coeff.dot(time);
  }

  return ret;
}

Vector3d TrajectoryGeneratorWaypoint::getAccPoly(MatrixXd polyCoeff, int k,
                                                 double t){
  Vector3d ret;
  int _poly_num1D = (int)polyCoeff.cols() / 3;
  for (int dim = 0; dim < 3; dim++) {
    VectorXd coeff = (polyCoeff.row(k)).segment(dim * _poly_num1D, _poly_num1D);
    VectorXd time = VectorXd::Zero(_poly_num1D);

    for (int j = 0; j < _poly_num1D; j++){
      if (j == 0 || j == 1){
        time(j) = 0.0;
      }
      else{
        time(j) = j * (j - 1) * pow(t, j - 2);
      }
    }
    ret(dim) = coeff.dot(time);
  }

  return ret;
}