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
  // cout << "[Debug] inside PolyQPGeneration, 111" << endl; // for test
  int m = Time.size();
  cout << "[Debug] m segments:" << m << endl; // for test
  MatrixXd PolyCoeff(m, 3 * p_num1d); //每行是一段polynomial，一行内，三个三个一组(x，y，z)，一共p_num1d组

  /**
   *
   * STEP 3.2:  generate a minimum-jerk piecewise monomial polynomial-based
   * trajectory
   *
   * **/
  //get C,A,Q
  // cout << "[Debug] inside PolyQPGeneration, 333" << endl; // for test
  Matrix3d I_3x3 = Matrix3d::Identity();
  MatrixXd I_12x12 = MatrixXd::Identity(12,12);
  MatrixXd One_3x3 = MatrixXd::Ones(3,3);
  double lambda = 1e-8;  // 正则化参数
  // cout << "[Debug] inside PolyQPGeneration, 444" << endl; // for test
  MatrixXd Q_Matrix = MatrixXd::Zero(3*p_num1d*m,3*p_num1d*m);// define size? //24mx24m
  MatrixXd A_Matrix = MatrixXd::Zero(3*p_num1d*m,3*p_num1d*m); //24mx24m
  MatrixXd C_Matrix_Trans = MatrixXd::Zero(24*m,12*m+12); //24mx12m+12
  // cout << "[Debug] inside PolyQPGeneration, 555" << endl; // for test
  C_Matrix_Trans.block(0,0,3*4,3*4) = I_12x12;
  C_Matrix_Trans.block(24*m-24+12,3*(m+3),3*4,3*4) = I_12x12;
  // cout << "[Debug] inside PolyQPGeneration, 666" << endl; // for test
  for(int i = 0;i<m;i++){
  double tm = Time(i);
  // cout << "[Debug] inside polynomial trajectory, constructing Qi" << endl; // for test
  //Construct Qi:
  MatrixXd Q_i = MatrixXd::Zero(3 * p_num1d, 3 * p_num1d);
  Q_i.block(21,21,3,3) = 705600*pow(tm,6)*One_3x3/pow(tm,6);
  Q_i.block(21,18,3,3) = 302400*pow(tm,5)*One_3x3/pow(tm,5);
  Q_i.block(21,15,3,3) = 100800*pow(tm,4)*One_3x3/pow(tm,4);
  Q_i.block(21,12,3,3) = 20160*pow(tm,3)*One_3x3/pow(tm,3);

  Q_i.block(18,21,3,3) = 302400*pow(tm,5)*One_3x3/pow(tm,5);
  Q_i.block(18,18,3,3) = 129600*pow(tm,4)*One_3x3/pow(tm,4);
  Q_i.block(18,15,3,3) = 43200*pow(tm,3)*One_3x3/pow(tm,3);
  Q_i.block(18,12,3,3) = 8640*pow(tm,2)*One_3x3/pow(tm,2);

  Q_i.block(15,21,3,3) = 100800*pow(tm,4)*One_3x3/pow(tm,4);
  Q_i.block(15,18,3,3) = 43200*pow(tm,3)*One_3x3/pow(tm,3);
  Q_i.block(15,15,3,3) = 14400*pow(tm,2)*One_3x3/pow(tm,2);
  Q_i.block(15,12,3,3) = 2880*pow(tm,1)*One_3x3/pow(tm,1);

  Q_i.block(12,21,3,3) = 20160*pow(tm,3)*One_3x3/pow(tm,3);
  Q_i.block(12,18,3,3) = 8640*pow(tm,2)*One_3x3/pow(tm,2);
  Q_i.block(12,15,3,3) = 2880*pow(tm,1)*One_3x3/pow(tm,1);
  Q_i.block(12,12,3,3) = 576*One_3x3;

  // 添加正则化
  // Q_i.diagonal().array() += lambda;

  // 或者对高阶项使用更大的正则化系数
  for(int i = 0; i < p_num1d; i++) {
      double reg = lambda * pow(10, i);  // 高阶项用更大的正则化系数
      Q_i.block(3*i,3*i,3,3).diagonal().array() += reg;
  }



  // cout<<"Q_i= "<<Q_i<<endl; // for test
  // Q_i <<  705600*pow(Time(i),6)*I_3x3, 302400*pow(Time(i),5)*I_3x3, 100800*pow(Time(i),4)*I_3x3, 20160*pow(Time(i),3)*I_3x3,
  //         302400*pow(Time(i),5)*I_3x3, 129600*pow(Time(i),4)*I_3x3, 43200*pow(Time(i),3)*I_3x3, 8640*pow(Time(i),2)*I_3x3,
  //         100800*pow(Time(i),4)*I_3x3, 43200*pow(Time(i),3)*I_3x3, 14400*pow(Time(i),2)*I_3x3, 2880*pow(Time(i),1)*I_3x3,
  //         20160*pow(Time(i),3)*I_3x3,  8640*pow(Time(i),2)*I_3x3,  2880*pow(Time(i),1)*I_3x3,  576;
  // cout << "[Debug] inside polynomial trajectory, constructing Q_Matrix" << endl; // for test
  Q_Matrix.block(3*p_num1d*i, 3*p_num1d*i, 3*p_num1d, 3*p_num1d) = Q_i;
  // cout << "[Debug] inside polynomial trajectory, constructing Ai" << endl; // for test
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
  // cout << "[Debug] inside polynomial trajectory, constructing A_Matrix" << endl; // for test
  A_Matrix.block(3*p_num1d*i, 3*p_num1d*i, 3*p_num1d, 3*p_num1d) = A_i; 
  // cout<< "A_i = " << A_i <<endl; // for test 
  
  // cout << "[Debug] inside polynomial trajectory, constructing C_Matrix"<<"m:"<<m<<"i:"<<i << endl; // for test
  //Construct C_Matrix: //24mx12m+12
  if(i>0){
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
  // cout << "[Debug] inside polynomial trajectory, constructing C_Trans" << endl; // for test
  MatrixXd C_Matrix = C_Matrix_Trans.transpose();
  // cout << "[Debug] inside polynomial trajectory, constructing A_inv" << endl; // for test
  MatrixXd A_inv(A_Matrix.rows(), A_Matrix.cols());
  Eigen::FullPivLU<Eigen::MatrixXd> lu(A_Matrix);
  A_inv = lu.inverse();
  cout<<"A_inv size: "<<A_inv.size()<< endl;
  // cout<<"A_inv = "<<A_inv<<endl; // for test
  // cout << "[Debug] inside polynomial trajectory, constructing A_invTrans" << endl; // for test
  MatrixXd A_invTrans = A_inv.transpose(); // selfadd:matrix inverse note!
  // cout << "[Debug] inside polynomial trajectory, constructing R_Matrix" << endl; // for test
  MatrixXd R_Matrix = C_Matrix*A_invTrans*Q_Matrix*A_inv*C_Matrix_Trans;
  // cout<<"R_Matrix: "<<R_Matrix<<endl; // for test

  //dP∗=− RPP^−1 * RFP^T *dF
  // cout << "[Debug] inside polynomial trajectory, constructing dF" << endl; // for test
  MatrixXd R_fp = R_Matrix.block(0,3*m+21,3*m+21,9*m-9);
  MatrixXd R_pp = R_Matrix.block(3*m+21,3*m+21,9*m-9,9*m-9);
  // cout<<"R_pp: "<<R_pp<<endl; // for test

  MatrixXd Rpp_inv(R_pp.rows(), R_pp.cols());
  Eigen::FullPivLU<Eigen::MatrixXd> lu_2(R_pp);
  Rpp_inv = lu_2.inverse();
  // MatrixXd Rpp_inv = R_pp.inverse();

  MatrixXd Rfp_trans = R_fp.transpose();
  // cout<<"Rpp_inv: "<<Rpp_inv<<endl; // for test
  // cout<<"Rfp_trans: "<<Rfp_trans<<endl; // for test

  VectorXd dF = VectorXd::Zero(3*(m+7));
  // fill in dF 
  dF.block(0,0,3,1) = Path.row(0).transpose(); //Position0
  dF.block(3*(m+3),0,3,1) = Path.row(m).transpose(); //PositionM
  for(int j =1;j<m;j++){
    dF.block(12+3*(j-1),0,3,1) =Path.row(j).transpose(); //position condition(except p0 and pM)
  }
  cout<<"dF: "<<dF; // for test
  // cout << "[Debug] inside polynomial trajectory, constructing dP_star" << endl; // for test
  VectorXd dP_star = -1*Rpp_inv*Rfp_trans*dF;
  cout<<"dP_star:"<<dP_star<<endl; // for test
  //P = A^-1*C^T*[dp,df]^T
  // cout << "[Debug] inside polynomial trajectory, constructing dPF" << endl; // for test
  VectorXd dPF(12*m+12); 
  dPF << dF,dP_star;
  // cout<<"dPF:"<<dPF<<endl; // for test

  VectorXd D;
  D = C_Matrix_Trans*dPF;
  cout<<"D_Matrix = "<<D<<endl;

  // cout << "[Debug] inside polynomial trajectory, constructing PolyCoeff" << endl; // for test
  PolyCoeff = A_inv*C_Matrix_Trans*dPF;
  // cout<<"Polycoeff  before resize:"<<PolyCoeff<<endl; // for test
  // PolyCoeff.resize(m, 3 * p_num1d);
  // MatrixXd PolyCoeff_resize = MatrixXd::Zero(3, 3);
  // Map<MatrixXd, RowMajor> row_major(PolyCoeff.data(), PolyCoeff.rows(), PolyCoeff.cols());
  // PolyCoeff.conservativeResize(m, 3 * p_num1d);

  //check continuous
  // VectorXd D_2;
  // D_2 = A_Matrix*PolyCoeff;
  // cout<<"D_2-D = "<<D_2-D<<endl;

  // Create new matrix with desired size
  MatrixXd PolyCoeff_resized(m, 3 * p_num1d);

  // Copy data in row-major order
  // for(int s = 0; s <  m; ++s) {
  //   for(int k = 0; k < 3 * p_num1d; ++k) { // modified: std::min(PolyCoeff.cols(), 3 * p_num1d)
  //     try{
  //       PolyCoeff_resized(s,k) = PolyCoeff(s*3 * p_num1d+k);
  //     }catch (const std::exception& e) {
  //     cout<<"resize failed.. i,j= "<< s<<" "<<k<<endl;
  //     }
  //   } 
  // }
  for(int s = 0; s < m; s++) {
    for(int dim = 0; dim < 3; dim++) {
        for(int j = 0; j < p_num1d; j++) {
          // cout<<"s,j,dim= "<<s<<j<<dim<<endl;
            PolyCoeff_resized(s, dim * p_num1d + j) = PolyCoeff(3 * s * p_num1d + j * 3 + dim);
        }
    }
  }
  cout<<"Polycoeff:"<<PolyCoeff_resized<<endl; // for test

  // 在得到PolyCoeff_resized后添加 for test
  for(int s = 0; s < m; s++) {
      std::cout << "段 " << s << " 的系数：\n";
      for(int dim = 0; dim < 3; dim++) {
          std::cout << "维度 " << dim << ": ";
          for(int j = 0; j < p_num1d; j++) {
              std::cout << PolyCoeff_resized(s, dim * p_num1d + j) << " ";
          }
          std::cout << "\n";
      }
  }

 // testing continuousty
  for(int i = 0; i < m-1; i++) {
    // 检查段i结束点和段i+1起始点的连续性
    Vector3d pos_end = getPosPoly(PolyCoeff_resized, i, Time(i));
    Vector3d pos_start = getPosPoly(PolyCoeff_resized, i+1, 0);
    Vector3d vel_end = getVelPoly(PolyCoeff_resized, i, Time(i));
    Vector3d vel_start = getVelPoly(PolyCoeff_resized, i+1, 0);
    Vector3d acc_end = getAccPoly(PolyCoeff_resized, i, Time(i));
    Vector3d acc_start = getAccPoly(PolyCoeff_resized, i+1, 0);
    
    std::cout << "段 " << i << " 和 " << i+1 << " 的连接点：\n";
    std::cout << "位置差： " << (pos_end - pos_start).norm() << "\n";
    std::cout << "速度差： " << (vel_end - vel_start).norm() << "\n";
    std::cout << "加速度差： " << (acc_end - acc_start).norm() << "\n";
}

// 检查A_Matrix的结构
std::cout << "A_Matrix 的条件数：" << 
    A_Matrix.jacobiSvd().singularValues()(0) / 
    A_Matrix.jacobiSvd().singularValues()(A_Matrix.jacobiSvd().singularValues().size()-1) 
    << "\n";

// 检查Q_Matrix的结构
std::cout << "Q_Matrix 非零元素的位置：\n";
for(int i = 0; i < Q_Matrix.rows(); i++) {
    for(int j = 0; j < Q_Matrix.cols(); j++) {
        if(abs(Q_Matrix(i,j)) > 1e-10) {
            std::cout << "(" << i << "," << j << "): " << Q_Matrix(i,j) << "\n";
        }
    }
}

// 检查R_Matrix的结构
std::cout << "R_fp的维度: " << R_fp.rows() << "x" << R_fp.cols() << "\n";
std::cout << "R_pp的维度: " << R_pp.rows() << "x" << R_pp.cols() << "\n";
std::cout << "R_pp的条件数: " << 
    R_pp.jacobiSvd().singularValues()(0) / 
    R_pp.jacobiSvd().singularValues()(R_pp.jacobiSvd().singularValues().size()-1) 
    << "\n";


  return PolyCoeff_resized;
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
    }
    ret(dim) = coeff.dot(time); //selfadd: 使用点积计算多项式在时间 t 的值。这等价于计算 a0 + a1*t + a2*t^2 + ...
      // cout << "dim:" << dim << " coeff:" << coeff << endl;
    
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