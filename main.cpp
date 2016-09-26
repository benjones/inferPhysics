#include <iostream>
#include <Eigen/Eigen>

#include "FADBAD++/fadiff.h"

using Eigen::MatrixXd;
using Eigen::Matrix;

using fadbad::FwdDiff;
using DiffMatrix = Eigen::Matrix<FwdDiff<double>, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
T square(T t){ return t*t;}

template <typename T>
std::ostream& operator <<( std::ostream& outs, FwdDiff<T> t){
  outs << '[' << t.x() << " ";
  for(auto i = 0; i < t.size(); i++){
	outs << t.d(i) << ' ';
  }
  outs << ']';
  return outs;
}


MatrixXd makeProjectileMotion(double x0, double v0, int nSteps){

  MatrixXd ret = MatrixXd::Zero(2, nSteps);

  ret(0,0) = x0;
  ret(1,0) = v0;

  for(int i = 1; i < nSteps; i++){
	//position update x_n +1 = x_n + v_n * dt
	ret(0, i) = ret(0, i -1) + ret(1, i -1); //assume dt = 1

	ret(1, i) = ret(1, i -1) - 9.81;
  }

  return ret;

}

/*
* ret(0, i) = -ret(0, i - 1) - (m / k)*(-9.81);
* ret(1, i) = -(k / (2 * m))*(ret(0, i - 1)*ret(0, i - 1)) - 9.81;
* I have a few questions and my work with me, I'm not sure if I setup the differential equation quite right for spring force.
* I had: X" = -k/m * x - b/m *x' and v' = -k/m * x - b/m * x'
*/
MatrixXd makeSpringForce(double x0, double v0, double m, double b, double k, int nSteps) {
	MatrixXd ret = MatrixXd::Zero(2, nSteps);

	ret(0, 0) = x0;
	ret(1, 0) = v0;

	// k - spring constant, b - damping, m - mass
	// Initially will assume damping will be = 0
	for (int i = 1; i < nSteps; i++) {
		ret(0, i) = ret(0,i-1) +  (k / (2*m)) * ret(1, i - 1);
		ret(1, i) = ret(1,i-1) - 9.81;
	}
	return ret;
}


//each column of targets is a snapshot we want to hit



template <typename Scalar>
Scalar computeEnergy(const MatrixXd& targets, const Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M){

	Scalar ret = 0;
	
	// Should give us a scalar value
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> guessI = targets.col(0).template cast<Scalar>();

	for (int i = 0; i < targets.cols(); i++) {
	  //this approach is probably more efficient and (arguably) clearer

	  ret +=  (guessI - targets.col(i).template cast<Scalar>()).squaredNorm();
	  
	  guessI = M*guessI;

	}

	return ret;
 
}


int main(){
  
  MatrixXd X, gradF;
  DiffMatrix M;
  M.setIdentity(2, 2);

  //X = makeProjectileMotion(0, 50, 10);

  X = makeSpringForce(0, 50, 20, 0, 1, 10);
  //std::cout << makeSpringForce(0, 50, 20, 0, 1, 10) << std::endl;
 
  double alpha = 1e-7;
  double tol = 0.00001;
  int i = 0;
  double gradNorm;
  do {
	for(auto r = 0; r < M.rows(); r++){
	  for(auto c = 0; c < M.cols(); c++){
		M(r, c).diff(r*M.cols() + c, M.size());
	  }
	}

	auto energyAndDerivatives = computeEnergy(X,M);
	gradNorm = 0;
	for(auto r = 0; r < M.rows(); r++){
	  for(auto c = 0; c < M.cols(); c++){
		M(r, c) -= alpha*energyAndDerivatives.d(r*M.cols() + c);
		gradNorm += square(energyAndDerivatives.d(r*M.cols() + c));
	  }
	}
	std::cout << energyAndDerivatives << " grad norm: " << gradNorm << std::endl ;
   
	i++;
  }while (gradNorm > tol && i < 100);
  

  return 0;

}
