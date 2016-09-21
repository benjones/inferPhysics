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

//todo add one for a spring force to test with
/*
* Do you mean add function in to add a spring force at a specific time step or did you want me to
* add a function that replicates spring force? Like we did with projectile motion above.
*/


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

  X = makeProjectileMotion(0, 50, 10);
 
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
