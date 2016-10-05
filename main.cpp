#include <iostream>
#include <Eigen/Eigen>
#include <cmath>
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

// Generate Spring Force Data
MatrixXd makeSpringForce(double x0, double v0, double m, double b, double k, double dt, int nSteps) {
	MatrixXd M = MatrixXd::Zero(2, 2);
	MatrixXd ret = MatrixXd::Zero(2, nSteps);

	M(0, 0) = 1 - (k / m)*(dt*dt);
	M(0, 1) = (1 - (b/m)*dt)*dt;
	M(1, 0) = (-k / m)*dt;
	M(1, 1) = (1 - (b/m)*dt);

	ret(0, 0) = x0;
	ret(1, 0) = v0;

	// k - spring constant, b - damping, m - mass
	for (int i = 1; i < nSteps; i++) {
		ret.col(i) = M*ret.col(i-1);
	}

	std::cout << "M: " << M << std::endl;

	std::cout << ret << std::endl;
	std::cin.get();

	return ret;
}



template <typename Scalar>
Scalar computeEnergy(const MatrixXd& targets, const Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M){

	Scalar ret = 0;
	
	// Should give us a scalar value
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> guessI = targets.col(0).template cast<Scalar>();

	for (int i = 0; i < targets.cols(); i++) {
	  //this approach is probably more efficient and (arguably) clearer
	  // sum_i  exp(-i)|| what we have now||^2
	  ret +=  exp(-i)*(guessI - targets.col(i).template cast<Scalar>()).squaredNorm();
	  
	  guessI = M*guessI;

	}

	return ret;
 
}



/*
* Create function that checks Mguess against M generated
*/
bool checkMguess(MatrixXd M, MatrixXd Mguess) {
	double tol = 0.00001;

	double normDiff = abs(Mguess.norm() - M.norm());
	if (tol > normDiff)
		return true;

	return false;
}



MatrixXd convertToMatrixXd(DiffMatrix M) {
	MatrixXd ret = MatrixXd::Zero(2, 2);
	ret(0, 0) = M(0, 0).val();
	ret(0, 1) = M(0, 1).val();
	ret(1, 0) = M(1, 0).val();
	ret(1, 1) = M(1, 1).val();

	return ret;
}


void predictedPath(MatrixXd M, MatrixXd X) {
	double *a = new double[10];
	a[0] = 0;
	for (int i = 1; i < X.cols(); i++) {
		a[i] = M.row(0).dot(X.col(i - 1));
	}

	FILE *f;
	f = fopen("data.txt", "wb");
	if (f == NULL) {
		std::cout << "Error" << std::endl;
	}
	for (int i = 0; i < 10; i++) {
		fwrite(&a[i], sizeof(a[i]), 1, f);
	}
	
	fclose(f);
	delete[] a;

}

int main(){
  
  MatrixXd X, gradF;
  DiffMatrix M;
  M.setIdentity(2, 2);

  //X = makeProjectileMotion(0, 50, 10);

  //damping b = 2 * sqrt(k*m)
  //makeSpringForce(x0,v0,mass,damping - b,Spring Constant - k,dt,nSteps)
  double k = 1;
  double m = 20;
  double b = 1e-2;
  //double b = 2 * sqrt(k*m);
  X = makeSpringForce(0, 25, m, b, k, 1, 10);
 
  
  double alpha = 1e-5;
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
  }while (gradNorm > tol && i < 500);
  
  std::cout << convertToMatrixXd(M) << std::endl;

  predictedPath(convertToMatrixXd(M), X);

  X = makeSpringForce(0, 25, m, b, k, 1, 10);
 

  return 0;

}
