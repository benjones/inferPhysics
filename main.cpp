#include <iostream>
#include <Eigen/Eigen>
#include <cmath>
#include <fstream>
#include "FADBAD++/fadiff.h"
#include <vector>
#include "Artist.h"
#include "CreatePath.h"

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



template <typename Scalar>
Scalar computeEnergy(const MatrixXd& targets, const Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M, Artist s){

  Scalar ret{0};
	
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> guessI = s.X.col(0).template cast<Scalar>();
	int j = 0;
	for (int i = 0; i < s.timeSteps; i++) {
		// Skips if i is not equal to next time step in the sequence	  
		if (i == s.time.at(j)) {
			ret += exp(-i)*(guessI.topRows(2) - s.X.col(j).topRows(2).template cast<Scalar>()).squaredNorm();
			j++;
		}
		guessI = M*guessI;
	}

	return ret;
 
}



/*
* Create function that checks Mguess against M generated
*/
bool checkMguess(MatrixXd M, MatrixXd Mguess) {
	double tol = 0.00001;

	double normDiff = std::abs(Mguess.norm() - M.norm());
	if (tol > normDiff)
		return true;

	return false;
}



MatrixXd convertToMatrixXd(DiffMatrix M) {
	MatrixXd ret = MatrixXd::Zero(M.rows(), M.cols());
	for (int i = 0; i < M.rows(); i++) {
		for (int j = 0; j < M.cols(); j++) {
			ret(i, j) = M(i, j).val();
		}
  }
	return ret; //M.unaryExpr([](const FwdDiff<double>& x)-> double { return x.val(); }).eval();
}


void predictedPath(MatrixXd M, MatrixXd X, Artist s) {
  std::vector<double> a(s.timeSteps);
  std::vector<double> b(s.timeSteps);
  a[0] = 0;
  b[0] = 0;
  int j = 1;
  // Fix bugs, associated with predicted and actual path, rip out synthetic data gen.
	for (int i = 1; i < s.timeSteps; i++) {
		if (i == s.time.at(j)) {
			b[i] = X(0, j);
			a[i] = M(0,0)*s.X(0,j);
			j++;
		}
		else {
			b[i] = std::numeric_limits<double>::quiet_NaN();
			a[i] = std::numeric_limits<double>::quiet_NaN();
		}
		
	}

	std::ofstream predictedStream("../Data/predictedPath", std::ios::binary);
	std::ofstream actualStream("../Data/actualPath", std::ios::binary);
	predictedStream.write(reinterpret_cast<const char*>(a.data()), sizeof(decltype(a)::value_type)*a.size());
	actualStream.write(reinterpret_cast<const char*>(b.data()), sizeof(decltype(b)::value_type)*b.size());

}


int main(){
	
  Artist s;
  s.loadJsonFile("../Data/RandomSnapShots.json");
  //s.loadJsonFile("../Data/ProjectileMotion.json");
  //s.loadJsonFile("../Data/SpringForce.json");

  MatrixXd gradF;
  DiffMatrix M;
  M.setIdentity(s.degreesOFreedom+s.hiddenDegrees, s.degreesOFreedom + s.hiddenDegrees);
  if (s.hiddenDegrees == 1) {
	  M(1, 2) = -9.81;
  }
  
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
	
	
	auto energyAndDerivatives = computeEnergy(s.X,M,s);
	gradNorm = 0;
	for(auto r = 0; r < M.rows(); r++){
	  for(auto c = 0; c < M.cols(); c++){
		M(r, c) -= alpha*energyAndDerivatives.d(r*M.cols() + c);
		gradNorm += square(energyAndDerivatives.d(r*M.cols() + c));
	  }
	}

	if (0 == i % 1000) {
		std::cout << energyAndDerivatives << " grad norm: " << gradNorm << std::endl;
	}
	
	i++;
  }while (gradNorm > tol && i < 20000);
 
  std::cout << "Value of i at termination: " << i << std::endl;
  std::cout << "grad norm: " << gradNorm << std::endl;
  std::cout << convertToMatrixXd(M) << std::endl;

  predictedPath(convertToMatrixXd(M), s.X, s);

 

  return 0;

}

