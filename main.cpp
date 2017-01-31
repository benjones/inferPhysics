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
T square(const T& t){ return t*t;}


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

  
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> guessI =
	  Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic>::Zero(
		  s.degreesOfFreedom + s.hiddenDegrees + s.collisionState, 1);
	guessI.head(s.degreesOfFreedom) =
	  s.snapshots.col(0).template cast<Scalar>();


	int j = 0;
	for (int i = 0; i < s.numFrames; i++) {
		// Skips if i is not equal to next time step in the sequence	  
		if (i == s.frameNumbers[j]) {
			ret += exp(-s.frameNumbers[j])*(guessI.topRows(s.degreesOfFreedom) -
				s.snapshots.col(j).topRows(
					s.degreesOfFreedom).template cast<Scalar>()).squaredNorm();
			j++;
		}
		//todo: turn this into a "handle collisions" function
		if (s.collisionState > 0 && guessI(0) < 0) {
			guessI(s.degreesOfFreedom + s.hiddenDegrees) = -guessI(0);
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


void predictedPath(MatrixXd M, MatrixXd path, Artist s) {
	std::vector<double> a(s.numFrames);
	std::vector<double> b(s.snapshots.cols());
	int k = 0;
	//a[0] = M.row(0)*s.X.col(0);
	for (int i = 0; i < s.numFrames; i++) {
		a[i] = path(0, i);
	}

	for (int j = 0; j < s.snapshots.cols(); j++) {
		b[j] = s.snapshots(0, j);
	}
	
	std::ofstream predictedStream("../Data/predictedPath", std::ios::binary);
	std::ofstream actualStream("../Data/actualPath", std::ios::binary);
	predictedStream.write(reinterpret_cast<const char*>(a.data()), sizeof(decltype(a)::value_type)*a.size());
	actualStream.write(reinterpret_cast<const char*>(b.data()), sizeof(decltype(b)::value_type)*b.size());

}


int main(){

  Artist s;
  //s.loadJsonFile("../Data/RandomSnapShots.json");
  //s.loadJsonFile("../Data/Random2.json");
  s.loadJsonFile("../Data/ProjectileMotion.json");
  //s.loadJsonFile("../Data/SpringForce.json");


  DiffMatrix M;
  // Adding 1 to account for collisions
  M.setIdentity(s.totalDOF, s.totalDOF);

  //replace with a guess function something more general
  /*if (s.hiddenDegrees == 1) {
	  M(1, 2) = -9.81;
	  }*/

  
  double alpha = 1e-5;	// Should allow for alpha value to change.
  double tol = 0.00001;
  int i = 0;

  double gradNorm;

  std::cout << "Number of frames: " << s.numFrames << std::endl;
  do {
	for(auto r = 0; r < M.rows(); r++){
	  for(auto c = 0; c < M.cols(); c++){
		M(r, c).diff(r*M.cols() + c, M.size());
	  }
	}
	
	
	auto energyAndDerivatives = computeEnergy(s.snapshots,M,s);
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

  auto realM = convertToMatrixXd(M);


  Eigen::VectorXd currentState = Eigen::VectorXd::Zero(s.totalDOF);
  currentState.head(s.degreesOfFreedom) = s.snapshots.col(0);
  for (int i = 0; i < s.numFrames; i++) {
	currentState = realM*currentState;
  }

 
  std::cout << std::endl;
  std::cout << "Value of i at termination: " << i << std::endl;
  std::cout << "grad norm: " << gradNorm << std::endl;
  std::cout << convertToMatrixXd(M) << std::endl;
  std::cout << s.snapshots << std::endl;

  predictedPath(convertToMatrixXd(M), next, s);

 

  return 0;

}

