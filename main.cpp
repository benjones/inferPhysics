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
		  s.totalDOF, 1);
	guessI.topRows(s.degreesOfFreedom) =
	  s.snapshots.col(0).template cast<Scalar>();	
	// using guessI.head(s.degreesOfFreedom) causes C2338 error, tried calling vector method on matrix. Changed head to
	// topRows(s.degreesOfFreedom)

	int j = 0;
	for (int i = 0; i <= s.numFrames; i++) {
		// Skips if i is not equal to next time step in the sequence	  
		if (i == s.frameNumbers[j]) {
			ret += (guessI.topRows(s.degreesOfFreedom) -
				s.snapshots.col(j).topRows(
					s.degreesOfFreedom).template cast<Scalar>()).squaredNorm();
			j++;
			//exp(-s.frameNumbers[j]), removed to allow for more than 1 iteration to occur.
		}
		//todo: turn this into a "handle collisions" function
		if (s.collisionState > 0 && guessI(0) < 0) {
			guessI(s.degreesOfFreedom + s.hiddenDegrees) = -guessI(0);
			guessI(0) = 0;
			//handleCollisions(&guessI, s.totalDOF);
		}
		guessI = M*guessI;
	}

	return ret;
 
}

// todo finish function to handle collisions
// Want to make a function that will take any type of vector / matrix and account for the collisions. 
template  <typename Scalar>
void handleCollisions(Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> c, int DoF) {
	c(DoF) = -c(0);
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


void predictedPath(std::vector<double> a, Artist s) {
	std::vector<double> b(s.snapshots.cols());

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
  //s.loadJsonFile("../Data/Random3.json");
  //s.loadJsonFile("../Data/ProjectileMotion.json");
  s.loadJsonFile("../Data/SpringForce.json");


  DiffMatrix M;
  M.setIdentity(s.totalDOF, s.totalDOF);

  //replace with a guess function something more general
  //need some clarification. Should I try and predict the type of value that would account for a hidde degree? 
  /*if (s.hiddenDegrees == 1) {
	  M(1, 2) = -9.81;
	  }*/

  
  double alpha = 1e-10; // Needs to be able to change...
  double tol = 0.00001;
  int i = 0;

  double gradNorm;

  //std::cout << "Number of frames: " << s.numFrames << std::endl;
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

  std::vector<double> a(s.numFrames+1);
  Eigen::VectorXd currentState = Eigen::VectorXd::Zero(s.totalDOF);
  currentState.head(s.degreesOfFreedom) = s.snapshots.col(0);
  int j = 0;
  for (int i = 0; i <= s.numFrames; i++) {
	currentState = realM*currentState;
	if (currentState(0) < 0) {
		currentState(s.totalDOF) = -currentState(0);
		currentState(0) = 0;
	}
	a[i] = currentState(0);
	// Used to double check value
	/*if (i == s.frameNumbers[j]) {
		std::cout << i << ": " << a[i] << std::endl;
		j++;
	}*/
  }
  //std::ofstream predictedStream("../Data/predictedPath", std::ios::binary);
  //predictedStream.write(reinterpret_cast<const char*>(a.data()), sizeof(decltype(a)::value_type)*a.size());

 
  std::cout << std::endl;
  std::cout << "Value of i at termination: " << i << std::endl;
  std::cout << "grad norm: " << gradNorm << std::endl;
  std::cout << convertToMatrixXd(M) << std::endl;
  std::cout << s.snapshots << std::endl;

  predictedPath(a, s);

  return 0;

}

