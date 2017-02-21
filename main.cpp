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
T square(const T& t) { return t*t; }


template <typename T>
std::ostream& operator <<(std::ostream& outs, FwdDiff<T> t) {
	outs << '[' << t.x() << " ";
	for (auto i = 0; i < t.size(); i++) {
		outs << t.d(i) << ' ';
	}
	outs << ']';
	return outs;
}


//const MatrixXd& targets, wasn't being used. Left commented out just incase I need to add it back in.
template <typename Scalar>
Scalar computeEnergy(const Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M, Artist s) {

	Scalar ret{ 0 };


	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> guessI =
		Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(
			s.totalDOF, 1);
	guessI.topRows(s.degreesOfFreedom) =
		s.snapshots.col(0).template cast<Scalar>();

	int j = 0;
	for (int i = 0; i <= s.numFrames; i++) {
		// Skips if i is not equal to next time step in the sequence	  
		if (i == s.frameNumbers[j]) {
			ret += exp(-((double)(s.frameNumbers[j] / s.numFrames)))*(guessI.topRows(s.degreesOfFreedom) -
				s.snapshots.col(j).topRows(
					s.degreesOfFreedom).template cast<Scalar>()).squaredNorm();
			j++;
		}
		guessI = M*guessI;
		//todo: turn this into a "handle collisions" function
		if (s.collisionState > 0 && guessI(0) < 0) {
			//handleCollisions(&guessI, s.totalDOF);  //once we have more complicaed collisions stuff, revisit this.
			guessI(s.degreesOfFreedom + s.hiddenDegrees) = -guessI(0);
		}

	}

	return ret;

}



// todo finish function to handle collisions
//Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M 
template <typename Scalar>
Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> handleEnergyCollisions(Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> XnplusOne, int DoF) {
	XnplusOne(DoF - 1) = -XnplusOne(0);

	return XnplusOne;
}

//Not being used.. Should I remove this  function?
Eigen::VectorXd handleCollisions(MatrixXd M, Eigen::VectorXd XnplusOne, int DoF) {
	XnplusOne(DoF - 1) = -XnplusOne(0);

	//add in -kx/m*dt for position and velocity
	//k - big? 
	//want me to add these to the json file?
	auto k = 3;
	auto m = 20;
	auto x = 0.5;
	auto dt = 0.1;

	//MatrixXd M;
	//M.setIdentity(DoF, DoF);
	M(0, 3) = -k / m * x * dt;
	M(1, 3) = -k / m * x * dt;
	XnplusOne = M*XnplusOne;

	return XnplusOne;
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


int main(int argc, char**argv) {

	Artist s;
	//s.loadJsonFile("../Data/Smallbounce.json");
	s.loadJsonFile("../Data/HitandThud.json");
	//s.loadJsonFile("../Data/EqualBounce.json");
	//s.loadJsonFile("../Data/ProjectileMotion.json");
	//s.loadJsonFile("../Data/SpringForce.json");

	/*if (argc < 2) {
		std::cout << "usage: inferphysics <json file>" << std::endl;
		exit(-1);
	}
	s.loadJsonFile(argv[1]);*/

	DiffMatrix M, Mprime;
	M.setIdentity(s.totalDOF, s.totalDOF);
	//M(1,2) = -9.81;

	//use good step size max of = 1, current alpa = max, do compute gradient update guess with doubles not auto diff
	//min alpha = 1e-12 if less break out
	const double minAlpha = 1e-12;
	const double maxAlpha = 1;
	double alpha = maxAlpha;
	double tol = 0.00001;
	int i = 0;
	int count = 0;

	double gradNorm;

	do {
		for (auto r = 0; r < M.rows(); r++) {
			for (auto c = 0; c < M.cols(); c++) {
				M(r, c).diff(r*M.cols() + c, M.size());
			}
		}

		auto energyAndDerivatives = computeEnergy(M, s);
		auto energy = energyAndDerivatives.val();
		gradNorm = 0;

		auto energyPrime = 0.0;
		Mprime = M;
		//Loop until energyPrime < energy
		do {
			for (auto r = 0; r < M.rows(); r++) {
				for (auto c = 0; c < M.cols(); c++) {
					Mprime(r, c) -= alpha*energyAndDerivatives.d(r*M.cols() + c);
					gradNorm += square(energyAndDerivatives.d(r*M.cols() + c));
				}
			}
			auto energyPrime = computeEnergy(Mprime, s).val();
			if (energy > energyPrime) {
				M = Mprime;
				break;
			}
			else {
				alpha /= 2.0;
			}

			if (energy == energyPrime) {
				count++;
			}
			/*if(count >= 3){
				M = Mprime;
				alpha *= 2.0;
				count = 0;
				break;
			}*/

		} while (energy <= energyPrime);

		
		/*for (auto r = 0; r < M.rows(); r++) {
			for (auto c = 0; c < M.cols(); c++) {
				M(r, c) -= alpha*(double)energyAndDerivatives.d(r*M.cols() + c);
				gradNorm += square(energyAndDerivatives.d(r*M.cols() + c));
			}
		}*/
		

		if (0 == i % 1000) {
			std::cout << energyAndDerivatives << " grad norm: " << gradNorm << std::endl;
		}

		i++;
	} while (gradNorm > tol && i < 20000);

	auto realM = convertToMatrixXd(M);

	std::vector<double> a(s.numFrames);
	Eigen::VectorXd currentState = Eigen::VectorXd::Zero(s.totalDOF);
	currentState.head(s.degreesOfFreedom) = s.snapshots.col(0);


	//currentState(0) -= 1;
	a[0] = currentState(0);
	std::cout << "0 : " << a[0] << std::endl;
	int j = 1;
	for (int i = 1; i < s.numFrames; i++) {
		currentState = realM*currentState;
		if (currentState(0) < 0 && s.collisionState == 1) {
			currentState(s.totalDOF - 1) = -currentState(0);
			//currentState = handleCollisions(realM, currentState, s.totalDOF);
		}
		a[i] = currentState(0);
		// Used to double check against snapshot values
		if (i == s.frameNumbers[j]) {
			std::cout << i << ": " << a[i] << std::endl;
			j++;
		}
	}


	std::cout << std::endl;
	std::cout << "Value of i at termination: " << i << std::endl;
	std::cout << "grad norm: " << gradNorm << std::endl;
	std::cout << convertToMatrixXd(M) << std::endl;
	std::cout << s.snapshots << std::endl;

	predictedPath(a, s);

	return 0;

}
