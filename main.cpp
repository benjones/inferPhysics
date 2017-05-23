#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <fstream>
#include "FADBAD++/fadiff.h"
#include "FADBAD++/fadbad.h"
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

/*	if (t == s.frameTimes[k] | t == s.frameNumbers[k]) {
			guessI = s.snapshots.col(k).topRows(s.degreesOfFreedom).template cast<Scalar>();
			k++;
		}
		else {
			guessI = allTimeSteps.col(t);
		}
*/
//const MatrixXd& targets, wasn't being used. Left commented out just incase I need to add it back in.
template <typename Scalar>
Scalar computeEnergy(const Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M, Artist s, double dt) {
	Scalar ret{ 0 };
	if (dt > 1) {
		dt = 1;
	}
	//Final time in sequence
	int timeSteps = ceil(s.frameTimes.back()/dt) + 1;
	(assert(timeSteps > s.frameTimes.back() / dt));
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> allTimeSteps = 
		Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(s.totalDOF, timeSteps+1);
	//std::cout <<"Number of timeSteps: " << timeSteps << std::endl;

	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> guessI =
	  Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(
			s.totalDOF, 1);
	allTimeSteps.col(0).topRows(s.degreesOfFreedom) =
					  s.snapshots.col(0).template cast<Scalar>();
	for (int t = 0; t < timeSteps - 1; t++) {
		guessI = allTimeSteps.col(t);
		allTimeSteps.col(t+1) = (M*guessI);	// Was: M*(M*guessI). This caused the energy of a known matrix to explode.
	}

	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> interpTimeStep = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(
		s.totalDOF, 1);
	//Eigen::Matrix<Scalar, Eigen::Dynamic, 1> checkVals = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(
		//s.totalDOF, 1);
	for (int i = 0; i < s.snapshots.cols(); i++) {
		//Time of current snapshot being looked at
		int j = (s.frameTimes[i] / dt);
		double tj = j * dt;
		double alpha = (s.frameTimes[i] - tj) / dt;
		double oppAlpha = 1.0 - alpha;
		for (auto r = 0; r < s.degreesOfFreedom; ++r) {
			interpTimeStep(r) = (oppAlpha*allTimeSteps(r, j) + alpha*allTimeSteps(r, j + 1));
		}
		/*std::cout << "oppAlpha*allTimeSteps(0, j): " << oppAlpha*allTimeSteps(0, j) << " oppAlpha*allTimeSteps(0, j): " << oppAlpha*allTimeSteps(0, j) << std::endl;
		std::cout << "alpha*allTimeSteps(0,j+1): " << alpha*allTimeSteps(0, j + 1) << " alpha*allTimeSteps(1,j+1): " << alpha*allTimeSteps(1, j + 1) << std::endl;
		std::cout << "interpTimeStep(0): " << interpTimeStep(0) << " interpTimeStep(1): " << interpTimeStep(1)<< std::endl;
		std::cout << "Snapshot: " << s.snapshots.col(i) << std::endl;
		
		// Sanity Check....
		checkVals = (interpTimeStep - s.snapshots.col(i).topRows(s.degreesOfFreedom).template cast<Scalar>());

		std::cout << "Sanity Check (ref to checkVals vector): " << checkVals(0) << " " << checkVals(1) << std::endl;*/

		ret += (interpTimeStep - s.snapshots.col(i).topRows(s.degreesOfFreedom).template cast<Scalar>()).squaredNorm();

		//std::cout << ret << std::endl;
		//std::cin.get();
		

		// Used to debug negative alpha value.
		/*if (alpha < 0) {
			std::cout << "alpha: " << alpha << " oppAlpha: " << oppAlpha << std::endl;
			std::cout << "tj: " << tj << std::endl;
			std::cout << "frameTime[i]: " << s.frameTimes[i] << std::endl;
			std::cout << "j: " << j << " j+1: " << (j + 1) << std::endl;
			std::cout << "j Frametime: " << j*dt << " j+1 FrameTime: " << (j + 1)*dt << std::endl;
			std::cin.get();
		}*/
	}

	return ret;
}

/*
// todo finish function to handle collisions
template <typename Scalar>
Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> handleEnergyCollisions(Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> XnplusOne, int DoF) {
	XnplusOne(DoF - 1) = -XnplusOne(0);

	return XnplusOne;
}
*/

/*
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
*/


// Create function that checks Mguess against M generated
bool checkMguess(MatrixXd M, MatrixXd Mguess) {
	double tol = 0.00001;

	double normDiff = std::abs(Mguess.norm() - M.norm());
	if (tol > normDiff)
		return true;

	return false;
}

// Compute output with momentum. 
std::vector<double> computeWithMomentum(MatrixXd M, Artist s) {
	double alpha = 0.02;
	double beta = 0.99;

	std::vector<double> a(s.numFrames);
	Eigen::VectorXd currentState = Eigen::VectorXd::Zero(s.totalDOF);
	Eigen::VectorXd Z = Eigen::VectorXd::Zero(s.totalDOF);
	Z.setOnes();
	currentState.head(s.degreesOfFreedom) = s.snapshots.col(0);

	std::cout << "Printing with Momentum" << std::endl;
	//currentState(0) -= 1;
	a[0] = currentState(0);
	std::cout << "0 : " << a[0] << std::endl;
	int j = 1;
	for (int i = 1; i < s.numFrames; i++) {
		Z = beta*Z + M*currentState;
		currentState = currentState - alpha*Z;
		if (currentState(0) < 0 && s.collisionState == 1) {
			currentState(s.totalDOF - 1) = -currentState(0);
			//currentState = handleCollisions(realM, currentState, s.totalDOF);
		}
		a[i] = currentState(0);
		// Used to double check against snapshot values
		if (i == s.frameNumbers[j]) {
			std::cout << i << ": actual: " << a[i] << " expected: " << s.snapshots(0, j) << std::endl;
			j++;
		}
	}

	return a;
}

MatrixXd computeMatrix(MatrixXd Mcomp, Artist s, double dt) {
	MatrixXd Mprime;
	DiffMatrix M = Mcomp.template cast<FwdDiff<double>>();

	//use good step size max of = 1, current alpa = max, do compute gradient update guess with doubles not auto diff
	//min alpha = 1e-12 if less break out
	const double minAlpha = 1e-12;
	const double maxAlpha = 1;

	double alpha = maxAlpha;


	double tol = 1e-6;
	int i = 0;
	int count = 0;

	double gradNorm;
	std::cout << "dt: " << dt << std::endl;

	
	do {

		for (auto r = 0; r < M.rows(); r++) {
			for (auto c = 0; c < M.cols(); c++) {
				M(r, c).diff(r*M.cols() + c, M.size());
			}
		}

		auto energyAndDerivatives = computeEnergy(M, s, dt);
		auto energy = energyAndDerivatives.val();

		double energyPrime = 0.0;
		//Loop until energyPrime < energy
		do {
			Mprime = M.template cast<double>();
			// I believe the implementation will go as such M(r, c) -= alpha*(gamma*energyAndDerivatives.d(r*M.cols() + c)... gamma = (1-beta^(k+1-i))/(1-beta)
			for (auto r = 0; r < M.rows(); r++) {
				for (auto c = 0; c < M.cols(); c++) {
					Mprime(r, c) -= alpha*energyAndDerivatives.d(r*M.cols() + c);
				}
			}
			energyPrime = computeEnergy(Mprime, s, dt);
			if (energy > energyPrime) {

				//success, energy went down
				M = Mprime.template cast<FwdDiff<double>>();
				count++;
			}
			else {
				//failure, need to halve alpha
				count = 0;
				alpha /= 2.0;
			}

			if (count >= 5) {
				//multiple consecutive successes, yay
				alpha *= 2.0;
				count = 0;

			}
			if (alpha < minAlpha) {
				alpha = minAlpha;
				//std::cout << "alpha is too small" << std::endl;
				break;
			}
			else if (alpha > maxAlpha) {
				alpha = maxAlpha;

			}

		} while (energy <= energyPrime);

		gradNorm = 0;
		for (auto r = 0; r < M.rows(); r++) {
			for (auto c = 0; c < M.cols(); c++) {
				gradNorm += square(energyAndDerivatives.d(r*M.cols() + c));
			}
		}

		if (0 == i % 1000) {
			std::cout << energyAndDerivatives << " grad norm: " << gradNorm << std::endl;
			std::cout << "Alpha: " << alpha << std::endl;
		}

		i++;

	} while (gradNorm > tol && i < 10000);
	std::cout << "Value of i at termination: " << i << std::endl;
	std::cout << "grad norm: " << gradNorm << std::endl;
	// Print out Eigenvalues
	MatrixXd Mdoubles = M.template cast<double>();
	// Clamping Eigenvalues to zero and taking squareroot
	/*Eigen::EigenSolver<MatrixXd> Meig(Mdoubles);

	MatrixXd eigValsDiag = Meig.eigenvalues().real();
	MatrixXd eigVectors = Meig.eigenvectors().real();

	// Clamping eigenvalues to zero
	for (int r = 0; r < eigValsDiag.rows(); r++) {
		if (eigValsDiag(r) < 0) {
			eigValsDiag(r) = 0;
		}
	}
	
	//std::cout << "Diagonal Eigenvalues matrix: \n" << eigValsDiag << std::endl;
	
	// sqrt of eigenvalues - Diagonal Matrix
	eigValsDiag = eigValsDiag.array().sqrt();
	
	//std::cout << "Sqrt diagonal Eigenvalues matrix: \n" << eigValsDiag << std::endl;
	
	// recomputing A = Q*D*Q^T
	Mdoubles = eigVectors*eigValsDiag.asDiagonal()*eigVectors.transpose();

	//std::cout << "Mdoubles.eigenvalues: \n" << Mdoubles.eigenvalues() << std::endl;

	std::cout << "Mdoubles: \n" << Mdoubles << std::endl;*/
	//std::cin.get();
	Mdoubles = Mdoubles*Mdoubles;
	if (dt > (1.0/s.fps)) {
		return computeMatrix(Mdoubles, s, dt / 2);
	}
	else {
		return Mdoubles;
	}
	
}

void predictedPath(std::vector<double> withMomentum, std::vector<double> withOutMomentum, Artist s) {
	std::vector<double> b(s.snapshots.cols());

	for (int j = 0; j < s.snapshots.cols(); j++) {
		b[j] = s.snapshots(0, j);
	}

	std::ofstream predictedStream("../Data/predictedPath", std::ios::binary);
	std::ofstream actualStream("../Data/actualPath", std::ios::binary);
	std::ofstream momentumStream("../Data/pathWithMomentum", std::ios::binary);
	predictedStream.write(reinterpret_cast<const char*>(withOutMomentum.data()), sizeof(decltype(withOutMomentum)::value_type)*withOutMomentum.size());
	momentumStream.write(reinterpret_cast<const char*>(withMomentum.data()), sizeof(decltype(withMomentum)::value_type)*withMomentum.size());
	actualStream.write(reinterpret_cast<const char*>(b.data()), sizeof(decltype(b)::value_type)*b.size());

}

// add something small to to the matrix and make sure they are correct before going into gradient descent loop, compute energy.
int main(int argc, char**argv) {

	Artist s;
	//s.loadJsonFile("../Data/Smallbounce.json");
	//s.loadJsonFile("../Data/HitandThud.json");
	//s.loadJsonFile("../Data/EqualBounce.json");
	//s.loadJsonFile("../Data/ProjectileMotion.json");
	s.loadJsonFile("../Data/SpringForce.json");

	// Creates exact set of snapshots trying to be duplicated
	double m = 25; double b = 1e-2; double k = 1; double dt1 = 0.2; int timeSteps = 96;
	CreatePath p;
	p.createSpringForcePath(0, 25, m, b, k, dt1, timeSteps);
	p.writeJsonFile("Test.json", dt1);
	

	/*if (argc < 2) {
		std::cout << "usage: inferphysics <json file>" << std::endl;
		exit(-1);
	}
	s.loadJsonFile(argv[1]);*/

	MatrixXd M;
	M.setIdentity(s.totalDOF, s.totalDOF);

	double dt = (1.0 / s.fps);
	while (dt < 1) {
		dt = dt * 2;
	}

	Eigen::MatrixXd realM = computeMatrix(M, s, dt);
	std::vector<double> momentumV = computeWithMomentum(M, s);
	
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
		  std::cout << i << ": actual: " << a[i] << " expected: " << s.snapshots(0, j) << std::endl;
			j++;
		}
	}

	CreatePath p2;
	p2.writeCreatedPathJsonFile("TestOuput.json", 0.2, a.size(), a);
	p2.writeCreatedPathJsonFile("TestOuputWithMomentum.json", 0.2, momentumV.size(), momentumV);

	std::cout << realM << std::endl;
	std::cout << s.snapshots << std::endl;

	predictedPath(momentumV, a, s);

	return 0;

}
