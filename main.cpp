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
	//guessI.topRows(s.degreesOfFreedom) =
		//s.snapshots.col(0).template cast<Scalar>();
	allTimeSteps.col(0).topRows(s.degreesOfFreedom) =
					  s.snapshots.col(0).template cast<Scalar>();
	//std::cout << allTimeSteps.col(0) << std::endl;
	
	for (int t = 0; t < timeSteps - 1; t++) {
		guessI = allTimeSteps.col(t);
		allTimeSteps.col(t+1) = M*guessI;
		//allTimeSteps.col(t + 1) = M*allTimeSteps.col(t); // This crashes the program, I feel like this is the same as the two lines above.
	}

	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> interpTimeStep = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(
		s.totalDOF, 1);
	for(int i = 0; i < s.snapshots.cols(); i++) {
		//Time of current snapshot being looked at
		int j = (s.frameTimes[i] / dt);
		double tj  = j * dt;
		Scalar alpha{(tj - s.frameTimes[i])/dt};
		Scalar oppAlpha{1.0 - alpha};
		//std::cout << "alpha: " << alpha << std::endl;
		//std::cout << "oppAlpha: " << oppAlpha << std::endl;
		/*
		Eigen::Matrix<Scalar, Eigen::Dynamic, 1> cj = allTimeSteps.col(j);
		Eigen::Matrix<Scalar, Eigen::Dynamic, 1> cjp1 = allTimeSteps.col(j + 1);
		
		interpTimeStep = cj*oppAlpha;
		interpTimeStep += cjp1*alpha;*/
		//std::cin.get();
		for(auto r = 0; r < interpTimeStep.rows(); ++r){
		  interpTimeStep(r) = (oppAlpha*allTimeSteps(r, j) + alpha*allTimeSteps(r, j+1));
		}
		
		ret += (interpTimeStep.topRows(s.degreesOfFreedom) - s.snapshots.col(i).topRows(s.degreesOfFreedom).template cast<Scalar>()).squaredNorm();
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
	Eigen::EigenSolver<MatrixXd> Meig(Mdoubles);
	// I'm not sure if using real values matters but it was the only way I was able to reconstruct Mdoubles with eigenvalues and eigenvectors.
	MatrixXd eigDiag = Meig.eigenvalues().real().asDiagonal();
	MatrixXd eigVectors = Meig.eigenvectors().real();

	// Print out eigValues of matrix and eigenVectors after converted to real values.
	//std::cout << "Diagonal Eigenvalues matrix: " << std::endl << eigDiag << std::endl;
	//std::cout << "Eigenvectors: " << std::endl << eigVectors << std::endl;
	
	// Clamping Eigenvalues to zero D + P, but not 100% sure what P should be, I guess i'm a touch confused if that should be the eigenvectors or if
	// it should be the bound limits we are trying to set. (lower) 
	// Or should I be doing something more like D + gamma*U for handling a lower bound?


	// However since I was struggling with the concept of what to add to D I tried multiplying the Diagonal Eigenvalue matrix together
	// it didn't provide 100% accurate results for solving the matrix, but it did seem like the snapshot values were starting to trend in the correct direction.
	std::cout << "eigDiag*eigDiag: " << std::endl << eigDiag*eigDiag << std::endl;


	// A = Q*D*Q^-1, what we want to do. Currently doing Q * (D*D) * Q^-1
	Mdoubles = (eigVectors*(eigDiag*eigDiag)*eigVectors.inverse());
	//Mdoubles = (Meig.eigenvectors()*Meig.eigenvalues().asDiagonal()*Meig.eigenvectors().inverse()).real();
	std::cout << "using real and imaginary: " << std::endl << Meig.eigenvectors()*Meig.eigenvalues().asDiagonal()*Meig.eigenvectors().inverse() << std::endl;
	//std::cout << "Meig.eigenvectors()*(Meig.eigenvalues()*Meig.eigenvalues())*Meig.eigenvectors().inverse(): " << std::endl << (eigVectors*(eigDiag*eigDiag)*eigVectors.inverse()) << std::endl;
	//std::cout << "Mdoubles Eigenvalues: " << std::endl << Mdoubles.eigenvalues() << std::endl;
	std::cout << "Mdoubles: " << std::endl << Mdoubles << std::endl;
	std::cin.get();
	if (dt > (1.0/s.fps)) {
		return computeMatrix(Mdoubles.sqrt(), s, dt / 2);
	}
	else {
		return Mdoubles;
	}
	
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
	//s.loadJsonFile("../Data/HitandThud.json");
	//s.loadJsonFile("../Data/EqualBounce.json");
	//s.loadJsonFile("../Data/ProjectileMotion.json");
	s.loadJsonFile("../Data/SpringForce.json");

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

	std::cout << realM << std::endl;
	std::cout << s.snapshots << std::endl;

	predictedPath(a, s);

	return 0;

}
