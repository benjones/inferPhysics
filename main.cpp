#include <iostream>
#include <Eigen/Eigen>
#include <cmath>
#include <fstream>
#include "FADBAD++/fadiff.h"
#include <vector>
#include "Artist.h"

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


MatrixXd makeProjectileMotion(double x0, double v0, double t0, int nSteps){

  MatrixXd ret = MatrixXd::Zero(3, nSteps);

  ret(0,0) = x0;
  ret(1,0) = v0;
  ret(2, 0) = t0;

  for(int i = 1; i < nSteps; i++){
	//position update x_n +1 = x_n + v_n * dt
	ret(0, i) = ret(0, i -1) + ret(1, i - 1); //assume dt = 1

	ret(1, i) = ret(1, i -1) - 9.81;
	ret(2, i) = i;
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
Scalar computeEnergy(const MatrixXd& targets, const Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M, Artist s){

  Scalar ret{0};
	
	// Should give us a scalar value
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> guessI = s.X.col(0).template cast<Scalar>();
	int j = 0;
	for (int i = 0; i < s.timeSteps; i++) {
	  //this approach is probably more efficient and (arguably) clearer
	  //sum_i  exp(-i)|| what we have now||^2
		if (i == s.time.at(j)) {
			ret += exp(-i)*(guessI.topRows(2) - s.X.col(j).topRows(2).template cast<Scalar>()).squaredNorm();
			guessI = M*guessI;
			j++;
		}

		
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
  std::vector<double> a(s.timeSteps+2);
  std::vector<double> b(s.timeSteps+2);
  a[0] = 0;
  b[0] = 0;
  int j = 1;
	for (int i = 1; i < s.timeSteps; i++) {
		if (i == s.time.at(j)) {
			b[i] = X(0, j);
			a[i] = M.row(0).dot(X.col(j-1));
			j++;
		}
		else {
			b[i] = std::numeric_limits<double>::quiet_NaN();
			a[i] = std::numeric_limits<double>::quiet_NaN();
		}
		
	}

	std::ofstream predictedStream("predictedPath", std::ios::binary);
	std::ofstream actualStream("actualPath", std::ios::binary);
	predictedStream.write(reinterpret_cast<const char*>(a.data()), sizeof(decltype(a)::value_type)*a.size());
	actualStream.write(reinterpret_cast<const char*>(b.data()), sizeof(decltype(b)::value_type)*b.size());

}


int main(){
	
  Artist s;
  s.loadJsonFile("ProjectileMotion.json");
  //s.loadJsonFile("SpringForce.json");

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


/* Left commented out for time being
//damping b = 2 * sqrt(k*m)
//makeSpringForce(x0,v0,mass,damping - b,Spring Constant - k,dt,nSteps)
double k = 1;
double m = 35;
double b = 1e-2;
double dt = 1;
int nSteps = 20;
X = makeSpringForce(0, 25, m, b, k, dt, nSteps);

X = makeProjectileMotion(0, 12.5, 0, nSteps);
}*/