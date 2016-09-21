#include <iostream>
#include <Eigen/Eigen>

using Eigen::MatrixXd;

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
Eigen::MatrixXd computeEnergy(const MatrixXd& targets, const MatrixXd& M){

	MatrixXd ret = MatrixXd::Zero(1,1);

	// Should give us a scalar value
	for (int i = 0; i < targets.cols(); i++) {
		ret += targets.col(i).transpose()*targets.col(i) - targets.col(i).transpose()*M.transpose()*targets.col(0) -
			targets.col(i).transpose()*M*targets.col(i) + targets.col(0).transpose()*M.transpose()*M*targets.col(0);
	}

	return ret;
 
}

// Mi is being used as Mi-1, I am assuming that I should be using the Matrix M from the previous iteration 
// of my while loop, so I added an argument to account for that.
Eigen::MatrixXd computeEnergyGradient(const MatrixXd& targets, const MatrixXd& M, const MatrixXd& Mi){

	MatrixXd ret = MatrixXd::Zero(2,2);

	
	// Should get a 2x2 Matrix, gradient has not been updated to use Mi-1 for first and last M. 
	// todo update first and last M to Mi-1 after function works.
	for (int i = 0; i < targets.cols(); i++) {
		ret += 2 * (M*targets.col(0) - targets.col(i))*(i*Mi*targets.col(0)).transpose();
		ret(0, 0) += -2 * i * (targets.col(i).transpose()*Mi)*targets.col(0);
		ret(0, 1) += -2 * i * (targets.col(i).transpose()*Mi)*targets.col(0);
		ret(1, 0) += -2 * i * (targets.col(i).transpose()*Mi)*targets.col(0);
		ret(1, 1) += -2 * i * (targets.col(i).transpose()*Mi)*targets.col(0);
	}

	return ret;

  
}


int main(){
  
  MatrixXd M, X, f, gradF, Mi;
  M.setIdentity(2, 2);
  Mi.setIdentity(2,2);
  X = makeProjectileMotion(0, 50, 10);

 
  double alpha = 0.1;
  double error = 0.00001;
	int i = 0;
	
	while (error < f.norm() || i < 10) {
		f = computeEnergy(X,M);
		gradF = computeEnergyGradient(X,M, Mi);
		Mi = M;
		M = M - alpha*gradF;
		std::cout << f << std::endl;
		i++;
	}


  return 0;

}
