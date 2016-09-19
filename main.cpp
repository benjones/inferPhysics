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

//each column of targets is a snapshot we want to hit
Eigen::MatrixXd comptueEnergy(const MatrixXd& targets, const MatrixXd& M){

  
  
  
}

//
Eigen::MatrixXd computeEnergyGradient(const MatrixXd& targets, const MatrixXd& M){


  
}


int main(){

  std::cout << makeProjectileMotion(0, 50, 10) << std::endl;
  return 0;
  
	//Declarations
	float alpha = 0.1;
	float error = 0.0001;
	Eigen::Matrix2f m;
	m.setIdentity(2, 2);
	
	
	Eigen::Vector2f x0(1, 1);
	
	/*
	Eigen::Matrix2f m1;
	m1(0, 0) = 1;
	m1(0, 1) = -9.81;
	m1(1, 0) = 1;
	m1(1, 1) = 1;
	Eigen::ConjugateGradient<Eigen::Matrix2f> cg;
	cg.compute(m1);
	std::cout << cg.solve(x0) << std::endl;
	*/

	
	float f = 1;
	float gradF = 1;
	int i = 0;
	Eigen::Vector2f xi(1, i);
	Eigen::Matrix2f gradM;
	gradM(0, 0) = 1;
	gradM(0, 1) = -9.81;
	gradM(1, 0) = 1;
	gradM(1, 1) = 1;

	
	while (error < f || i < 100) {
		f = (xi - m*x0).transpose()*(xi - m*x0);
		gradF = ((-2 * x0.transpose())*(xi - m*x0))*((xi - m*x0).transpose()*(xi - m*x0));
		m = m - (alpha*gradF)*gradM;
		std::cout << f << std::endl;
		i++;
	}


  return 0;
}
