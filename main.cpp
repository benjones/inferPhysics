#include <iostream>
#include <Eigen/Eigen>


int main(){

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
