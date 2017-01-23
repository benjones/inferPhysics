#pragma once

#include <iostream>
#include <Eigen/Eigen>
#include <vector>

using Eigen::MatrixXd;

class Artist {
public:
	std::string filename;
	std::vector<double> time;
	MatrixXd X;

	void loadJsonFile(const std::string fileName);
	int degreesOFreedom;
	int hiddenDegrees;
	int timeSteps;
private:
	int collisionState = 1;

	

};