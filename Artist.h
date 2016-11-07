#pragma once

#include <iostream>
#include <Eigen/Eigen>
#include <vector>

using Eigen::MatrixXd;

class Artist {
public:
	std::string filename;
	int numImages;
	std::vector<double> data;
	std::vector<double> time;
	std::vector<double> velocity;
	MatrixXd X;



	void loadJsonFile(const std::string fileName);
	int degreesOFreedom;
	int hiddenDegrees;
	int timeSteps;

	

};