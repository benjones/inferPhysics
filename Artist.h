#pragma once

#include <iostream>
#include <Eigen/Eigen>
#include <vector>

using Eigen::MatrixXd;

class Artist {
public:
	std::string filename;
	std::vector<double> time;
	std::vector<int> frame;
	MatrixXd X;

	void loadJsonFile(const std::string fileName);
	int degreesOFreedom;
	int hiddenDegrees;
	int timeSteps;
	int fps;
	int numFrames;

private:
	int collisionState = 1;
};