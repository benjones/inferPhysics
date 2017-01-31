#pragma once

#include <iostream>
#include <Eigen/Eigen>
#include <vector>

using Eigen::MatrixXd;

class Artist {
public:
	std::string filename;
  std::vector<int> frameNumbers;
	MatrixXd snapshots;

	void loadJsonFile(const std::string fileName);
	int degreesOfFreedom;
	int hiddenDegrees;
	int fps;
	int numFrames;
  int collisionState;
  int totalDOF;
};
