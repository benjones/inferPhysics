#pragma once

#include <iostream>
#include <Eigen/Eigen>
#include <vector>

using Eigen::MatrixXd;

/*
NOT FINISHED OR TESTED YET, ADDING IN CLASS TO CREATE A PATH, CAN BE COMPARED AGAINST PATH CREATED FROM SNAPSHOTS
*/

class CreatePath {
public:
	MatrixXd X;

	void createSpringForcePath(double x0, double v0, double m, double b, double k, double dt, int nSteps);
	void createProjectileMotionPath(double x0, double v0, double t0, int nSteps);
	void writeJsonFile(const std::string fileName, double dt);
	void writeCreatedPathJsonFile(const std::string fileName, double dt, int numFrames, std::vector<double> s);
	void silly(double x0, double v0, double t0, int nSteps);
	int degreesOFreedom;
	int hiddenDegrees;
	int timeSteps;
	double mass, damping, springConstant, dt;



};
