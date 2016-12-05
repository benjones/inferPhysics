#include "CreatePath.h"
#include <Eigen/Eigen>
#include "json.h"
#include <fstream>

using Eigen::MatrixXd;


void CreatePath::createProjectileMotionPath(double x0, double v0, double t0, int nSteps) {

	X = MatrixXd::Zero(3, nSteps);

	timeSteps = nSteps;

	X(0, 0) = x0;
	X(1, 0) = v0;
	X(2, 0) = t0;

	for (int i = 1; i < nSteps; i++) {
		//position update x_n +1 = x_n + v_n * dt
		X(0, i) = X(0, i - 1) + X(1, i - 1); //assume dt = 1

		X(1, i) = X(1, i - 1) - 9.81;
		X(2, i) = i;
	}

}


void CreatePath::createSpringForcePath(double x0, double v0, double m, double b, double k, double dt, int nSteps) {

	MatrixXd M = MatrixXd::Zero(2, 2);
	X = MatrixXd::Zero(2, nSteps);

	timeSteps = nSteps;

	M(0, 0) = 1 - (k / m)*(dt*dt);
	M(0, 1) = (1 - (b / m)*dt)*dt;
	M(1, 0) = (-k / m)*dt;
	M(1, 1) = (1 - (b / m)*dt);

	X(0, 0) = x0;
	X(1, 0) = v0;

	// k - spring constant, b - damping, m - mass
	for (int i = 1; i < nSteps; i++) {
		X.col(i) = M*X.col(i - 1);
	}


}


void CreatePath::writeJsonFile(const std::string fileName) {
	std::ofstream file;
    file.open("../Data/" + fileName);

    Json::Value obj;
	double s;
	Json::Value v(Json::arrayValue);
	std::cout << "Enter Degrees of Freedom: " << std::endl;
	std::cin >> degreesOFreedom;
	obj["degreesFreedom"] = degreesOFreedom;
	std::cout << "Enter Hidden Degrees of Freedom: " << std::endl;
	std::cin >> hiddenDegrees;
	obj["hiddenDegrees"] = hiddenDegrees;
	obj["nSteps"] = timeSteps;
	
	
	for (int i = 0; i < timeSteps; i++) {
		obj["SnapShot"][i]["time"] = i;
		for (int j = 0; j < degreesOFreedom; j++) {
			obj["SnapShot"][i]["data"][j] = X(j, i);
		}
		
	}

	

    Json::StyledWriter styledWriter;
    file << styledWriter.write(obj);


    file.close();
}



/* Original functions used to construct and test
//damping b = 2 * sqrt(k*m)
//makeSpringForce(x0,v0,mass,damping - b,Spring Constant - k,dt,nSteps)
double k = 1;
double m = 35;
double b = 1e-2;
double dt = 1;
int nSteps = 20;
X = makeSpringForce(0, 25, m, b, k, dt, nSteps);

X = makeProjectileMotion(0, 12.5, 0, nSteps);
*/


/* In here for referene as I create a class that creates predicted paths and then outputs them to a JSON file.
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

*/