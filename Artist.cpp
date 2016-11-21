
#include "json.h"
#include "Artist.h"
#include <fstream>


void Artist::loadJsonFile(const std::string filename) {
	std::ifstream input(filename);

	Json::Value root;
	Json::Reader reader;

	if (!reader.parse(input, root)) {
		std::cout << "Unable to read input file: " << filename << "\n"
			<< reader.getFormatedErrorMessages() << std::endl;

		exit(1);
	}

	degreesOFreedom = root.get("degreesFreedom", 0).asInt();
	hiddenDegrees = root.get("hiddenDegress", 0).asInt();
	timeSteps = root.get("nSteps", 0).asInt();
	

	X = MatrixXd::Zero(degreesOFreedom + hiddenDegrees, root["SnapShot"].size());
	for (int i = 0; i < root["SnapShot"].size(); i++) {
		time.push_back(root["SnapShot"][i]["time"].asDouble());

		X(0, i) = root["SnapShot"][i]["data"][0].asDouble();
		X(1, i) = root["SnapShot"][i]["data"][1].asDouble();
		if (hiddenDegrees == 1) {
			X(2, i) = time.at(i);
		}
	}
	
}
