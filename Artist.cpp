
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
	hiddenDegrees = root.get("hiddenDegrees", 0).asInt();
	timeSteps = root.get("nSteps", 0).asInt();
	fps = root.get("fps", 0).asInt();
	numFrames = fps*timeSteps;
	X = MatrixXd::Zero(degreesOFreedom + hiddenDegrees + collisionState, root["SnapShot"].size());
	for (int i = 0; i < root["SnapShot"].size(); i++) {
		time.push_back(root["SnapShot"][i]["time"].asDouble());
		frame.push_back(root["SnapShot"][i]["frame"].asInt());

		X(0, i) = root["SnapShot"][i]["data"][0].asDouble();
		X(1, i) = root["SnapShot"][i]["data"][1].asDouble();
		if (hiddenDegrees == 1) {
			X(2, i) = time.at(i);
			X(3, i) = 0;
		}
		else {
			X(2, i) = 0;
		}
	}
	
}
