
#include "jsoncpp/include/json/json.h"
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

	numImages = root.get("snapshots", 0).asInt();
	degreesOFreedom = root.get("degreesFreedom", 0).asInt();
	hiddenDegrees = root.get("hiddenDegress", 0).asInt();
	timeSteps = root.get("nSteps", 0).asInt();
	
	X = MatrixXd::Zero(degreesOFreedom + hiddenDegrees, numImages);
	for (int i = 0; i < numImages; i++) {
		time.push_back(root["time"][i].asDouble());
		data.push_back(root["position"][i].asDouble());
		velocity.push_back(root["velocity"][i].asDouble());
		
		X(0, i) = data.at(i);
		X(1, i) = velocity.at(i);
		if (hiddenDegrees == 1) {
			X(2, i) = time.at(i);
		}
	}
	
}
