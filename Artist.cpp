
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

	//crash if not provided
	degreesOfFreedom = root["degreesFreedom"].asInt();
	hiddenDegrees = root["hiddenDegrees"].asInt();
	collisionState = root.get("collisionState", 0).asInt();
	//use sensible default
	fps = root.get("fps", 10).asInt();
	//std::cout << " Degrees of Freedom: " << degreesOfFreedom << std::endl;
	//hidden degrees and collision stuff aren't part of the snapshot
	snapshots = MatrixXd::Zero(degreesOfFreedom, root["Snapshot"].size());
	for (int i = 0; i < root["Snapshot"].size(); i++) {
	  frameTimes.push_back(root["Snapshot"][i]["time"].asDouble());
	  frameNumbers.push_back(static_cast<int>(std::round(frameTimes.back()*fps)));
	  for (auto j = 0; j < degreesOfFreedom; j++) {
		  snapshots(j, i) = root["Snapshot"][i]["data"][j].asDouble();
		}

	}

	numFrames = frameNumbers.back();
	totalDOF = degreesOfFreedom + hiddenDegrees + collisionState;
}
