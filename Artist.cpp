
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
	degreesOFreedom = root["degreesFreedom"].asInt();
	hiddenDegrees = root["hiddenDegrees"].asInt();
	collisionState = root.get("collisionState", 0).asInt();
	//use sensible default
	fps = root.get("fps", 10).asInt();
	//hidden degrees and collision stuff aren't part of the snapshot
	snapshots = MatrixXd::Zero(degreesOFreedom, root["SnapShot"].size());
	for (int i = 0; i < root["SnapShot"].size(); i++) {
	  frameNumbers.push_back(
		  static_cast<int>(std::round(root["Snapshot"][i]["time"].asDouble()*fps)));

		
		for(auto j = 0; j < degreesOfFreedom){
		  snapshots(j, i) = root["Snapshot"][i]["data"][j].asDouble();
		}

	}

	numFrames = frameNumbers.back();
	totalDOF = degreesOfFreedom + hiddenDegrees + collisionState;
}
