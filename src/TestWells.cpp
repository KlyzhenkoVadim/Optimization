#include "TestWells.h"

void getOptData(PsoValueType op) {
	std::cout << "gBestCost: " << op.cost << std::endl;
	std::cout << "gBestPos: ";
	std::cout<< std::setprecision(7);
	for (auto x : op.optVec)
		std::cout << x << ", ";
	std::cout << std::endl;
}

void writeDataCartesian(std::vector<Eigen::Vector3d>& pointsCartesian,
                        std::string filename) {
	std::ofstream output;
	output.open(filename);
	output << std::setprecision(9);
	for (size_t i = 0; i < pointsCartesian.size(); ++i) {
		output << pointsCartesian[i][0] << "," << pointsCartesian[i][1] << ","
				<< pointsCartesian[i][2] << std::endl;
	}
	output.close();
}

void writeInclinometry(const std::vector<Eigen::Vector4d>& pMD,
                       std::string filename) {
	std::ofstream output;
	output.open(filename);
	output << "md, inc, azi\n";
	for (size_t i = 0; i < pMD.size(); ++i) {
		std::pair<double, double> incAzi = cartesianToSpherical(
			Eigen::Vector3d{pMD[i][1], pMD[i][2], pMD[i][3]});
		output << pMD[i][0] << "," << incAzi.first << ","
				<< incAzi.second << std::endl;
	}
	output.close();
}
