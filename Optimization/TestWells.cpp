#include "TestWells.h"

void getOptData(PSOvalueType op) {
	std::cout << "gBestCost: " << op.second << std::endl;
	std::cout << "gBestPos: ";
	std::cout<< std::setprecision(7);
	for (auto x : op.first)
		std::cout << x << ", ";
	std::cout << std::endl;
}

void writeDataCartesian(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename) {
	std::ofstream output;
	output.open(filename);
	output << std::setprecision(9);
	for (size_t i = 0; i < pointsCartesian.size(); ++i) {
		output << pointsCartesian[i][0] << "," << pointsCartesian[i][1] << "," << pointsCartesian[i][2] << std::endl;
	}
	output.close();
}

void writeDataMD(std::vector<Eigen::Vector4d>& pointsMD, std::string filename) {
	std::ofstream output;
	output.open(filename);
	for (size_t i = 0; i < pointsMD.size(); ++i) {
		output << pointsMD[i][0] << "," << pointsMD[i][1] << "," << pointsMD[i][2] << "," << pointsMD[i][3] << std::endl;
	}
	output.close();
}

void writeInclinometry(const std::vector<Eigen::Vector4d>& pMD, std::string filename) {
	std::ofstream output;
	output.open(filename);
	output << "md, inc, azi\n";
	for (size_t i = 0; i < pMD.size(); ++i) {
		std::pair<double, double> incAzi = CartesianToSpherical(Eigen::Vector3d{ pMD[i][1],pMD[i][2],pMD[i][3] });
		output << pMD[i][0] << "," << incAzi.first << "," << incAzi.second << std::endl;
	}
	output.close();
}