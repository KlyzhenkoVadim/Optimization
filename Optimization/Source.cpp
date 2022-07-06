#include "Eigen/Dense"
#include <iostream>
#include <vector>
#include "TrajectoryTemplate.h"
#include "CurveHold.h"
#include <fstream>
#include "CostFuncs.h"
class Calculator {
public:
	virtual void calc() = 0;
};
class Calculater : public Calculator {
public:
	void calc() override {
		std::cout << "Calculater " << std::endl;
	}
};
class Calculatir : public Calculator {
public:
	void calc() override {
		std::cout << "Calculatir " << std::endl;

	}
};

void writeData(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename);

int main()
{
	//std::cout << -p1;

	Calculator* c1 = new Calculater();
	Calculator* c2 = new Calculatir();

	std::vector<Calculator*> calcs;

	calcs.push_back(c1);
	calcs.push_back(c2);

	calcs[0]->calc();
	calcs[1]->calc();

	Eigen::Vector3d pInitial{ 0,0,0 };
	double theta = 0;
	double phi = 0;
	Eigen::Vector3d pTarget{ 100,100,500 };
	double R = 100;
	TrajectoryTemplate * ch1 = new CurveHold(pInitial, pTarget, theta, phi, R);
	ch1->fit();
	ch1->points(CoordinateSystem::CARTESIAN);
	ch1->points(CoordinateSystem::MD);
	//writeData(ch1->pointsCartesian, "output.txt");

	Eigen::Vector3d pInitial2{ 0,15,0 };
	Eigen::Vector3d pTarget2{ -100,-100,500 };
	TrajectoryTemplate* ch2 = new CurveHold(pInitial2, pTarget2, 0, 0, 100);
	ch2->fit();
	ch2->points(CoordinateSystem::CARTESIAN);
	ch2->points(CoordinateSystem::MD);
	double ch1DDI = DDI(ch1->pointsCartesian, ch1->pointsMd, false, 10);
	std::cout << "DDI is: " << ch1DDI << std::endl;
	std::vector<TrajectoryTemplate*> tT;
	tT.push_back(ch1);
	std::cout << "allLength: " << allLength(tT) << std::endl;
	
	double sf = sepFactor(ch1->pointsCartesian, ch1->pointsMd, ch2->pointsCartesian, false, 100);
	std::cout << "sepFactor: " << sf << std::endl;

	return 0;
}

void writeData(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename) {
	std::ofstream output;
	output.open(filename);
	for (size_t i = 0; i < pointsCartesian.size(); ++i) {
		output << pointsCartesian[i][0] << "," << pointsCartesian[i][1] << "," << pointsCartesian[i][2] << std::endl;
	}
	output.close();
}
