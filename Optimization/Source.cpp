#include "Eigen/Dense"
#include <iostream>
#include <vector>
#include "TrajectoryTemplate.h"
#include "CurveHold.h"
#include "CurveHoldCurveHold.h"
#include "Hold.h"
#include <fstream>
#include "CostFuncs.h"
#include "PSO.h"
#include <cmath>


//std::cout << -p1;

	//Calculator* c1 = new Calculater();
	//Calculator* c2 = new Calculatir();

	//std::vector<Calculator*> calcs;

	//calcs.push_back(c1);
	//calcs.push_back(c2);

	//calcs[0]->calc();
	//calcs[1]->calc();
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

void writeDataCartesian(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename);
void writeDataMD(std::vector<Eigen::Vector4d>& pointsMD, std::string filename);
double Rastrigin(Eigen::VectorXd x) {
	size_t n = x.size();
	double result = 10. * n;
	for (size_t i = 0; i < n; ++i)
		result += x(i) * x(i) - 10 * cos(2 * PI * x(i));
	return result;

}

double Sphere(Eigen::VectorXd x) {
	size_t n = x.size();
	double result = 0.;
	for (size_t idx = 0; idx < n; ++idx)
		result += x(idx) * x(idx);
	return result;
}

void getOptData(PSOvalueType op);


int main()
{
	double m = 1800 / PI;
	std::vector<double> minValues = { -1000., -1000., 100., 60., 0., 100. };
	std::vector<double> maxValues = { 1000., 1000., 2000., 70., 360., 1000. };
	std::vector<double> inert(100, .9);
	// x = x,y,z,inc,azi,beta
	std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd& x)> Well1 = [&](Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> well1;
		Eigen::Vector3d pIHold = { 0,0,0 };
		Eigen::Vector3d pTHold = { 0,0,x(5)};
		Eigen::Vector3d pChch1 = {x(0),x(1),x(2) };
		Eigen::Vector3d pT1Chch2 = { 400.0,8.0,3000.0 };
		Eigen::Vector3d pT3Chch2 = {1500 ,8 ,3010 };
		well1.push_back(new Hold(pIHold, pTHold));
		well1.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pChch1, x[3], x[4], 110.));
		//well1.push_back(new CurveHoldCurveHold(pChch1,x[3], x[4], m, m, pT3Chch2, 89.479, 0., 1100.04545));
		well1.push_back(new CurveHoldCurveHold(pChch1, x[3], x[4], m, m, pT1Chch2, pT3Chch2));
		return well1;
	};
	std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd& x)> Well2 = [&](Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> well2;
		Eigen::Vector3d pIHold = { 0,15,0 };
		Eigen::Vector3d pTHold = { 0,15,x[5] };
		Eigen::Vector3d pChch1 = { x[0],x[1],x[2] };
		Eigen::Vector3d pTChch2 = { -1500,8,3010 };
		well2.push_back(new Hold(pIHold, pTHold));
		well2.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pChch1, x[3], x[4], 110.));
		well2.push_back(new CurveHoldCurveHold(pChch1, x[3], x[4], m, m, pTChch2, 89.28, 180., 800.06));
		return well2;
	};
	std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd& x)> Well3 = [&](Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> well3;
		Eigen::Vector3d pIHold = { 0,20,0 };
		Eigen::Vector3d pTHold = { 0,20,x[5] };
		Eigen::Vector3d pChch1 = { x[0],x[1],x[2] };
		Eigen::Vector3d pTChch2 = { -1500,420,3010 };
		well3.push_back(new Hold(pIHold, pTHold));
		well3.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pChch1, x[3], x[4], 110.));
		well3.push_back(new CurveHoldCurveHold(pChch1, x[3], x[4], m, m, pTChch2, 89.326, 180., 850.06));
		return well3;
	};
	std::vector < std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd&)>> mainWell;
	std::vector < std::vector<TrajectoryTemplate*>> trajectories;
	mainWell.push_back(Well1);
	mainWell.push_back(Well2);
	mainWell.push_back(Well3);
	Eigen::VectorXd arg{{-381.313, -82.6248, 1523.74, 60.8638, 192.592, 839.657}};
	std::vector<TrajectoryTemplate *> well = Well1(arg);
	Eigen::Vector3d pIchch = { 0,0,0 };
	Eigen::Vector3d PTchch = { 100,100,1000};
	TrajectoryTemplate* chch = well[2];
	chch->fit();
	chch->points(CoordinateSystem::CARTESIAN);
	std::vector<Eigen::Vector3d>pCchch1 = chch->pointsCartesian;
	int cond = solve(well);
	std::vector<Eigen::Vector3d> pS = allPointsCartesian(well);
	writeDataCartesian(pCchch1,"output.txt");
	
	std::time_t start = time(NULL);
	/*
	for (size_t i = 0; i < mainWell.size(); ++i) {
		std::function<double(Eigen::VectorXd)> score = [&](Eigen::VectorXd x) {
			std::vector<TrajectoryTemplate*> tmp = mainWell[i](x);
			return orderScore(tmp, trajectories); };
		PSOvalueType opt = PSO(score, minValues, maxValues,7 , 6, inert);
		trajectories.push_back(mainWell[i](opt.first));
		getOptData(opt);
	}
	std::cout << "Time: " << time(NULL) - start << std::endl;
	*/
	return 0;
}

void getOptData(PSOvalueType op) {
	std::cout << "gBestCost: " << op.second << std::endl;
	std::cout << "gBestPos: ";
	for (auto x : op.first)
		 std::cout<< x << ", ";
	std::cout << std::endl;
}

void writeDataCartesian(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename) {
	std::ofstream output;
	output.open(filename);
	for (size_t i = 0; i < pointsCartesian.size(); ++i) {
		output << pointsCartesian[i][0] << "," << pointsCartesian[i][1] << "," << pointsCartesian[i][2] << std::endl;
	}
	output.close();
}

void writeDataMD(std::vector<Eigen::Vector4d>& pointsMD, std::string filename) {
	std::ofstream output;
	output.open(filename);
	for (size_t i = 0; i < pointsMD.size(); ++i) {
		output << pointsMD[i][0] << "," << pointsMD[i][1] << "," << pointsMD[i][2] <<"," << pointsMD[i][3] << std::endl;
	}
	output.close();
}
