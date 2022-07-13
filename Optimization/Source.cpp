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
void writeDataOpt(std::vector<size_t> order, size_t wellNum, PSOvalueType Opt);

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
		Eigen::Vector3d pT3Chch2 = {1500. ,8.0 ,3010.0 };
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
		Eigen::Vector3d pT1Chch2 = { -700.,8.0,3000. };
		Eigen::Vector3d pT3Chch2 = { -1500.,8,3010. };
		well2.push_back(new Hold(pIHold, pTHold));
		well2.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pChch1, x[3], x[4], 110.));
		well2.push_back(new CurveHoldCurveHold(pChch1, x[3], x[4], m, m,pT1Chch2 ,pT3Chch2));
		//well2.push_back(new CurveHoldCurveHold(pChch1,x[3], x[4], m, m, pT3Chch2, 89.28384005452959, 180.0,800.0624975587845));
		return well2;
	};
	std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd& x)> Well3 = [&](Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> well3;
		Eigen::Vector3d pIHold = { 0,20,0 };
		Eigen::Vector3d pTHold = { 0,20,x[5] };
		Eigen::Vector3d pChch1 = { x[0],x[1],x[2] };
		Eigen::Vector3d pT1Chch2{ -650.,420.,3000. };
		Eigen::Vector3d pT3Chch2 = { -1500.,420.,3010. };
		well3.push_back(new Hold(pIHold, pTHold));
		well3.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pChch1, x[3], x[4], 110.));
		well3.push_back(new CurveHoldCurveHold(pChch1, x[3], x[4], m, m, pT1Chch2,pT3Chch2));
		//well3.push_back(new CurveHoldCurveHold(pChch1,x[3], x[4], m, m, pT3Chch2, 89.32596310201549, 180.0,850.0588214941364));

		return well3;
	};
	std::vector < std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd&)>> mainWell;
	std::vector < std::vector<TrajectoryTemplate*>> trajectories;
	mainWell.push_back(Well1);
	mainWell.push_back(Well2);
	mainWell.push_back(Well3);
	std::vector<std::vector<size_t>> order = {{1,2,3},{1,3,2},{2,1,3},{2,3,1},{3,1,2},{3,2,1}};
	/*for (size_t idx = 0; idx < 6; ++idx) {
		for (auto ord : order[idx]) {
			std::function<double(Eigen::VectorXd)> score = [&](Eigen::VectorXd x) {
				std::vector<TrajectoryTemplate*> tmp = mainWell[ord - 1](x);
				return orderScore(tmp, trajectories); };
			PSOvalueType opt = PSO(score, minValues, maxValues, 12, 6, inert,0.3,0.5,100);
			trajectories.push_back(mainWell[ord - 1](opt.first));
			int cond = solve(trajectories.back());
			if (cond > 0) {
				break;
				std::cout << "Unbuild Trajectory\n";
			}
			writeDataOpt(order[idx], ord, opt);
		}
		trajectories.clear();
	}*/
	Eigen::VectorXd arg1{ {-346.702, 190.616, 1582.23, 60.6974 ,149.54, 974.232} },
		arg2{ {84.8843, -411.808, 856.599, 64.9406, 279.866 ,247.553} }, arg3{ {266.87, 324.09, 1470.35, 60.8709, 51.2091, 862.812} };
	std::vector<std::vector<TrajectoryTemplate*>> wells = { mainWell[0](arg1),mainWell[1](arg2),mainWell[2](arg3) };
	for (size_t i = 0; i < wells.size(); ++i)
		int cond = solve(wells[i]);
	std::vector<Eigen::Vector3d> pCartesian = allPointsCartesian(wells[2]);
	writeDataCartesian(pCartesian, "output.txt");
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

void writeDataOpt(std::vector<size_t> order,size_t wellNum, PSOvalueType Opt) {
	std::fstream file;
	file.open("C:/Users/klyzhenko.vs/0CET/optDataCpp.csv", std::ios::out | std::ios::app);
	file << "[";
	for (size_t i = 0; i < order.size(); ++i) {
		file << order[i];
		if(i!=order.size()-1)
			file << " ";
	}
	file << "]," << wellNum << ",[";
	for (size_t i = 0; i < Opt.first.size(); ++i) {
		file << Opt.first[i];
		if (i != Opt.first.size() - 1)
			file << " ";
	}
	file << "]," << Opt.second << std::endl;
	file.close();
}