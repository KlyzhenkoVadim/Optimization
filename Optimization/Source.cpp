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
#include "API.h"
#include "Curve.h"
#include "TestWells.h"
#include "OptimizeWells.h"


void outputWellOptimization(const std::vector<PSOvalueType>& opts, std::vector<std::vector<Eigen::Vector3d>>& pCWells, std::vector<std::vector<Eigen::Vector4d>>& pMDWells,
	std::vector<Eigen::Vector3d>& targets, std::vector<Constraint>& cs);

void setMinMaxValuesWell(std::vector<double>& minValues, std::vector<double>& maxValues, const std::vector<Constraint> cs) {
	minValues.push_back(400); // 100
	maxValues.push_back(900); // 400
	for (size_t i = 0; i < 3; ++i) {
		minValues.push_back(cs[i].lMin.theta);
		minValues.push_back(0); //minValues.push_back(cs[i].lMin.phi);
		maxValues.push_back(cs[i].lMax.theta); //maxValues.push_back(cs[i].lMax.phi);
		maxValues.push_back(1);
	}
	minValues.push_back(0); // 400
	maxValues.push_back(1); // 900
	for (size_t j = 1; j < 3; ++j) {
		minValues.push_back(cs[j - 1].lMin.TVD);
		maxValues.push_back(cs[j].lMax.TVD);
	}
	minValues.push_back(0); // R
	maxValues.push_back(1.5); // 1e8
};

int main()
{
	double R = 1200 / PI; //For targets4 THE BEST 5802898,683790,0 // 5804304, 682287.9,0
	Eigen::Vector3d pinit2 = {5802508,680718,0}, pinit4 = { 5802898,683790,0 }, pinit5{ 5804519,684699,0 }, target40R = { 5803236,682857,2900 },//  5805450, 682009.6,0 - GOOD  //5802610, 681358
		target4001 = { 5803529,682498,2900 }, target4003 = { 5803536,683257,2900 }, target4000 = { 5803409,683700,2900 },
		target5008 = { 5804725,684117,2900 }, target5007 = { 5805418,683633,2900 }, target50R = { 5806589,684341,2900 },
		target5006 = { 5806299,684056,2900 }, target5PO = { 5805071,683838,2900 },
		target4002 = { 5802919,680249,2900 }, target4PO = { 5802826,680679,2900 };
		//5802987,680237 target2s
		//40R 4001 4003 4000
	std::vector<Constraint> cs = { { {900,30,110},{900,90,170} }, { {2000,0,20},{2000,50,80} },{ {2600,0,20},{2600,40,80} } };
	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector < std::vector<Eigen::Vector4d>> pMDWells;
	std::vector<Eigen::Vector3d>  targets = { target40R,target4001,target4003,target4000 }, target5s = {target5008,target5007,target50R,target5006,target5PO};
	std::vector<Eigen::Vector3d> target2s = { target4002,target4PO };
	std::vector<PSOvalueType> opts;
	std::vector<double> minValues{400,0,0,0,0,0}, maxValues{900,1.5,1.5,40,360,2900};
	std::vector<Eigen::Vector3d>pC;
	std::vector<Eigen::Vector4d>pMD;
	Eigen::VectorXd x{{}}; // 5804304, 682287.9
	//OptimizeCHCHWells(pinit2, target2s, cs, minValues, maxValues, pCWells, pMDWells, opts);
	Eigen::Vector3d pInit{ 0,0,0 }, target{ 100,500,1000 };
	WellTrajectoryConstraints wc{400,1.5,5000,2000,2000};
	TestAPI(pinit4, target40R, wc);
	//testSolver(x,pInit, target, wc, pC, pMD);
	//outputWellOptimization(opts, pCWells, pMDWells, target2s, cs);
	return 0;
}

void outputWellOptimization(const std::vector<PSOvalueType>& opts, std::vector<std::vector<Eigen::Vector3d>>& pCWells, std::vector<std::vector<Eigen::Vector4d>>& pMDWells, 
	std::vector<Eigen::Vector3d>&targets, std::vector<Constraint>& cs) 
{
	for (size_t i = 0; i < targets.size(); ++i) {
		//ShowOptData(opts[i], cs);
		getOptData(opts[i]);
	}
	std::vector<std::string>  names, name5s = { "5008", "5007", "50R", "5006", "5PO"},name4s = { "40R", "4001", "4003", "4000" },
		name2s = {"4002","4PO"};//
	names = name2s;
	std::string tmp = "outputWell(", tmpInc = "inclinometry(";
	for (size_t i = 0; i < pCWells.size(); ++i) {
		writeDataCartesian(pCWells[i], tmp + names[i] + ").txt");
		writeInclinometry(pMDWells[i], tmpInc + names[i] + ").txt");
	}
}