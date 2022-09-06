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

void setMinMaxValuesWell(std::vector<double>& minValues, std::vector<double>& maxValues,const std::vector<Constraint> cs) {
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
}

int main()
{
	double R = 1200 / PI;
	Eigen::Vector3d pinit = { 5803419,683949,0 }, target1 = { 5803236,682857,2900 },
		target2 = { 5803529,682498,2900 }, target3 = { 5803536,683257,2900 }, target4 = { 5803409,683700,2900 }; //40R 4001 4003 4000
	std::vector<Constraint> cs = { { {900,30,110},{900,90,170} },
		{ {2000,0,20},{2000,50,80} },{ {2600,0,20},{2600,40,80} } };
	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector < std::vector<Eigen::Vector4d>> pMDWells;
	std::vector<Eigen::Vector3d>  targets = { target1,target2,target3,target4 };
	std::vector<PSOvalueType> opts;
	std::vector<double> minValues{400,0,0,0,0}, maxValues{900,1.5,1.5,40,360};
	std::vector<Eigen::Vector3d>pC;
	std::vector<Eigen::Vector4d>pMD;
	Eigen::VectorXd x { {690.35,0.63081,0.536736,15.8902,80.8191} };
	testWellCHCH(x, pinit, target1, pC, pMD);
	//testFindTVD(pC, pMD, 2600);
	//testPenaltyConstraints(cs, pC, pMD);
	testScoreCHCH(x,pinit, target1, cs, minValues,maxValues, true);
	targets.clear();
	targets.push_back(target1);
	//targets.push_back(target2);
	//OptimizeCHCHWells(pinit, targets, cs, minValues, maxValues, pCWells, pMDWells, opts);
	//outputWellOptimization(opts, pCWells, pMDWells,targets, cs);

	return 0;
}

void outputWellOptimization(const std::vector<PSOvalueType>& opts, std::vector<std::vector<Eigen::Vector3d>>& pCWells, std::vector<std::vector<Eigen::Vector4d>>& pMDWells, 
	std::vector<Eigen::Vector3d>&targets, std::vector<Constraint>& cs) 
{
	for (size_t i = 0; i < targets.size(); ++i) {
		//ShowOptData(opts[i], cs);
		getOptData(opts[i]);
	}
	std::vector<std::string> names = { "40R", "4001", "4003", "4000" };
	std::string tmp = "outputWell(", tmpInc = "inclinometry";
	for (size_t i = 0; i < pCWells.size(); ++i) {
		writeDataCartesian(pCWells[i], tmp + names[i] + ").txt");
		writeInclinometry(pMDWells[i], tmpInc + names[i] + ").txt");
	}
}