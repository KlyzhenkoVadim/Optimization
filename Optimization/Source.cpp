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
#include "WellTrajectorySolver.h"
#include "Curve.h"
#include "TestWells.h"
#include "OptimizeWells.h"
#include <filesystem>





enum class platformNum { four, five, two };
namespace fs = std::filesystem;
const std::string path = "Results/";

void outputWellOptimization(const std::vector<PSOvalueType>& opts, std::vector<std::vector<Eigen::Vector3d>>& pCWells, std::vector<std::vector<Eigen::Vector4d>>& pMDWells,
	std::vector<std::string> names);

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

void setTargetsSlotting(std::vector<Eigen::Vector3d>& targets, platformNum num, bool newData)
{
	if (num == platformNum::four)
	{
		Eigen::Vector3d target40R = { 5803236,682857,2900 },target4001 = { 5803529,682498,2900 },
			target4003 = { 5803536,683257,2900 }, target4000 = { 5803409,683700,2900 };
		if (newData) {
			target40R[2] = 2462;
			target4001[2] = 2468;
			target4003[2] = 2462;
			target4000[2] = 2463;
		}
		targets = { target40R,target4001,target4003,target4000 };
	}
	if (num == platformNum::five)
	{
		Eigen::Vector3d target5008 = { 5804725,684117,2900 }, target5007 = { 5805418,683633,2900 },
			target50R = { 5806589,684341,2900 }, target5006 = { 5806299,684056,2900 }, target5PO = { 5805071,683838,2900 };
		if (newData)
		{
			target5008[2] = 2453;
			target5007[2] = 2446;
			target50R[2] = 2451;
			target5006[2] = 2452;
			target5PO[2] = 2448;
		}
		targets = {target5008,target5007,target50R,target5006,target5PO};
	}
	if (num == platformNum::two)
	{
		Eigen::Vector3d target4002 = { 5802919,680249,2900 }, target4PO = { 5802826,680679,2900 };
		if (newData)
		{
			target4002[2] = 2468;
			target4PO[2] = 2457;
		}
		targets = { target4002,target4PO };
	}
}

int main(int argc, char* argv[])
{
	std::string filename = argv[1];//"C:/Users/klyzhenko.vs/0CET/input.json";//
	Eigen::Vector3d pinit2 = { 5802508,680718,0 }, pinit4 = { 5802898,683790,0 }, pinit5{ 5804519,684699,0 };
	std::vector<Constraint> cs;// = { { {900,30,110},{900,90,170} }, { {2000,0,20},{2000,50,80} },{ {2600,0,20},{2600,40,80} } };
	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector < std::vector<Eigen::Vector4d>> pMDWells;
	std::vector<Eigen::Vector3d>  targets;
	std::vector<PSOvalueType> opts;
	std::vector<double> minValues, maxValues;
	Eigen::Vector3d pinit;
	AllConstraints ac;
	
	
	if (!fs::is_directory(path) || !fs::exists(path)) {
		fs::create_directory(path);
	}
	

	std::vector<std::string> names;
	parseData(pinit, targets, ac, minValues, maxValues, names,filename);
	OptimizeCHCHWells(pinit, targets, ac, minValues, maxValues, pCWells, pMDWells, opts);
	outputWellOptimization(opts, pCWells, pMDWells, names);
	return 0;
}

void outputWellOptimization(const std::vector<PSOvalueType>& opts, std::vector<std::vector<Eigen::Vector3d>>& pCWells, std::vector<std::vector<Eigen::Vector4d>>& pMDWells,
	std::vector<std::string> names)
{
	for (size_t i = 0; i < opts.size(); ++i) {
		getOptData(opts[i]);
	}
	std::string tmp = path + "cartesian", 
		tmpInc = path + "inclinometry";
	for (size_t i = 0; i < pCWells.size(); ++i) {
		writeDataCartesian(pCWells[i], tmp + names[i] + ".txt");
		writeInclinometry(pMDWells[i], tmpInc + names[i] + ".csv");
	}
}