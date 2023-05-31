#pragma once
#include "TrajectoryTemplate.h"
#include "CurveHoldCurveHold.h"
#include "Hold.h"
#include "CostFuncs.h"
#include <fstream>
#include "C:\Users\klyzhenko.vs\Desktop\ArchiveWProjects\optimization\Optimization\packages\nlohmann.json.3.11.2\build\native\include\nlohmann\json.hpp"
#include "PSO.h"
#include "Penalties.h"
#include "TestWells.h"

struct AllConstraints
{
	WellTrajectoryConstraints wc{ 0,1.5,1 / EPSILON / EPSILON,1 / EPSILON / EPSILON };
	std::vector<Constraint> cs;
	double MaxIncTarget = 90;
};

class DirectionalDrillingSolver
{
private:
	AllConstraints constraints;
	Eigen::Vector3d pinit;
	std::vector<Eigen::Vector3d> targets;
	std::vector<double> minValues, maxValues;
	std::vector<std::string> names;
	std::vector<TrajectoryTemplate*> wellfunction(const std::vector<double>& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target);
	std::vector<PSOvalueType> opts;
	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector<std::vector<Eigen::Vector4d>> pMDWells;

public:
	DirectionalDrillingSolver();
	void setData(std::string filename);
	void Optimize();
	void writeResults(std::string outputdir);

};
