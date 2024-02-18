#pragma once
#include <fstream>

#include <nlohmann/json.hpp>

#include <CostFuncs.h>
#include <CurveHoldCurveHold.h>
#include <Hold.h>
#include <PSO.h>
#include <Penalties.h>
#include <TestWells.h>
#include <TrajectoryTemplate.h>

struct AllConstraints
{
	WellTrajectoryConstraints wc{0, 1.5, 1 / EPSILON / EPSILON,
                               1 / EPSILON / EPSILON};
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
    std::vector<TrajectoryTemplate*> wellfunction(
            const std::vector<double>& x, const Eigen::Vector3d& pinit,
            const Eigen::Vector3d& target);
	std::vector<PsoValueType> opts;
	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector<std::vector<Eigen::Vector4d>> pMDWells;

public:
	DirectionalDrillingSolver();
	void setData(const std::string& filename);
	void Optimize();
	void writeResults(const std::string& outputdir);

};