#pragma once
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
#include "C:\Users\klyzhenko.vs\Desktop\optimization\Optimization\packages\nlohmann.json.3.11.2\build\native\include\nlohmann\json.hpp"

struct Layer {
	double TVD, theta, phi;
};

struct Constraint {
	Layer lMin, lMax;
};

struct AllConstraints
{
	WellTrajectoryConstraints wc{0,1.5,1/EPSILON,1/EPSILON};
	std::vector<Constraint> cs;
	double MaxIncTarget = 90;
};

double OneWellScore(std::vector<TrajectoryTemplate*>& mainWell, double penalty = 1000);

void ShowOptData(PSOvalueType opt, std::vector<Constraint> cs);

int findTVD(const std::vector<Eigen::Vector3d>& pC, double TVD);

std::vector<TrajectoryTemplate*> wellCHCH(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target);

std::vector<TrajectoryTemplate*> Well(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs);

std::vector<TrajectoryTemplate*> wellSolverDirectional(const Eigen::VectorXd& x, const Eigen::Vector3d pInit, const Eigen::Vector3d& target);

std::vector<TrajectoryTemplate*> wellSolverHorizontal(const Eigen::VectorXd& x, const Eigen::Vector3d& pInit, const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3);

void OptimizeWells(const Eigen::Vector3d& pinit,const std::vector<Eigen::Vector3d>& targets,const std::vector<Constraint>& cs,
	const std::vector<double>& minValues, const std::vector<double>& maxValues, std::vector<std::vector<Eigen::Vector3d>>& pCWells,
	std::vector<std::vector<Eigen::Vector4d>>& pMDWells, std::vector<PSOvalueType>& opts);

void OptimizeCHCHWells(const Eigen::Vector3d& pinit, const std::vector<Eigen::Vector3d>& targets, const AllConstraints& constraints,
	const std::vector<double>& minValues, const std::vector<double>& maxValues, std::vector<std::vector<Eigen::Vector3d>>& pCWells,
	std::vector<std::vector<Eigen::Vector4d>>& pMDWells, std::vector<PSOvalueType>& opts);

void FindBestWellHead(const std::vector<Eigen::Vector3d>& targets);

void testPenaltyConstraints(const std::vector<Constraint>& cs, std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD);

void testScoreCHCH(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs,
	const std::vector<double>& minValues, const std::vector<double>& maxValues);

void TestAPI(const Eigen::Vector3d& pInit, const Eigen::Vector3d& target1, const Eigen::Vector3d& target3, const WellTrajectoryConstraints& wc);

void testWellTrajectory(std::vector<TrajectoryTemplate*>& tmp, std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD, const std::string& filename);

void parseData(Eigen::Vector3d& pinit, std::vector<Eigen::Vector3d>& targets, AllConstraints& constraints,
	std::vector<double>& minValues,  std::vector<double>& maxValues, std::vector<std::string>& names, const std::string& filename);
