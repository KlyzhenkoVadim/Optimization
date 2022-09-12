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
#include "API.h"
#include "Curve.h"
#include "TestWells.h"

struct Layer {
	double TVD, theta, phi;
};

struct Constraint {
	Layer lMin, lMax;
};


double OneWellScore(std::vector<TrajectoryTemplate*>& mainWell, double penalty = 1000);



void ShowOptData(PSOvalueType opt, std::vector<Constraint> cs);

int findTVD(const std::vector<Eigen::Vector3d>& pC, double TVD);

std::vector<TrajectoryTemplate*> wellCHCH(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target,const Constraint& c);

std::vector<TrajectoryTemplate*> Well(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs);

void OptimizeWells(const Eigen::Vector3d& pinit,const std::vector<Eigen::Vector3d>& targets,const std::vector<Constraint>& cs,
	const std::vector<double>& minValues, const std::vector<double>& maxValues, std::vector<std::vector<Eigen::Vector3d>>& pCWells,
	std::vector<std::vector<Eigen::Vector4d>>& pMDWells, std::vector<PSOvalueType>& opts);

void OptimizeCHCHWells(const Eigen::Vector3d& pinit, const std::vector<Eigen::Vector3d>& targets, const std::vector<Constraint>& cs,
	const std::vector<double>& minValues, const std::vector<double>& maxValues, std::vector<std::vector<Eigen::Vector3d>>& pCWells,
	std::vector<std::vector<Eigen::Vector4d>>& pMDWells, std::vector<PSOvalueType>& opts);

void FindBestWellHead(const std::vector<Eigen::Vector3d>& targets);

void testFindTVD(const std::vector<Eigen::Vector3d>& pC,const std::vector<Eigen::Vector4d>& pMD, double TVD);

void testWellCHCH(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target,
	std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD);

void testPenaltyConstraints(const std::vector<Constraint>& cs, std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD);

void testCH(double dls, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target);

void testScoreCHCH(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs,
	const std::vector<double>& minValues, const std::vector<double>& maxValues, bool argKnow = true);

void testNewCurveHold(const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, double inc, double azi, size_t nums);

void testWell(const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const Eigen::VectorXd& x, const std::vector<Constraint>& cs,
	std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD);

void testPenaltyInc(double incHold, double aziHold, double incMin, double incMax);

void testWellSolver(const Eigen::VectorXd& x, const Eigen::Vector3d& pInit, const Eigen::Vector3d& target);

void testSolver(const Eigen::VectorXd& x, const Eigen::Vector3d& pInit, const Eigen::Vector3d& target,
	const WellTrajectoryConstraints& cs, std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD);

void TestAPI(const Eigen::Vector3d& pInit, const Eigen::Vector3d& target, const WellTrajectoryConstraints& wc);

