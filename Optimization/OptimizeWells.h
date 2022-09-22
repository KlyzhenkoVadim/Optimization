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
#include "Penalties.h"


double OneWellScore(std::vector<TrajectoryTemplate*>& mainWell, double penalty = 1000);

std::vector<TrajectoryTemplate*> Well(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs);

std::vector<TrajectoryTemplate*> wellSolverHorizontal(const Eigen::VectorXd& x, const Eigen::Vector3d& pInit, const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3);

void OptimizeHorizontal(const Eigen::Vector3d& pinit, const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3);

void OptimizeHorizontals(const Eigen::Vector3d& pinit, const std::vector<Eigen::Vector3d>& targets1, const std::vector<Eigen::Vector3d>& targets3);

void testHorizontal(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3);

void OptimizeWells(const Eigen::Vector3d& pinit,const std::vector<Eigen::Vector3d>& targets,const std::vector<Constraint>& cs,
	const std::vector<double>& minValues, const std::vector<double>& maxValues, std::vector<std::vector<Eigen::Vector3d>>& pCWells,
	std::vector<std::vector<Eigen::Vector4d>>& pMDWells, std::vector<PSOvalueType>& opts);


void MakePointsEvo(const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::string& filename, const std::string& numtarget);