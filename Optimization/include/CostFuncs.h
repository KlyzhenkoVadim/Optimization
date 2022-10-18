#pragma once
#include "TrajectoryTemplate.h"
#include <cmath>
#include <iostream>
#include "Penalties.h"

double sepFactor(std::vector<Eigen::Vector3d>& pCartesianW1,std::vector<Eigen::Vector4d>& pMDW1,
				std::vector<Eigen::Vector3d>& pCartesianW2, std::vector<Eigen::Vector4d>& pMDW2,
				double TVDstart = 0 ,bool actFunc = true , double penalty = 20);

double DDI(std::vector<Eigen::Vector3d>& pCartesian, std::vector<Eigen::Vector4d>& pMD, bool actFunc = true, double penalty = 5);

double orderScore1(std::vector<TrajectoryTemplate*>& mainWell, std::vector<std::vector<Eigen::Vector3d>>& pCTrajectories,
	std::vector<std::vector<Eigen::Vector4d>>& pMDTrajectories,double SepFactorShift = 0,double penalty = 1000);

double scoreSolver(std::vector<TrajectoryTemplate*>& tmp, const WellTrajectoryConstraints& cs, double penalty = 1000);
