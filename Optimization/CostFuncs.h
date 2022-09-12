#pragma once
#include "TrajectoryTemplate.h"
#include <cmath>
#include <iostream>


struct WellTrajectoryConstraints {
	double minDepthFirstHold;
	double maxDLS;
	double maxMD;
	double maxDistEastWest;
	double maxDistNorthSouth;
};

double sigmoid(double x, double penalty, double alpha, double x0);

int signum(double x);

double PenaltyLength(double length, double MaxLength,double penalty = 100);

double PenaltyAHDNSEW(const std::vector<Eigen::Vector3d>& pC, double EWMAX, double NSMAX, double penalty = 100);

double sepFactor(std::vector<Eigen::Vector3d>& pCartesianW1,std::vector<Eigen::Vector4d>& pMDW1,
				std::vector<Eigen::Vector3d>& pCartesianW2, std::vector<Eigen::Vector4d>& pMDW2,
				double TVDstart = 0 ,bool actFunc = true , double penalty = 20);

double AHD(std::vector<Eigen::Vector3d>& pX, std::vector<Eigen::Vector3d>& pY);

double dls(Eigen::Vector3d& tangent1, Eigen::Vector3d& tangent2);

double Tortuosity(std::vector<Eigen::Vector4d>& pointsMD);

double DDI(std::vector<Eigen::Vector3d>& pCartesian, std::vector<Eigen::Vector4d>& pMD, bool actFunc = true, double penalty = 5);

double orderScore1(std::vector<TrajectoryTemplate*>& mainWell, std::vector<std::vector<Eigen::Vector3d>>& pCTrajectories,
	std::vector<std::vector<Eigen::Vector4d>>& pMDTrajectories,double SepFactorShift = 0,double penalty = 1000);

double scoreSolver(std::vector<TrajectoryTemplate*>& tmp, const WellTrajectoryConstraints& cs, double penalty = 1000);