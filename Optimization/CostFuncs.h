#pragma once
#include "TrajectoryTemplate.h"

double sigmoid(double x, double penalty, double alpha, double x0);

int signum(double x);

double sepFactor(std::vector<Eigen::Vector3d>& pCartesianW1,std::vector<Eigen::Vector4d>& pMDW1,
				std::vector<Eigen::Vector3d>& pCartesianW2, std::vector<Eigen::Vector4d>& pMDW2,
				double TVDstart = 0 ,bool actFunc = true , double penalty = 20);

double AHD(std::vector<Eigen::Vector3d>& pX, std::vector<Eigen::Vector3d>& pY);

double dls(Eigen::Vector3d& tangent1, Eigen::Vector3d& tangent2);

double Tortuosity(std::vector<Eigen::Vector4d>& pointsMD);

double DDI(std::vector<Eigen::Vector3d>& pCartesian, std::vector<Eigen::Vector4d>& pMD, bool actFunc = true, double penalty = 5);

double OneWellScore(std::vector<TrajectoryTemplate*>& mainWell, double penalty = 1000);

double orderScore1(std::vector<TrajectoryTemplate*>& mainWell, std::vector<std::vector<Eigen::Vector3d>>& pCTrajectories,
	std::vector<std::vector<Eigen::Vector4d>>& pMDTrajectories,double SepFactorShift = 0,double penalty = 1000);