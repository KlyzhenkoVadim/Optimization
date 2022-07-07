#pragma once
#include "TrajectoryTemplate.h"

double sigmoid(double x, double penalty, double alpha, double x0);

double sepFactor(std::vector<Eigen::Vector3d>& pCartesianW1,std::vector<Eigen::Vector4d>& pMDW1,
				std::vector<Eigen::Vector3d>& pCartesianW2, bool actFunc = true , double penalty = 100);

double AHD(std::vector<Eigen::Vector3d>& pX, std::vector<Eigen::Vector3d>& pY);

double dls(Eigen::Vector3d& tangent1, Eigen::Vector3d& tangent2);

double Tortuosity(std::vector<Eigen::Vector4d>& pointsMD);

double DDI(std::vector<Eigen::Vector3d>& pCartesian, std::vector<Eigen::Vector4d>& pMD, bool actFunc = true, double penalty = 10.);

double orderScore(std::vector<TrajectoryTemplate*>& mainWell, std::vector<std::vector<TrajectoryTemplate*>>& Trajectories, double penalty = 1000);
