#pragma once
#include "TrajectoryTemplate.h"
#include <cmath>
#include <iostream>
#include "Penalties.h"

double sepFactor(const std::vector<Eigen::Vector3d>& pCartesianW1,
                 const std::vector<Eigen::Vector4d>& pMDW1,
                 const std::vector<Eigen::Vector3d>& pCartesianW2,
                 const std::vector<Eigen::Vector4d>& pMDW2, double TVDstart = 0,
                 bool actFunc = true, double penalty = 20);

double DDI(std::vector<TrajectoryTemplate*>& well,
           const std::vector<Eigen::Vector3d>& pCartesian, bool actFunc = true,
           double penalty = 5);

double ERD(const std::vector<Eigen::Vector3d>& points);

double orderScore1(
    std::vector<TrajectoryTemplate*>& mainWell,
    const std::vector<std::vector<Eigen::Vector3d>>& pCTrajectories,
    const std::vector<std::vector<Eigen::Vector4d>>& pMDTrajectories,
    double SepFactorShift = 0, double penalty = 1000);

double scoreSolver(std::vector<TrajectoryTemplate*>& tmp,
                   const WellTrajectoryConstraints& cs, double penalty = 1000);

double dls(Eigen::Vector3d& tangent1, Eigen::Vector3d& tangent2);

double Tortuosity(std::vector<TrajectoryTemplate*>& well);
