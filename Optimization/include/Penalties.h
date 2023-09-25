#pragma once
#include "TrajectoryTemplate.h"

struct Layer {
	double TVD, theta, phi;
};

struct Constraint {
	Layer lMin, lMax;
};

struct WellTrajectoryConstraints {
	double minDepthFirstHold;
	double maxDLS;
	double maxMD;
	double maxDistEastWest;
	double maxDistNorthSouth;
};


double PenaltyDLS(std::vector<TrajectoryTemplate*>& Well, double penalty);

double PenaltyIncWell(const std::vector<Eigen::Vector4d>& pMD, double incMin,
                      double incMax, double penalty);

double PenaltyIncTarget(std::vector<TrajectoryTemplate*>& Well, double penalty);

double PenaltyAlphaTarget(std::vector<TrajectoryTemplate*>& Well, double penalty);

double PenaltyConstraint(const std::vector<Constraint>& cs,
                         const std::vector<Eigen::Vector3d>& pC,
                         const std::vector<Eigen::Vector4d>& pMD, double penalty);

double PenaltyLength(double length, double MaxLength, double penalty = 100);

double PenaltyAHDNSEW(const std::vector<Eigen::Vector3d>& pC, double EWMAX,
                      double NSMAX, double penalty = 100);
