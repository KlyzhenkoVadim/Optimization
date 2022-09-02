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
	void change() {
		lMax.phi = (lMax.phi + 180.);// -360 > EPSILON ? lMax.phi - 180 : lMax.phi + 180;
		lMin.phi = (lMin.phi + 180.);// -360 > EPSILON ? lMin.phi - 180 : lMin.phi + 180;
	}
};

double AzimuthChooser(Constraint cs, double r);

double PenaltyLength(double length, double penalty = 100);

double PenaltyDLS(std::vector<TrajectoryTemplate*>& Well, double penalty = 100);

double PenaltyIncTarget(std::vector<TrajectoryTemplate*>& Well, double penalty = 7);

double PenaltyAlphaTarget(std::vector<TrajectoryTemplate*>& Well, double penalty = 7);

std::vector<TrajectoryTemplate*> Well(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs);

Eigen::Vector3d OptimizeWellHead(std::vector<double> minValuesWellHead, std::vector<double> maxValuesWellHead);