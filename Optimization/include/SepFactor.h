#pragma once
#include "OptimizeWells.h"

double SeparationFactor(const std::vector<Eigen::Vector3d>& pC1,
	const std::vector<Eigen::Vector4d>& pMD1,
	const std::vector<Eigen::Vector3d>& pC2,
	const std::vector<Eigen::Vector4d>& pMD2,
	double TVDShift = 0);

void test_separationFactor();

void test_sepfactor_speed();