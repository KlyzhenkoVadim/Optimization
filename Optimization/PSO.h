#pragma once
#include "Eigen\Dense"
#include <random>
#include <iostream>

using PSOvalueType =  std::pair<Eigen::VectorXd , double>;

PSOvalueType PSO(std::function<double(const Eigen::VectorXd&)> func, const std::vector<double>& minValues, const std::vector<double>& maxValues,
	size_t numAgents, size_t dimension,size_t numIterations = 100, const std::vector<double>& inertia= std::vector<double>(500,0.9),
	double socCoef = 0.3, double indCoef = 0.5);
