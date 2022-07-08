#pragma once
#include "Eigen\Dense"
#include <random>
#include <iostream>

using PSOvalueType =  std::pair<Eigen::VectorXd , double>;

PSOvalueType PSO(std::function<double(Eigen::VectorXd)> func, std::vector<double>& minValues, std::vector<double>& maxValues,
	size_t numAgents, size_t dimension,std::vector<double>& inertia, double socCoef = 0.3, double indCoef = 0.5, size_t numIterations = 100);
