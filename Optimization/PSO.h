#pragma once
#include "Eigen\Dense"
#include <random>

using PSOvalueType =  std::pair<Eigen::VectorXd , double>;

PSOvalueType PSO(std::function<double(Eigen::VectorXd)> func, std::vector<double> minValues, std::vector<double> maxValues,
	size_t numAgents, size_t dimension, double socCoef, double indCoef, size_t numIterations,  std::vector<double> inertia);
