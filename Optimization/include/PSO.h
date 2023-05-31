#pragma once
#include "Eigen\Dense"
#include "C:\Users\klyzhenko.vs\Desktop\ArchiveWProjects\optimization\Optimization\packages\nlohmann.json.3.11.2\build\native\include\nlohmann\json.hpp"
#include <random>
#include <iostream>

namespace pso 
{
	constexpr double EPSILON = 1e-10;
} 
//using PSOvalueType =  std::pair<std::vector<double> , double>;

struct PSOvalueType
{
	std::vector<double> optVec;
	double cost;
	std::vector<double> costHist;
};

PSOvalueType PSO(std::function<double(const std::vector<double>&)> func, const std::vector<double>& minValues, const std::vector<double>& maxValues,
	size_t numAgents, size_t dimension,size_t numIterations = 100, const std::vector<double>& inertia= std::vector<double>(500,0.9),
	double socCoef = 0.3, double indCoef = 0.5,bool history = false);