#pragma once
#include<iomanip>
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

void writeInclinometry(const std::vector<Eigen::Vector4d>& pMD, std::string filename);

void writeDataCartesian(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename);

void writeDataMD(std::vector<Eigen::Vector4d>& pointsMD, std::string filename);

void getOptData(PSOvalueType op);

void writeDataOpt(std::vector<size_t> order, size_t wellNum, PSOvalueType Opt, std::vector<std::vector<TrajectoryTemplate*>>& trajs);

void writeDataOptSep(std::vector<std::vector<Eigen::Vector3d>>& pCTrajs, std::vector<std::vector<Eigen::Vector4d>>& pMDTrajs, size_t& num);
