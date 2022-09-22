#pragma once
#include<iomanip>
#include <iostream>
#include <fstream>
#include "PSO.h"
#include "TrajectoryTemplate.h"

void writeInclinometry(const std::vector<Eigen::Vector4d>& pMD, std::string filename);

void writeDataCartesian(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename);

void writeDataMD(std::vector<Eigen::Vector4d>& pointsMD, std::string filename);

void getOptData(PSOvalueType op);

