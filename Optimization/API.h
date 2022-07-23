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

class someAPI {
private:
	std::vector<Eigen::Vector3d> Wellheads, Targets1, Targets3; // ¬екторы координат усть€ и целей с учЄтом пор€дка !!!
	//double intensivity, R; // интенсивность набора угла.
	std::vector<Eigen::VectorXd> args;
	std::vector<TrajectoryTemplate*> trajectories;
	std::vector<Eigen::Vector3d> pCtrajectories;
	std::vector<Eigen::Vector4d> pMDtrajectories;
	std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd& x)> typeWell;
	std::string type;
public:
	someAPI(const std::string type);
	void setData();
	void setPSOdata();
	void getData();
	void getPSOdata();
	void Optimize();
	void TypeTrajectory();
};