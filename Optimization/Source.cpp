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
#include "OptimizeWells.h"


void ShowOptData(PSOvalueType opt,std::vector<Constraint> cs) {
	std::cout << "Optimization result:\n_____________________________________\n";
	std::cout << "CostValue: " << opt.second << '\n';
	std::cout << "Hold0TVD: " << opt.first[0] << '\n';
	for (size_t i = 0; i < 3; ++i) {
		std::cout << "Curve" << i + 1 << " (inc" << i + 1 << ", azi" << i + 1 << "): " << opt.first[2 * i + 1] << ", " << AzimuthChooser(cs[i],opt.first[2 * (i + 1)]) << "\n";
		std::cout << "Hold" << i + 1 << "TVDstart: " << opt.first[i + 7] << "\n";
	}
	std::cout << "CurveHoldDLS: " << 1800 / PI / opt.first[10] << "\n";
}

std::vector<TrajectoryTemplate*> wellCHCH(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target) {
	std::vector<TrajectoryTemplate*> tmp;
	tmp.push_back(new Hold(pinit, 0, 0, x[0]));
	tmp.back()->getTarget1Point();
	double R1, R2;
	if (abs(x[1]) < EPSILON) {
		R1 = 1 / EPSILON;
	}
	if (abs(x[2]) < EPSILON) {
		R2 = 1 / EPSILON;
	}
	R1 = 1800 / PI / x[1];
	R2 = 1800 / PI / x[2];
	tmp.push_back(new CurveHoldCurveHold(tmp.back()->pointT1, 0, 0, R1, R2,target,x[3],x[4]));
	return tmp;
}

int main()
{
	double R = 1200 / PI;
	//Eigen::Vector3d pinit = { 5805626 - 2500,683999,0 };
	//5803868, 682438.2 // 5803148, 682359.4
	Eigen::Vector3d pinit = { 5805626 - 2500,683999 - 500,0 }, target1 = { 5803236,682857,2900 },
		target2 = { 5803529,682498,2900 }, target3 = { 5803536,683257,2900 }, target4 = { 5803409,683700,2900 };;
	std::vector<Constraint> cs = { { {900,30,110},{900,90,170} },
		{ {2000,0,200},{2000,50,260} },{ {2600,0,20},{2600,40,80} } };
	// 4001Init (NS,EW,TVD) : (5805626-500+15,683999-500,0)
	//4001Target:	(NS,EW,TVD) : (5803529,682498,2900)
	//Target4003:	5803536	683257	2900
	// Target4000	5803409	683700	2900
	
	std::vector<double> minValues,maxValues;
	minValues.push_back(400); // 100
	maxValues.push_back(700); // 400
	for (size_t i = 0; i < 3; ++i) {
		minValues.push_back(cs[i].lMin.theta);
		minValues.push_back(0); //minValues.push_back(cs[i].lMin.phi);
		maxValues.push_back(cs[i].lMax.theta); //maxValues.push_back(cs[i].lMax.phi);
		maxValues.push_back(1);
	}
	minValues.push_back(0); // 400
	maxValues.push_back(1); // 900
	for (size_t j = 1; j < 3; ++j) {
		minValues.push_back(cs[j - 1].lMin.TVD);
		maxValues.push_back(cs[j].lMax.TVD);
	}
	minValues.push_back(R);
	maxValues.push_back(1500);

	std::function<double(const Eigen::VectorXd&)> scoreCHCH = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> tmp = wellCHCH(x, pinit, target1);
		int s = solve(tmp);
		std::vector<Eigen::Vector3d> pCtmp = allPointsCartesian(tmp);
		std::vector<Eigen::Vector4d>  pMDtmp = allPointsMD(tmp);
		return 0;
	};

	std::function<double(const Eigen::VectorXd&)> scoreOne = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> tmp = Well(x,pinit,target1,cs);
		int condition = solve(tmp);
		if (condition != 0) {
			return 1000. * condition / tmp.size(); // penalty * percent of incorrect templates.
		}
		double  length = allLength(tmp), dlspenalty = PenaltyDLS(tmp);
		double IdealLength = length ;
		tmp[0]->getInitPoint();
		tmp.back()->getTarget1Point();
		tmp.back()->getTarget3Point();
		IdealLength = (tmp.back()->pointT1 - tmp[0]->pointInitial).norm() + (tmp.back()->pointT3 - tmp.back()->pointT1).norm();
		for (auto x : tmp) {
			delete x;
		}
		return length/IdealLength + PenaltyLength(length) + dlspenalty;
	};

	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector < std::vector<Eigen::Vector4d>> pMDWells;
	Eigen::Vector3d target = target1;
	double TVDShift = 0;
	std::vector<Eigen::Vector3d>  targets = { target1,target2,target3,target4 };
	std::function<double(const Eigen::VectorXd&)> CollRiskScore = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> currWell = Well(x, pinit, target, cs);
		double score = orderScore1(currWell, pCWells, pMDWells,TVDShift);
		double length = allLength(currWell), dlsPen = PenaltyDLS(currWell,1000), incPen = PenaltyIncTarget(currWell,20);
		double alphaCH = PenaltyAlphaTarget(currWell,100);
		for (auto x : currWell)
			delete x;
		return score + PenaltyLength(length, 100) + dlsPen;// +incPen;
	};

	bool StartOptimization;
	std::cout << "Start Optimization?";
	std::cin >> StartOptimization;
	if (StartOptimization) {
		for (size_t idd = 0; idd < 4; ++idd) {
			PSOvalueType opt = PSO(CollRiskScore, minValues, maxValues, 30, 11, 200);
			ShowOptData(opt, cs);
			getOptData(opt);
			std::vector<TrajectoryTemplate*> well = Well(opt.first,pinit,target,cs);
			int cond = solve(well);
			double sc = CollRiskScore(opt.first);
			pCWells.push_back(allPointsCartesian(well));
			pMDWells.push_back(allPointsMD(well));
			TVDShift = std::max(TVDShift, opt.first[0]);
			target = targets[idd+1];
			std::cout << "Length: " << allLength(well) << "\n";
		}
		writeDataCartesian(pCWells[0], "outputWell(A).txt");
		writeInclinometry(pMDWells[0], "inclinometry40R.txt");
		writeDataCartesian(pCWells[1], "outputWell(B).txt");
		writeInclinometry(pMDWells[1], "inclinometry4001.txt");
		writeDataCartesian(pCWells[2], "outputWell(C).txt");
		writeInclinometry(pMDWells[2], "inclinometry4003.txt");
		writeDataCartesian(pCWells[3], "outputWell(D).txt");
		writeInclinometry(pMDWells[3], "inclinometry4000.txt");
	}
	// Hold1 inc1 azi1 inc2 azi2 inc3 azi3 tvd1 tvd2 tvd3 R
	std::vector<TrajectoryTemplate*> tTest;
	
	Eigen::VectorXd arg{ { 169.974, 30.1904, 0.591368, 49.7777, 0.00832002, 23.6241, 0.565303, 763.154, 1358.08, 2427.16, 418.177} };
	if (!StartOptimization) {
		//NS: (5803236, 5803536)
		//EW: (682498, 683700)
		/*std::vector<double> minVs{5803236,682498}, maxVs{5803536 , 683700};
		Eigen::Vector3d wellhead = OptimizeWellHead(minVs, maxVs);
		for (auto x : wellhead) {
			std::cout << x << ",";
		}*/
		std::vector<double> minVs{ 5800236,680000}, maxVs{ 5810000,685000};
		std::function<double(const Eigen::VectorXd&)> headScore = [&](const Eigen::VectorXd& x) {
			Eigen::Vector3d InitialPoint = { x[0],x[1],0. };
			double azi, inc, cost = 0;
			for (size_t i = 0; i < 4; ++i) {
				inc = CartesianToSpherical(targets[i] - InitialPoint).first;
				azi = CartesianToSpherical(targets[i] - InitialPoint).second;
				if (azi - 110 > EPSILON and azi - 170 < EPSILON and inc - 0 > EPSILON and inc - 40 < EPSILON) {
					cost += (targets[i] - InitialPoint).norm();
					continue;
				}
				/*if (azi - 110 - 180 > EPSILON and azi - 170 - 180 < EPSILON and inc - 0 > EPSILON and inc - 40 < EPSILON) {
					cost += (targets[i] - InitialPoint).norm();
					continue;
				}*/
				cost+= 1/EPSILON/EPSILON;
			}
			return cost;
		};
		PSOvalueType optHead = PSO(headScore, minVs, maxVs, 50, 2,500);
		getOptData(optHead);

	}
	return 0;
}