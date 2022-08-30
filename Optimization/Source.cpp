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

std::pair<double, double> CartesianToSpherical(Eigen::Vector3d t) {
	t.normalize();
	double theta = acos(t[2]), phi;
	if (theta < EPSILON)
		return { 0,0 };
	theta = theta * 180 / PI;
	phi = atan(t[1] / t[0]);
	if (phi > 0) {
		if (t[0] >= 0 and t[1] >= 0)
			return { theta,phi * 180 / PI };
		else if (t[0] < 0 and t[1] < 0) {
			phi += PI;
			return { theta,180 / PI * phi };
		}
	}
	else {
		if (t[0] >= 0 and t[1] <= 0)
			return { theta,phi * 180 / PI };
		else if (t[0] < 0 and t[1] > 0) {
			phi += PI;
			return{ theta,180 / PI * phi };
		}
		}
}

void showVector(Eigen::Vector3d vec) {
	for (auto x : vec)
		std::cout << x << ",";
	std::cout << "\n";
};

double PenaltyLength(double length, double penalty = 100) {
	double answer = length > 5000 ? penalty : 0;
	return answer;
};

double PenaltyDLS(std::vector<TrajectoryTemplate*> &Well, double penalty = 100) {
	double alpha, dotprod, dls;
	int s = solve(Well);
	for (size_t i = 0; i < Well.size()-1; ++i) {
		Well[i]->getInitPoint(CoordinateSystem::MD);
		Well[i]->getTarget1Point(CoordinateSystem::MD);
		Eigen::Vector3d tprev{ Well[i]->pointInitialMD[1],Well[i]->pointInitialMD[2],Well[i]->pointInitialMD[3] };
		Eigen::Vector3d tnext{ Well[i]->pointMDT1[1],Well[i]->pointMDT1[2],Well[i]->pointMDT1[3] };
		double length = Well[i]->pointMDT1[0] - Well[i]->pointInitialMD[0];
		dotprod = tprev.dot(tnext) / tprev.norm() / tnext.norm();
		if (length < EPSILON)
			continue;
		dotprod = abs(dotprod - 1) < EPSILON ? 1 : dotprod;
		alpha = 180/PI*acos(dotprod);
		dls = 10 * alpha / length;
		if (dls > 1.5) {
			return penalty;
		}
	}
	return 0;
}

//double PenaltyAHD(std::vector<TrajectoryTemplate*>& Well, double penalty = 100) {
//	std::vector<Eigen::Vector3d> pC = allPointsCartesian(Well);
//	double ahd = 0;
//	for (size_t i = 1; i < pC.size(); ++i) {
//		Eigen::Vector3d tmp = pC[i] - pC[i - 1];
//		tmp[2] = 0;
//		ahd += tmp.norm();
//	}
//	if (ahd - 2000 > EPSILON)
//		return penalty;
//	return 0;
//}

double AzimuthChooser(Constraint cs, double r) {
	// r \in U[0,1]
	if (r < 0.5) {
		return cs.lMin.phi + 2 * r * (cs.lMax.phi - cs.lMin.phi);
	}
	else {
		cs.change();
		return cs.lMin.phi + 2 * (r - 0.5) * (cs.lMax.phi - cs.lMin.phi);
	}
}

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

int main()
{	
	double R = 1200 / PI;
	Eigen::Vector3d pinit = { 5805626,683999,0 }, target = { 5803236,682857,2900 };
	std::vector<Constraint> cs = { { {900,30,110},{900,90,170} },
		{ {2000,0,200},{2000,50,260} },{ {2600,0,20},{2600,40,80} } };

	std::function<std::vector<TrajectoryTemplate*>(const Eigen::VectorXd& x)> Well = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> tr;
		tr.push_back(new Hold(pinit, 0, 0, x[0])); // holdLength
		tr.back()->getTarget1Point() ;
		std::vector<Layer> layers{ {x[0],0,0}, { 900,x[1],AzimuthChooser(cs[0],x[2])},
			{2000,x[3],AzimuthChooser(cs[1],x[4])},{2600,x[5],AzimuthChooser(cs[2],x[6])}};
		std::vector<double> tvds{ x[7],x[8],x[9] };
		// tvd (x[0],900)
		// tvd (900,2000)
		// tvd (2000,2600)
		double RCurveHold = x[10];
		//  theta (30,90) phi (110-170)
		//	theta (0,50) phi (200-260)	 выбрал  было (20-80)
		//	theta (0,40) phi (200-260)	выбрал	 было (20-80)
		for (size_t idx = 1; idx < layers.size(); ++idx) {
			tr.push_back(new Curve(tr.back()->pointT1, layers[idx - 1].theta, layers[idx - 1].phi, layers[idx].theta, layers[idx].phi, tvds[idx - 1], TypeCurve::TVD));
			tr.back()->getTarget1Point();
			tr.push_back(new Hold(tr.back()->pointT1, layers[idx].theta, layers[idx].phi, layers[idx].TVD, typeHold::TVD));
			tr.back()->getTarget1Point();
		}
		tr.push_back(new CurveHold(tr.back()->pointT1, target, layers.back().theta, layers.back().phi, RCurveHold));
		return tr;
	};
	
	std::vector<double> minValues,maxValues;
	minValues.push_back(100);
	maxValues.push_back(400);
	for (size_t i = 0; i < 3; ++i) {
		minValues.push_back(cs[i].lMin.theta);
		//minValues.push_back(cs[i].lMin.phi);
		minValues.push_back(0);
		maxValues.push_back(cs[i].lMax.theta);
		//maxValues.push_back(cs[i].lMax.phi);
		maxValues.push_back(1);
	}
	minValues.push_back(400); // tvd промежуточный в точке x[0] - 900 ????
	maxValues.push_back(900);
	for (size_t j = 1; j < 3; ++j) {
		minValues.push_back(cs[j - 1].lMin.TVD);
		maxValues.push_back(cs[j].lMax.TVD);
	}
	minValues.push_back(R);
	maxValues.push_back(1500);



	std::function<double(const Eigen::VectorXd&)> scoreOne = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> tmp = Well(x);
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
	std::vector<double> inert(500, 0.9);
	bool StartOptimization;
	std::cout << "Start Optimization?";
	std::cin >> StartOptimization;
	if (StartOptimization) {
		for (size_t idd = 0; idd < 15; ++idd) {
			PSOvalueType opt = PSO(scoreOne, minValues, maxValues, 50, 11, inert, 0.3, 0.5, 100);
			ShowOptData(opt, cs);
			getOptData(opt);
			std::vector<TrajectoryTemplate*> well = Well(opt.first);
			int cond = solve(well);
			std::cout << "Length: " << allLength(well) << "\n\n\n";
		}
	}
	// Hold1 inc1 azi1 inc2 azi2 inc3 azi3 tvd1 tvd2 tvd3 R
	std::vector<TrajectoryTemplate*> tTest;
	
	Eigen::VectorXd arg{ { 109.931, 39.0791, 0.495777, 49.8637, 0.188651, 39.9068, 0.552619, 427.504, 1541.77, 2588.57, 382.529} };
	tTest = Well(arg);
	PSOvalueType optV = {arg, 0};
	optV.second = scoreOne(optV.first);
	int s = solve(tTest);
	if (!StartOptimization) {
		/*std::vector<Eigen::Vector3d> pC = allPointsCartesian(tTest);
		ShowOptData(optV, cs);
		std::cout << "Length: " << allLength(tTest);
		writeDataCartesian(pC, "output10.txt");*/
		Eigen::Vector3d pii{ 0,0,0 }, ptt{ 0,0,1000 };
		std::vector<TrajectoryTemplate*> t11{ new CurveHoldCurveHold(pii,0,0,700,700,ptt,90,0) };
		t11[0]->fit();
		std::vector<Eigen::Vector3d> pcc = allPointsCartesian(t11);
		writeDataCartesian(pcc, "outputkust.txt");
	
	}
	return 0;
}