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
		lMax.phi = (lMax.phi + 180.) - 360 > EPSILON ? lMax.phi - 180 : lMax.phi + 180;
		lMin.phi = (lMin.phi + 180.) - 360 > EPSILON ? lMin.phi - 180 : lMin.phi + 180;
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
/*
bool isSatisfiedConstraint(Constraint constraint, Eigen::Vector3d tangent) {
	std::pair<double,double> thetaPhi = CartesianToSpherical(tangent);
	bool thetaGood, phiGood;
	if (constraint.thetaMin - thetaPhi.first < EPSILON and constraint.thetaMax - thetaPhi.first > EPSILON)
		thetaGood = true;
	else thetaGood = false;
	if (constraint.phiMin - thetaPhi.second < EPSILON and constraint.phiMax - thetaPhi.second > EPSILON)
		phiGood = true;
	else phiGood = false;

	if (thetaPhi.first < EPSILON)
		return thetaGood;
	return thetaGood * phiGood;
}*/
void showVector(Eigen::Vector3d vec) {
	for (auto x : vec)
		std::cout << x << ",";
	std::cout << "\n";
};

double PenaltyLength(double length, double penalty = 100) {
	double answer = length > 5000 ? penalty : 0;
	return answer;
};

double PenaltyDLS(std::vector<Eigen::Vector4d> pMD, double penalty = 100) {
	double alpha, dotprod, dls;
	for (size_t i = 1; i < pMD.size(); ++i) {
		Eigen::Vector3d tprev{ pMD[i - 1][1],pMD[i - 1][2],pMD[i - 1][3] }, tnext{ pMD[i][1],pMD[i][2],pMD[i][3] };
		dotprod = tprev.dot(tnext) / tprev.norm() / tnext.norm();
		dotprod = abs(dotprod - 1) < EPSILON ? 1 : dotprod;
		alpha = 180/PI*acos(dotprod);
		dls = 10 * alpha / (pMD[i][0] - pMD[i][1]);
		if (dls > 1.5) {
			return penalty;
		}
	}
	return 0;
}

int main()
{	
	double R = 1200 / PI;
	Eigen::Vector3d pinit = { 5805626,683999,0 }, target = { 5803236,682857,2900 };
	std::function<std::vector<TrajectoryTemplate*>(const Eigen::VectorXd& x)> Well = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> tr;
		tr.push_back(new Hold(pinit, 0, 0, x[0])); // holdLength
		tr.back()->getTarget1Point();
		std::vector<Layer> layers{ {x[0],0,0}, { 900,x[1],x[2] },{2000,x[3],x[4]},{2600,x[5],x[6]} };
		std::vector<double> tvds{ x[7],x[8],x[9] };
		// tvd (x[0],900)
		// tvd (900,2000)
		// tvd (2000,2600)
		double RCurveHold = x[10];
		//  theta (30,90) phi (110-170)
		//	theta (0,50) phi (200-260)	 выбрал  было (20-80)
		//	theta (0,40) phi (200-260)	выбрал	 было (20-80)
		for (size_t idx = 1; idx < layers.size(); ++idx) {
			tr.push_back(new Curve(tr.back()->pointT1, layers[idx - 1].theta, layers[idx - 1].phi, layers[idx].theta, layers[idx].phi, tvds[idx-1], TypeCurve::TVD));
			tr.back()->getTarget1Point();
			tr.push_back(new Hold(tr.back()->pointT1, layers[idx].theta, layers[idx].phi, layers[idx].TVD,typeHold::TVD));
			tr.back()->getTarget1Point();
		}
		tr.push_back(new CurveHold(tr.back()->pointT1, target, layers.back().theta, layers.back().phi, RCurveHold));
		return tr;
	};
	std::vector<Constraint> cs= { { {900,30,110},{900,90,170} },{ {2000,0,200},{2000,50,260} },{ {2600,0,200},{2600,40,260} } };
	
	std::vector<double> minValues,maxValues;
	minValues.push_back(100);
	maxValues.push_back(400);
	for (size_t i = 0; i < 3; ++i) {
		minValues.push_back(cs[i].lMin.theta);
		minValues.push_back(cs[i].lMin.phi);
		maxValues.push_back(cs[i].lMax.theta);
		maxValues.push_back(cs[i].lMax.phi);
	}
	minValues.push_back(700); // tvd промежуточный в точке x[0] - 900 ????
	maxValues.push_back(900);
	for (size_t j = 1; j < 3; ++j) {
		minValues.push_back(cs[j - 1].lMin.TVD);
		maxValues.push_back(cs[j].lMax.TVD);
	}
	minValues.push_back(R);
	maxValues.push_back(1500);



	std::function<double(const Eigen::VectorXd&)> scoreOne = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> tmp = Well(x);
		double oneScore = OneWellScore(tmp), length = allLength(tmp);
		std::vector<Eigen::Vector4d> pMD = allPointsMD(tmp);
		for (auto x : tmp) {
			delete x;
		}
		return oneScore + PenaltyLength(length) + PenaltyDLS(pMD);
	};
	std::vector<double> inert(500, 0.9);

	/*for (size_t idd = 0; idd < 50; ++idd) {
		PSOvalueType opt = PSO(scoreOne, minValues, maxValues, 50, 11, inert, 0.3, 0.5, 150);
		getOptData(opt);
		std::vector<TrajectoryTemplate*> well = Well(opt.first);
		int cond = solve(well);
		std::cout <<"Length: " << allLength(well) << "\n";
	}*/
	std::vector<TrajectoryTemplate*> tTest = Well(Eigen::VectorXd{ {302.286, 30.7503, 169.889, 49.827, 202.575, 39.7179, 211.694, 764.171, 906.212, 2580.87, 495.589} });
	int s = solve(tTest);
	std::vector<Eigen::Vector3d> pC = allPointsCartesian(tTest);
	writeDataCartesian(pC, "output1.txt");
	return 0;
}