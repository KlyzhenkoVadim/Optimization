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
			return { theta,360 + phi * 180 / PI };
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

std::vector<TrajectoryTemplate*> Well(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs) { 
	// const Eigen::Vector3d & pInit, const Eigen::Vector3d & pTarget
	std::vector<TrajectoryTemplate*> tr;
	tr.push_back(new Hold(pinit, 0, 0, x[0])); // holdLength
	tr.back()->getTarget1Point();
	std::vector<Layer> layers{ {x[0],0,0}, { 900,x[1],AzimuthChooser(cs[0],x[2])},
		{2000,x[3],AzimuthChooser(cs[1],x[4])},{2600,x[5],AzimuthChooser(cs[2],x[6])} };
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


void writeInclinometry(const std::vector<Eigen::Vector4d>& pMD, std::string filename) {
	std::ofstream output;
	output.open(filename);
	output << "md, inc, azi\n";
	for (size_t i = 0; i < pMD.size(); ++i) {
		std::pair<double, double> incAzi = CartesianToSpherical(Eigen::Vector3d{ pMD[i][1],pMD[i][2],pMD[i][3] });
		output << pMD[i][0] << "," << incAzi.first << "," << incAzi.second << std::endl;
	}
	output.close();
}

int main()
{
	double R = 1200 / PI;
	Eigen::Vector3d pinit1 = { 5805626 - 500,683999 - 500,0 }, target1 = { 5803236,682857,2900 };
	std::vector<Constraint> cs = { { {900,30,110},{900,90,170} },
		{ {2000,0,200},{2000,50,260} },{ {2600,0,20},{2600,40,80} } };
	Eigen::Vector3d pinit2 = { 5805626 - 500+15, 683999 - 500, 0 }, target2 = { 5803529,682498,2900 };
	Eigen::Vector3d pinit3 = { 5805626 - 500 - 15,683999 - 500-15,0 }, target3 = { 5803536,683257,2900 };
	Eigen::Vector3d pinit4 = { 5805626 - 500,683999 - 500 +15,0 }, target4 = { 5803409,683700,2900 };
	// 4001Init (NS,EW,TVD) : (5805626-500+15,683999-500,0)
	//4001Target:	(NS,EW,TVD) : (5803529,682498,2900)
	//Target4003:	5803536	683257	2900
	// Target4000	5803409	683700	2900


	//std::function<std::vector<TrajectoryTemplate*>(const Eigen::VectorXd& x)> Well = [&](const Eigen::VectorXd& x) { // const Eigen::Vector3d & pInit, const Eigen::Vector3d & pTarget
	//	std::vector<TrajectoryTemplate*> tr; 
	//	tr.push_back(new Hold(pinit, 0, 0, x[0])); // holdLength
	//	tr.back()->getTarget1Point();
	//	std::vector<Layer> layers{ {x[0],0,0}, { 900,x[1],AzimuthChooser(cs[0],x[2])},
	//		{2000,x[3],AzimuthChooser(cs[1],x[4])},{2600,x[5],AzimuthChooser(cs[2],x[6])}};
	//	std::vector<double> tvds{x[7],x[8],x[9]};
	//	// tvd (x[0],900)
	//	// tvd (900,2000)
	//	// tvd (2000,2600)
	//	double RCurveHold = x[10];
	//	//  theta (30,90) phi (110-170)
	//	//	theta (0,50) phi (200-260)	 выбрал  было (20-80)
	//	//	theta (0,40) phi (200-260)	выбрал	 было (20-80)
	//	for (size_t idx = 1; idx < layers.size(); ++idx) {
	//		tr.push_back(new Curve(tr.back()->pointT1, layers[idx - 1].theta, layers[idx - 1].phi, layers[idx].theta, layers[idx].phi, tvds[idx - 1], TypeCurve::TVD));
	//		tr.back()->getTarget1Point();
	//		tr.push_back(new Hold(tr.back()->pointT1, layers[idx].theta, layers[idx].phi, layers[idx].TVD, typeHold::TVD));
	//		tr.back()->getTarget1Point();
	//	}
	//	tr.push_back(new CurveHold(tr.back()->pointT1, target, layers.back().theta, layers.back().phi, RCurveHold));
	//	return tr;
	//};
	
	std::vector<double> minValues,maxValues;
	minValues.push_back(100);
	maxValues.push_back(400);
	for (size_t i = 0; i < 3; ++i) {
		minValues.push_back(cs[i].lMin.theta);
		minValues.push_back(0); //minValues.push_back(cs[i].lMin.phi);
		maxValues.push_back(cs[i].lMax.theta); //maxValues.push_back(cs[i].lMax.phi);
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
		std::vector<TrajectoryTemplate*> tmp = Well(x,pinit1,target1,cs);
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
	Eigen::Vector3d pinit = pinit1, target = target1;
	std::vector<Eigen::Vector3d> inits = { pinit1,pinit2,pinit3,pinit4 }, targets = { target1,target2,target3,target4 };
	std::function<double(const Eigen::VectorXd&)> CollRiskScore = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> currWell = Well(x, pinit, target, cs);
		double score = orderScore1(currWell, pCWells, pMDWells);
		double length = allLength(currWell), dlsPen = PenaltyDLS(currWell);
		for (auto x : currWell)
			delete x;
		return score + PenaltyLength(length) + dlsPen;
	};

	std::vector<double> inert(500, 0.9);
	bool StartOptimization;
	std::cout << "Start Optimization?";
	std::cin >> StartOptimization;
	if (StartOptimization) {
		for (size_t idd = 0; idd < 4; ++idd) {
			PSOvalueType opt = PSO(CollRiskScore, minValues, maxValues, 50, 11, inert, 0.3, 0.5, 100);
			ShowOptData(opt, cs);
			getOptData(opt);
			std::vector<TrajectoryTemplate*> well = Well(opt.first,pinit,target,cs);
			int cond = solve(well);
			pCWells.push_back(allPointsCartesian(well));
			pMDWells.push_back(allPointsMD(well));
			pinit = inits[idd+1];
			target = targets[idd+1];
			std::cout << "Length: " << allLength(well) << "\n\n\n";
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
	tTest = Well(arg,pinit2,target2,cs);
	//Eigen::Vector3d{0,0,0}
	
	
	PSOvalueType optV = {arg, 0};
	optV.second = scoreOne(optV.first);
	int s = solve(tTest);
	if (!StartOptimization) {
		std::vector<Eigen::Vector3d> pC = allPointsCartesian(tTest);
		ShowOptData(optV, cs);
		std::cout << "Length: " << allLength(tTest);
		writeDataCartesian(pC, "outputCURVE.txt");
		std::vector<Eigen::Vector4d> pMD1 = allPointsMD(tTest);
		writeInclinometry(pMD1, "inclinometry.txt");
	}
	return 0;
}



/*Start Optimization? FOR 3 WELLS
Optimization result:
_____________________________________
CostValue: 1.04008
Hold0TVD: 348.832
Curve1 (inc1, azi1): 33.9519, 169.986
Hold1TVDstart: 760.08
Curve2 (inc2, azi2): 44.8093, 204.357
Hold2TVDstart: 1281.01
Curve3 (inc3, azi3): 39.7848, 206.152
Hold3TVDstart: 2244.49
CurveHoldDLS: 0.779893
gBestCost: 1.04008
gBestPos: 348.832, 33.9519, 0.499887, 44.8093, 0.0363084, 39.7848, 0.551264, 760.08, 1281.01, 2244.49, 734.663,
Length: 3661.65


Optimization result:
_____________________________________
CostValue: 1.08319
Hold0TVD: 183.391
Curve1 (inc1, azi1): 30.7394, 290.001
Hold1TVDstart: 884.946
Curve2 (inc2, azi2): 49.0219, 201.502
Hold2TVDstart: 1374.46
Curve3 (inc3, azi3): 39.4466, 204.633
Hold3TVDstart: 2088.23
CurveHoldDLS: 1.45948
gBestCost: 1.08319
gBestPos: 183.391, 30.7394, 0.500008, 49.0219, 0.0125135, 39.4466, 0.538611, 884.946, 1374.46, 2088.23, 392.576,
Length: 3753.75


Optimization result:
_____________________________________
CostValue: 1.01598
Hold0TVD: 140.255
Curve1 (inc1, azi1): 33.8687, 169.752
Hold1TVDstart: 400.9
Curve2 (inc2, azi2): 32.5051, 200.034
Hold2TVDstart: 1279.28
Curve3 (inc3, azi3): 29.4436, 200.006
Hold3TVDstart: 2212.82
CurveHoldDLS: 0.896738
gBestCost: 1.01598
gBestPos: 140.255, 33.8687, 0.497935, 32.5051, 0.000285118, 29.4436, 0.500051, 400.9, 1279.28, 2212.82, 638.935,
Length: 3361.85*/



/*Start Optimization? 4 WELLS
Optimization result:
_____________________________________
CostValue: 1.04218
Hold0TVD: 336.151
Curve1 (inc1, azi1): 58.5904, 170
Hold1TVDstart: 899.719
Curve2 (inc2, azi2): 40.9049, 201.209
Hold2TVDstart: 1027.85
Curve3 (inc3, azi3): 39.6638, 215.827
Hold3TVDstart: 2254.31
CurveHoldDLS: 1.15657
gBestCost: 1.04218
gBestPos: 336.151, 58.5904, 0.499997, 40.9049, 0.0100758, 39.6638, 0.631888, 899.719, 1027.85, 2254.31, 495.393,
Length: 3669.05


Optimization result:
_____________________________________
CostValue: 1.09659
Hold0TVD: 110.961
Curve1 (inc1, azi1): 30.0207, 290.156
Hold1TVDstart: 857.542
Curve2 (inc2, azi2): 35.4675, 206.155
Hold2TVDstart: 1151
Curve3 (inc3, azi3): 39.9698, 200.522
Hold3TVDstart: 2396.91
CurveHoldDLS: 1.41415
gBestCost: 1.09659
gBestPos: 110.961, 30.0207, 0.501299, 35.4675, 0.0512904, 39.9698, 0.504352, 857.542, 1151, 2396.91, 405.16,
Length: 3800.37


Optimization result:
_____________________________________
CostValue: 1.02549
Hold0TVD: 113.302
Curve1 (inc1, azi1): 47.3368, 167.221
Hold1TVDstart: 795.645
Curve2 (inc2, azi2): 33.022, 200.017
Hold2TVDstart: 1163.06
Curve3 (inc3, azi3): 22.9157, 200.079
Hold3TVDstart: 2253.71
CurveHoldDLS: 1.14772
gBestCost: 1.02549
gBestPos: 113.302, 47.3368, 0.476839, 33.022, 0.000141088, 22.9157, 0.500662, 795.645, 1163.06, 2253.71, 499.213,
Length: 3392.86


Optimization result:
_____________________________________
CostValue: 1.52425
Hold0TVD: 331.609
Curve1 (inc1, azi1): 75.5996, 131.282
Hold1TVDstart: 871.431
Curve2 (inc2, azi2): 34.5653, 216.639
Hold2TVDstart: 1490.88
Curve3 (inc3, azi3): 23.7186, 207.123
Hold3TVDstart: 2512.73
CurveHoldDLS: 0.98891
gBestCost: 1.52425
gBestPos: 331.609, 75.5996, 0.177349, 34.5653, 0.138658, 23.7186, 0.559356, 871.431, 1490.88, 2512.73, 579.383,
Length: 3772.22*/