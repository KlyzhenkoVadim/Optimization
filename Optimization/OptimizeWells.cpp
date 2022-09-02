#include "OptimizeWells.h"


double AzimuthChooser(Constraint cs, double r) {
	// r \in U[0,1]
	if (r < 0.5) {
		return cs.lMin.phi + 2 * r * (cs.lMax.phi - cs.lMin.phi);
	}
	else {
		cs.change();
		return cs.lMin.phi + 2 * (r - 0.5) * (cs.lMax.phi - cs.lMin.phi);
	}
};

double PenaltyLength(double length, double penalty) {
	double answer = length > 5000 ? penalty : 0;
	return answer;
};

double PenaltyDLS(std::vector<TrajectoryTemplate*>& Well, double penalty) {
	double alpha, dotprod, dls;
	int s = solve(Well);
	for (size_t i = 0; i < Well.size() - 1; ++i) {
		Well[i]->getInitPoint(CoordinateSystem::MD);
		Well[i]->getTarget1Point(CoordinateSystem::MD);
		Eigen::Vector3d tprev{ Well[i]->pointInitialMD[1],Well[i]->pointInitialMD[2],Well[i]->pointInitialMD[3] };
		Eigen::Vector3d tnext{ Well[i]->pointMDT1[1],Well[i]->pointMDT1[2],Well[i]->pointMDT1[3] };
		double length = Well[i]->pointMDT1[0] - Well[i]->pointInitialMD[0];
		dotprod = tprev.dot(tnext) / tprev.norm() / tnext.norm();
		if (length < EPSILON)
			continue;
		dotprod = abs(dotprod - 1) < EPSILON ? 1 : dotprod;
		alpha = 180 / PI * acos(dotprod);
		dls = 10 * alpha / length;
		if (dls > 1.5) {
			return penalty;
		}
	}
	return 0;
};

double PenaltyIncTarget(std::vector<TrajectoryTemplate*>& Well, double penalty) {
	int s = solve(Well);
	Well.back()->getTarget1Point(CoordinateSystem::MD);
	Eigen::Vector3d tangent = {Well.back()->pointMDT1[1],Well.back()->pointMDT1[2],Well.back()->pointMDT1[3] };
	std::pair<double, double> incAzi = CartesianToSpherical(tangent);
	if (incAzi.first - 75 > EPSILON) {
		return penalty;
	}
	return 0;
};

double PenaltyAlphaTarget(std::vector<TrajectoryTemplate*>& Well, double penalty) {
	int s = solve(Well);
	Well.back()->getInitPoint(CoordinateSystem::MD);
	Well.back()->getTarget1Point(CoordinateSystem::MD);
	Eigen::Vector3d t1 = { Well.back()->pointInitialMD[1],Well.back()->pointInitialMD[2],Well.back()->pointInitialMD[3] };
	Eigen::Vector3d t2 = { Well.back()->pointMDT1[1] ,Well.back()->pointMDT1[2] ,Well.back()->pointMDT1[3] };
	double dotprod = t1.dot(t2) / t1.norm() / t2.norm();
	dotprod = fabs(dotprod - 1) < EPSILON ? 1 : dotprod;
	double alpha = acos(dotprod);
	if (alpha - PI / 3 > EPSILON) {
		return penalty;
	}
	return penalty;

}

std::vector<TrajectoryTemplate*> Well(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs) {
	// const Eigen::Vector3d & pInit, const Eigen::Vector3d & pTarget
	std::vector<TrajectoryTemplate*> tr;
	tr.push_back(new Hold(pinit, 0, 0, x[0])); // holdLength
	tr.back()->getTarget1Point();
	std::vector<Layer> layers{ {x[0],0,0}, { 900,x[1],AzimuthChooser(cs[0],x[2])},
		{2000,x[3],AzimuthChooser(cs[1],x[4])},{2600,x[5],AzimuthChooser(cs[2],x[6])} };
	double tvd0 = x[0] + (900 - x[0]) * x[7];
	std::vector<double> tvds{ tvd0,x[8],x[9] };
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

Eigen::Vector3d OptimizeWellHead(std::vector<double> minValuesWellHead, std::vector<double> maxValuesWellHead) {
	std::vector<double> inert(500, 0.9);
	double R = 1200 / PI;

	Eigen::Vector3d target1 = { 5803236,682857,2900 }, target2 = { 5803529,682498,2900 }, target3 = { 5803536,683257,2900 }, target4 = { 5803409,683700,2900 };;
	std::vector<Constraint> cs = { { {900,30,110},{900,90,170} },{ {2000,0,200},{2000,50,260} },{ {2600,0,20},{2600,40,80} } };
	std::vector<Eigen::Vector3d> targets{ target1,target2,target3,target4 };

	std::vector<double> minValues, maxValues;
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

	std::function<double(const Eigen::VectorXd&)> scoreWellhead = [&](const Eigen::VectorXd& x) {
		std::vector<std::vector<Eigen::Vector3d>> pCWells;
		std::vector<std::vector<Eigen::Vector4d>> pMDWells;
		double TVDShift;
		Eigen::Vector3d target = target1;
		Eigen::Vector3d pinit = { x[0],x[1],0 };
		double SumCost = 0;
		std::function<double(const Eigen::VectorXd&)> CollRiskScore = [&](const Eigen::VectorXd& y) {
			std::vector<TrajectoryTemplate*> currWell = Well(y, pinit, target, cs);
			double score = orderScore1(currWell,pCWells, pMDWells, TVDShift);
			double length = allLength(currWell), dlsPen = PenaltyDLS(currWell);
			for (auto y : currWell)
				delete y;
			return score + PenaltyLength(length) + dlsPen;
		};
		for (size_t idd = 0; idd < 4; ++idd) {
			PSOvalueType opt = PSO(CollRiskScore, minValues, maxValues, 25, 11, 50);
			//ShowOptData(opt, cs);
			//getOptData(opt);
			SumCost += opt.second;
			std::vector<TrajectoryTemplate*> well = Well(opt.first, pinit, target, cs);
			int cond = solve(well);
			pCWells.push_back(allPointsCartesian(well));
			pMDWells.push_back(allPointsMD(well));
			TVDShift = std::min(TVDShift, opt.first[0]);
			target = targets[idd + 1];
			std::cout << "Length: " << allLength(well) << "\n\n\n";
		}
		return SumCost;
	};

	PSOvalueType optWellhead = PSO(scoreWellhead, minValuesWellHead, maxValuesWellHead, 3, 2,  25);
	return optWellhead.first;
}