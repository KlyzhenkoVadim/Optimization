#include "API.h"

Solver::Solver() {};

void Solver::setPSOdata() {};

void Solver::setData(Point2d& pInitial, GeoPoint& Targets) {
	pointInitial = { pInitial.north,pInitial.east,0. };
	pointsT1 = { Targets.northT1, Targets.eastT1, Targets.tvdT1 };
	pointsT3 = { Targets.northT3, Targets.eastT3, Targets.tvdT3 };
	mainWell = [&](Eigen::VectorXd& x) {
		double m = 1800 / PI;
		std::vector<TrajectoryTemplate*> well;
		Eigen::Vector3d pTHold = { pointInitial[0],pointInitial[1],x[5] };
		Eigen::Vector3d pT3Chch1 = { x[0],x[1],x[2] };
		Eigen::Vector3d Tangent = calcTangentVector(x[4], x[3]);
		Eigen::Vector3d pT1Chch1 = pT3Chch1 - 110. * Tangent;
		well.push_back(new Hold(pointInitial, pTHold));
		well.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pT1Chch1, pT3Chch1));
		well.push_back(new CurveHoldCurveHold(pT3Chch1, x[3], x[4], m, m, pointsT1, pointsT3));
		return well;
	};
}

PSOvalueType Solver::getPSOdata() {
	return optData;
};

void Solver::Optimize() {
	size_t numIterations = 150;
	std::vector<double> inertia(numIterations, 0.9);
	std::vector<double> minValues{ pointInitial[0] - 1000.,pointInitial[1] - 1000.,pointsT1[2]-200.,60.,-180.,100.};
	std::vector<double> maxValues{ pointInitial[0] + 1000.,pointInitial[1] + 1000.,pointsT1[2]-50,70.,180.,1000. };
	std::function<double(Eigen::VectorXd)> score = [&](Eigen::VectorXd x) {
		std::vector<TrajectoryTemplate*> tmpwell = mainWell(x);
		return OneWellScore(tmpwell);
	};
	optData = PSO(score, minValues, maxValues, 30, 6, inertia,0.3,0.5,numIterations);
	trajectory = mainWell(optData.first);
};

double Solver::getTrajectoryLength() {
	int condition = solve(trajectory);
	if (condition != 0)
		return -1;
	return allLength(trajectory);
};



