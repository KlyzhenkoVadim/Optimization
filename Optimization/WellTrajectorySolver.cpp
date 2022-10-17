#include "WellTrajectorySolver.h"

using namespace well_trajectory;

WellTrajectorySolver::WellTrajectorySolver() {};

void WellTrajectorySolver::setPSOdata() {};

void WellTrajectorySolver::setData(Point2d& pInitial, Target& Targets) {
	pointInitial = { pInitial.north,pInitial.east,0. };
	pointsT1 = { Targets.northT1, Targets.eastT1, Targets.tvdT1 };
	pointsT3 = { Targets.northT3, Targets.eastT3, Targets.tvdT3 };
	if ((pointsT1 - pointsT3).norm() < EPSILON)
	{
		horizontal = false;
	}
	if (horizontal) {
		mainWell = [&](const Eigen::VectorXd& x) {
			double R1 = x[1] < EPSILON ? 1 / EPSILON : 1800 / x[1] / PI;
			double R2 = x[2] < EPSILON ? 1 / EPSILON : 1800 / x[2] / PI;
			std::vector<TrajectoryTemplate*> well;
			well.push_back(new Hold(pointInitial, 0, 0, x[0], typeHold::TVD));
			well.back()->getTarget1Point();
			well.push_back(new CurveHoldCurveHold(well.back()->pointT1, 0, 0, R1, R2, pointsT1, pointsT3));
			return well;
		};
	}
	else {
		mainWell = [&](const Eigen::VectorXd& x)
		{
			double R = x[1] < EPSILON ? 1 / EPSILON : 1800 / PI / x[1];
			std::vector<TrajectoryTemplate*> well;
			well.push_back(new Hold(pointInitial, 0, 0, x[0], typeHold::TVD));
			well.back()->getTarget1Point();
			well.push_back(new CurveHold(well.back()->pointT1, pointsT1, 0, 0, R));
			return well;
		};
	}
}

void WellTrajectorySolver::setConstraints(const WellTrajectoryConstraints& cs)
{
	if (!std::isnan(cs.minDepthFirstHold))
		OptimizeConstraints.minDepthFirstHold = cs.minDepthFirstHold;
	if (!std::isnan(cs.maxMD))
		OptimizeConstraints.maxMD = cs.maxMD;
	if (!std::isnan(cs.maxDLS))
		OptimizeConstraints.maxDLS = cs.maxDLS;
	if (!std::isnan(cs.maxDistNorthSouth))
		OptimizeConstraints.maxDistNorthSouth = cs.maxDistNorthSouth;
	if (!std::isnan(cs.maxDistEastWest))
		OptimizeConstraints.maxDistEastWest = cs.maxDistEastWest;
}

void WellTrajectorySolver::optimize() {
	size_t numIterations = 150;
	std::vector<double> minValues{ OptimizeConstraints.minDepthFirstHold,0 };
	std::vector<double> maxValues{ pointsT1[2],OptimizeConstraints.maxDLS };
	std::function<double(Eigen::VectorXd)> score = [&](Eigen::VectorXd x) {
		std::vector<TrajectoryTemplate*> tmpwell = mainWell(x);
		double oneScore = scoreSolver(tmpwell, OptimizeConstraints);
		for (auto x : tmpwell) {
			delete x;
		}
		return oneScore;
	};
	if (horizontal) {
		minValues.push_back(0);
		maxValues.push_back(OptimizeConstraints.maxDLS);
	}
	optData = PSO(score, minValues, maxValues, 5 * minValues.size(), minValues.size(), numIterations);
	trajectory = mainWell(optData.first);
	condition = solve(trajectory);
};

PSOvalueType WellTrajectorySolver::getPSOdata() {
	return optData;
};

double WellTrajectorySolver::getTrajectoryLength() {
	if (condition != 0 or optData.second > 99)
	{
		std::cout << "Warning! Impossible to build the Trajectory or satisfy the Constraints with:\n HoldDepth:" << optData.first[0] << "\nDLS: " << optData.first[1] << "\n";
		return 1 / EPSILON;
	}
	return allLength(trajectory);
};

std::vector<Eigen::Vector3d> WellTrajectorySolver::getTrajectoryPoints() {
	// ??? if solve > 0 ???
	if (condition != 0 or optData.second > 99)
	{
		std::cout << "Warning! Impossible to build the Trajectory or satisfy the Constraints with:\n HoldDepth:" << optData.first[0] << "\nDLS: " << optData.first[1] << "\n";
		pCtrajectory.push_back(pointInitial);
	}
	else
	{
		pCtrajectory = allPointsCartesian(trajectory);
	}
	return pCtrajectory;
};

std::vector<Eigen::Vector3d> WellTrajectorySolver::getInclinometry()
{
	std::vector<Eigen::Vector3d> inclinometry;
	if (condition != 0 or optData.second > 99)
	{
		std::cout << "Warning! Impossible to build the Trajectory or satisfy the Constraints with:\n HoldDepth:" << optData.first[0] << "\nDLS: " << optData.first[1] << "\n";
		trajectory[0]->getInitPoint(CoordinateSystem::MD);
		pMDtrajectory.push_back(trajectory[0]->pointInitialMD);
	}
	pMDtrajectory = allPointsMD(trajectory);
	std::pair<double, double> tmpIncAzi;
	for (size_t i = 0; i < pMDtrajectory.size(); ++i)
	{
		tmpIncAzi = CartesianToSpherical(Eigen::Vector3d{ pMDtrajectory[i][1],pMDtrajectory[i][2],pMDtrajectory[i][3] });
		inclinometry.push_back({ pMDtrajectory[i][0],tmpIncAzi.first,tmpIncAzi.second });
	}
	return inclinometry;
}
