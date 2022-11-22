#include "WellTrajectorySolver.h"

double well_trajectory::scoreSolver(std::vector<TrajectoryTemplate*>& well, const WellTrajectoryConstraints& cs, double penalty)
{
	int condition = solve(well);
	if (condition != 0)
	{
		return penalty * (-condition) / well.size(); // penalty * percent of incorrect templates.
	}
	double md = allLength(well);
	if (md > cs.maxMD)
	{
		return 100.;
	}
	well[0]->getInitPoint();
	well.back()->getTarget1Point();
	well.back()->getTarget3Point();
	double IdealLength = (well.back()->pointT1 - well[0]->pointInitial).norm() + (well.back()->pointT3 - well.back()->pointT1).norm();
	if (IdealLength == 0)
		return 0.;

	double tortuosity = 0;
	for (size_t i = 0; i < well.size(); ++i)
	{
		tortuosity += well[i]->getTortuosity();
	}
	tortuosity = 180. / PI * tortuosity;
	double ahd = 0;
	std::vector<Eigen::Vector3d> points3d = allPointsCartesian(well);
	double NS = 0, EW = 0;
	for (size_t i = 1; i < points3d.size(); ++i)
	{
		NS += fabs(points3d[i][0] - points3d[i - 1][0]);
		EW += fabs(points3d[i][1] - points3d[i - 1][1]);
		ahd += sqrt((points3d[i][0] - points3d[i - 1][0]) * (points3d[i][0] - points3d[i - 1][0])
			+ (points3d[i][1] - points3d[i - 1][1]) * (points3d[i][1] - points3d[i - 1][1]));
	}
	if (EW - cs.maxDistEastWest > EPSILON || NS - cs.maxDistNorthSouth > EPSILON)
	{
		return 100.;
	}
	
	double tvd = points3d.back()[2] - points3d[0][2];
	double ddi = log10((1. / 0.305) * ahd * md * tortuosity / tvd);
	ddi  = 5. / (1 + exp(-2.5 * (ddi - 6.25)));
	return md / IdealLength + ddi;
}

void well_trajectory::WellTrajectorySolver::setPSOdata() {};

void well_trajectory::WellTrajectorySolver::setData(Point2d& pInitial, Target& Targets) {
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

void well_trajectory::WellTrajectorySolver::setConstraints(const WellTrajectoryConstraints& cs)
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

void well_trajectory::WellTrajectorySolver::optimize() {
	size_t numIterations = 150;
	std::vector<double> minValues{ OptimizeConstraints.minDepthFirstHold,0 };
	std::vector<double> maxValues{ pointsT1[2],OptimizeConstraints.maxDLS };
	auto score = [&](const Eigen::VectorXd& x) {
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

PSOvalueType well_trajectory::WellTrajectorySolver::getPSOdata() {
	return optData;
};

double well_trajectory::WellTrajectorySolver::getTrajectoryLength() {
	if (condition != 0 or optData.second > 99)
	{
		std::cout << "Warning! Impossible to build the Trajectory or satisfy the Constraints with:\n HoldDepth:" << optData.first[0] << "\nDLS: " << optData.first[1] << "\n";
		return 1 / EPSILON;
	}
	return allLength(trajectory);
};

std::vector<Eigen::Vector3d> well_trajectory::WellTrajectorySolver::getTrajectoryPoints() {
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

std::vector<Eigen::Vector3d> well_trajectory::WellTrajectorySolver::getInclinometry()
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
