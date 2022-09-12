#include "API.h"

Solver::Solver() {};

void Solver::setPSOdata() {};

void Solver::setData(Point2d& pInitial, GeoPoint& Targets) {
	pointInitial = { pInitial.north,pInitial.east,0. };
	pointsT1 = { Targets.northT1, Targets.eastT1, Targets.tvdT1 };
	pointsT3 = { Targets.northT3, Targets.eastT3, Targets.tvdT3 };
	if ((pointsT1 - pointsT3).norm() < EPSILON) 
	{
			horizontal = false;
	}
	if (horizontal) {
		mainWell = [&](const Eigen::VectorXd& x) {
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

void Solver::setConstraints(const WellTrajectoryConstraints& cs)
{
	OptimizeConstraints = cs;
}

void Solver::Optimize() {
	size_t numIterations = 150;
	
	std::function<double(Eigen::VectorXd)> score = [&](Eigen::VectorXd x) {
		std::vector<TrajectoryTemplate*> tmpwell = mainWell(x);
		double oneScore = scoreSolver(tmpwell,OptimizeConstraints);
		for (auto x : tmpwell) {
			delete x;
		}
		return oneScore;
	};
	if (horizontal) {
		std::vector<double> minValues{ pointInitial[0] - 1000.,pointInitial[1] - 1000.,100.,0.,-180.,100. };
		std::vector<double> maxValues{ pointInitial[0] + 1000.,pointInitial[1] + 1000.,pointsT1[2],90.,180.,1000. };
		optData = PSO(score, minValues, maxValues, 30, 6, numIterations);
	}
	else {
		std::vector<double> minValues{OptimizeConstraints.minDepthFirstHold,0};
		std::vector<double> maxValues{pointsT1[2],OptimizeConstraints.maxDLS};
		optData = PSO(score, minValues, maxValues, 5*minValues.size(), minValues.size(), numIterations);
	}
	trajectory = mainWell(optData.first);
	condition = solve(trajectory);
};

PSOvalueType Solver::getPSOdata() {
	return optData;
};

double Solver::getTrajectoryLength() {
	if (condition != 0)
	{
		std::cout << "Warning! Impossible to build Trajectory with:\n HoldDepth:" << optData.first[0] << "\nDLS: " << optData.first[1] << "\n";
		return 1 / EPSILON;
	}
	return allLength(trajectory);
};

std::vector<Eigen::Vector3d> Solver::getTrajectoryPoints() {
	// ??? if solve > 0 ???
	if (condition != 0)
	{
		std::cout << "Warning! Impossible to build Trajectory with:\n HoldDepth:" << optData.first[0] << "\nDLS: " << optData.first[1] << "\n";
		pCtrajectory.push_back(pointInitial);
	}
	else
	{
		pCtrajectory = allPointsCartesian(trajectory);
	}
	return pCtrajectory;
};

std::vector<Eigen::Vector3d> Solver::getInclinometry()
{
	std::vector<Eigen::Vector3d> inclinometry;
	if (condition != 0)
	{
		std::cout << "Warning! Impossible to build Trajectory with:\n HoldDepth:" << optData.first[0] << "\nDLS: " << optData.first[1] << "\n";
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



