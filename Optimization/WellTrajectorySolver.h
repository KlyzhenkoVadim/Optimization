#pragma once
#include "TrajectoryTemplate.h"
#include "CurveHold.h"
#include "CurveHoldCurveHold.h"
#include "Curve.h"
#include "Hold.h"
#include "CostFuncs.h"
#include "PSO.h"
#include <string>

namespace well_trajectory
{
	struct Point2d
	{
		double north;
		double east;
	};

	struct Target {
		double northT1;
		double eastT1;
		double tvdT1;

		double northT3;
		double eastT3;
		double tvdT3;

		double H;

		double frGas;
		double frLiq;

		std::string name;
		std::string date;
	};

	using wellType = std::function<std::vector<TrajectoryTemplate*>(const Eigen::VectorXd& x)>;

	class WellTrajectorySolver {
	private:
		Eigen::Vector3d pointInitial; // Координата устья траектории
		Eigen::Vector3d pointsT1, pointsT3; // Координаты целей траектории.
		std::vector<TrajectoryTemplate*> trajectory;
		int condition = 0;
		std::vector<Eigen::Vector3d> pCtrajectory;
		std::vector<Eigen::Vector4d> pMDtrajectory;
		wellType mainWell;
		WellTrajectoryConstraints OptimizeConstraints{ 400,1.5,1 / EPSILON,1 / EPSILON,1 / EPSILON };
		PSOvalueType optData;
		bool horizontal = true;

	public:
		WellTrajectorySolver();
		void setPSOdata();
		void setData(Point2d& pInitial, Target& Targets);
		void setConstraints(const WellTrajectoryConstraints& cs);
		void optimize();
		PSOvalueType getPSOdata();
		std::vector<Eigen::Vector3d> getTrajectoryPoints();
		std::vector<Eigen::Vector3d> getInclinometry();
		double getTrajectoryLength();
	};
}
