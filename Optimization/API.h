#pragma once
#include "TrajectoryTemplate.h"
#include "CurveHold.h"
#include "CurveHoldCurveHold.h"
#include "Curve.h"
#include "Hold.h"
#include "CostFuncs.h"
#include "PSO.h"

struct GeoPoint {
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

struct Point2d
{
	double north;
	double east;
	
};

struct WellPad
{
	bool isWellPad;

	Point2d coord; // координаты центра кустовой площадки

	std::vector<GeoPoint> geoAims; // координаты геологических целей
};



using wellType = std::function<std::vector<TrajectoryTemplate*>(const Eigen::VectorXd& x)>;

class Solver {
private:
	Eigen::Vector3d pointInitial; // Координата устья траектории
	Eigen::Vector3d pointsT1, pointsT3; // Координаты целей траектории.
	std::vector<TrajectoryTemplate*> trajectory;
	int condition = 0;
	std::vector<Eigen::Vector3d> pCtrajectory;
	std::vector<Eigen::Vector4d> pMDtrajectory;
	wellType mainWell;
	WellTrajectoryConstraints OptimizeConstraints{ 400,1.5,5000,2000,2000 };
	PSOvalueType optData;
	bool horizontal = true;

public:
	Solver();
	void setPSOdata();
	void setData(Point2d& pInitial,GeoPoint& Targets);
	void setConstraints(const WellTrajectoryConstraints& cs);
	void Optimize();
	PSOvalueType getPSOdata();
	std::vector<Eigen::Vector3d> getTrajectoryPoints();
	std::vector<Eigen::Vector3d> getInclinometry();
	double getTrajectoryLength();
};