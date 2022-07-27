#pragma once
#include "TrajectoryTemplate.h"
#include "CurveHold.h"
#include "CurveHoldCurveHold.h"
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
	double east;
	double north;
};

struct WellPad
{
	bool isWellPad;

	Point2d coord; // ���������� ������ �������� ��������

	std::vector<GeoPoint> geoAims; // ���������� ������������� �����
};

using wellType = std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd& x)>;

class Solver {
private:
	Eigen::Vector3d pointInitial; // ���������� ����� ����������
	Eigen::Vector3d pointsT1, pointsT3; // ���������� ����� ����������.
	std::vector<TrajectoryTemplate*> trajectory;
	std::vector<Eigen::Vector3d> pCtrajectory;
	std::vector<Eigen::Vector4d> pMDtrajectory;
	wellType mainWell;
	PSOvalueType optData;

public:
	Solver();
	void setPSOdata();
	void setData(Point2d& pInitial,GeoPoint& Targets);
	PSOvalueType getPSOdata();
	double getTrajectoryLength();
	void Optimize();
};