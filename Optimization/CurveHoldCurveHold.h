#pragma once
#include "TrajectoryTemplate.h"
#include "CurveHold.h"

struct PointsCartesian {
	std::vector<Eigen::Vector3d> pointsArc1;
	std::vector<Eigen::Vector3d> pointsHold1;
	std::vector<Eigen::Vector3d> pointsArc2;
	std::vector<Eigen::Vector3d> pointsHold2;
};

struct PointsMD {
	std::vector<Eigen::Vector4d> pointsArc1;
	std::vector<Eigen::Vector4d> pointsHold1;
	std::vector<Eigen::Vector4d> pointsArc2;
	std::vector<Eigen::Vector4d> pointsHold2;
};

class CurveHoldCurveHold : public TrajectoryTemplate
{
private:
	Eigen::Vector3d p1; // координата начальной точки траектории (координаты N,E,V)
	Eigen::Vector3d p4; // координата конечной точки траектории (координаты N,E,V)
	Eigen::Vector3d t1; // единичный касательный вектор начальной точки
	Eigen::Vector3d t4;
	Eigen::Vector3d pInter;

	Eigen::Vector3d r1;
	Eigen::Vector3d r4;

	double tetta1; // значения зенитного угла (в градусах) касательного вектора t1 начальной точке
	double phi1; // значения азимутального угла (в градусах) касательного вектора t1 начальной точке
	double tetta4; // значения зенитного угла (в градусах) касательного вектора t1 конечной точке
	double phi4; // значения азимутального угла (в градусах) касательного вектора t1 конечной точке

	double R1; // радиус кривизны дуги к p1
	double R2; // радиус кривизны дуги к p4
	size_t nums;// число точек траектории

	PointsCartesian pointsCartesian;
	PointsMD pointsMD;

	double betta;
	double eps;

	Eigen::Vector3d t;
	Eigen::Vector3d p1Inter;
	Eigen::Vector3d p4Inter;
	double alpha1;
	double alpha2;
	double holdLength;

public:
	CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1, double phi1, double R1, double R2, const Eigen::Vector3d& p4, 
		double tetta4, double phi4, double betta = 0.0, double eps = 10e-4, size_t nums = 100);

	void fit() override;
	void points(CoordinateSystem coordinateSystem) override;
	double length() override;
};

