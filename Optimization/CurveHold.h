#pragma once
#include "TrajectoryTemplate.h"
class CurveHold : public TrajectoryTemplate
{
private:
	Eigen::Vector3d p1; // координата начальной точки траектории (координаты N,E,V)
	Eigen::Vector3d p3; // координата конечной точки траектории (координаты N,E,V)
	Eigen::Vector3d t1; // единичный касательный вектор начальной точки
	Eigen::Vector3d t2; // единичный касательный вектор конечной точки дуги

	// condition??
	

	double teta; // значения зенитного угла (в градусах) касательного вектора t1 начальной точке
	double phi; // значения азимутального угла (в градусах) касательного вектора t1 начальной точке
	double R; // радиус кривизны траектории (R>0)
	size_t nums;// число точек траектории

	std::vector<Eigen::Vector3d> pointsCartesian;
	std::vector<Eigen::Vector4d> pointsMD;
	double alpha;
	double betta;

public:
	CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3, double teta, double phi, double R, size_t nums = 100);
	void fit() override;
	void points(CoordinateSystem coordinateSystem) override;
	double length() override;
	double getAlpha();
	
	Eigen::Vector3d getInitPoint();
	Eigen::Vector3d getTargetPoint();
};

