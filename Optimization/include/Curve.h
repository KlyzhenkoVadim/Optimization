#pragma once
#include "TrajectoryTemplate.h"

enum class TypeCurve { DLS,TVD };

class Curve : public TrajectoryTemplate {
private:
	Eigen::Vector3d pi,pf;
	double inc1, azi1,R, RTVD;
	Eigen::Vector3d t1, t2; // касательные вектора начала и конца дуги.
	double alpha;
	double inc2, azi2;
	size_t nums;
	TypeCurve type;
	int condtition;
	int fit();

public:
	Curve(const Eigen::Vector3d& pi, double inc1, double azi1, double inc2, double azi2, double RTVD, TypeCurve type = TypeCurve::DLS, size_t nums = 20);// nums = 20
	int getCondition() override;
	void points(CoordinateSystem coordinateSystem) override;
	double length() override;
	double getTortuosity() override;
	void getInitPoint(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
	void getTarget1Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
	void getTarget3Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
};
