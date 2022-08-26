#pragma once
#include "TrajectoryTemplate.h"

enum class TypeCurve { DLS,TVD };

class Curve : public TrajectoryTemplate {
private:
	Eigen::Vector3d pi,pf;
	double inc1, azi1,R;
	Eigen::Vector3d t1, t2; // ����������� ������� ������ � ����� ����.
	double alpha;
	double inc2, azi2;
	size_t nums;

public:
	Curve(const Eigen::Vector3d& pi, double inc1, double azi1, double inc2, double azi2, double RTVD, TypeCurve type = TypeCurve::DLS, size_t nums = 15);
	void fit() override;
	void points(CoordinateSystem coordinateSystem) override;
	double length() override;
	void getInitPoint() override;
	void getTarget1Point() override;
	void getTarget3Point() override;
};
