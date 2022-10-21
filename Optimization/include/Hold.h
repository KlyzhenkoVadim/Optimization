#pragma once
#include "TrajectoryTemplate.h"

enum class typeHold {md,TVD,pEnd};
class Hold : public TrajectoryTemplate
{
private:
	Eigen::Vector3d pi; // ƒекартова координата начала Hold'a
	Eigen::Vector3d pf; // ƒекартова координата конца Hold'a
	Eigen::Vector3d direction; // направление Hold'a - единичный вектор 
	double Length, MDTVD;
	typeHold type = typeHold::pEnd;
	size_t nums;
	int condition;
	int fit();

public:
	Hold(const Eigen::Vector3d& pi, const Eigen::Vector3d& pf, size_t nums = 5);
	Hold(const Eigen::Vector3d& pi, double inc, double azi, double L,typeHold type = typeHold::md, size_t nums = 5);
	int getCondition() override;
	void points(CoordinateSystem coordinateSystem) override;
	double length() override;
	double getTortuosity() override;
	void getInitPoint(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
	void getTarget1Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
	void getTarget3Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
	Eigen::Vector3d FunctionPoint(double md) override;
	Eigen::Vector3d FunctionTangent(double md) override;
};
