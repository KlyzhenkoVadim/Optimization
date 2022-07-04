#pragma once
#include "TrajectoryTemplate.h"

class Hold : public TrajectoryTemplate
{
private:
	Eigen::Vector3d pi;
	Eigen::Vector3d pf;
	Eigen::Vector3d delta;
	size_t nums;

	std::vector<Eigen::Vector3d> pointsCartesian;
	std::vector<Eigen::Vector4d> pointsMd;

public:
	Hold(const Eigen::Vector3d& pi, const Eigen::Vector3d& pf, size_t nums = 10);

	void fit() = 0;
	void points(CoordinateSystem coordinateSystem) override;
	double length() override;

	Eigen::Vector3d getInitPoint();
	Eigen::Vector3d getTargetPoint();

};
