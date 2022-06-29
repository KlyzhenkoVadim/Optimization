#pragma once
#include <vector>
#include "Eigen/Dense"

// ��� ����������� ��������� �� ������: A compendium of directional calculations based on minimum curvature method. Sawaryn, Thorogood.

enum class CoordinateSystem { CARTESIAN, MD };

constexpr double PI = 3.14159265358979323846;
constexpr double EPSILON = 10e-10;

using interpolatedValueType = std::pair<std::vector<double>, std::vector<Eigen::Vector3d>>;

class TrajectoryTemplate
{
private:
	interpolatedValueType calcInterpolPoints(const Eigen::Vector3d& t1, const Eigen::Vector3d& t2, double alpha, size_t nums = 100);
public:
	virtual void fit() = 0;
	virtual void points(CoordinateSystem coordinateSystem) = 0;
	virtual double length() = 0;

	Eigen::Vector3d calcTangentVector(double azimuth, double inclination);

	std::vector<Eigen::Vector3d> calcInterpolCartesianPoints(const Eigen::Vector3d& p1,
		const Eigen::Vector3d& t1,
		const Eigen::Vector3d& t2,
		double R, double alpha, size_t nums = 100);

	std::vector<Eigen::Vector4d> calcInterpolMDPoints(const Eigen::Vector3d& p1,
		const Eigen::Vector3d& t1,
		const Eigen::Vector3d& t2,
		double R, double alpha, size_t nums = 100);
};	

