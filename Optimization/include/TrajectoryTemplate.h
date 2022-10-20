#pragma once
#include <vector>
#include "Eigen/Dense"
#include <iostream>

// все обозначения сохранены из статьи: A compendium of directional calculations based on minimum curvature method. Sawaryn, Thorogood.

Eigen::Vector3d calcTangentVector(double azimuth, double inclination);

std::pair<double, double> CartesianToSpherical(Eigen::Vector3d t);

enum class CoordinateSystem { CARTESIAN, MD };

constexpr double PI = 3.14159265358979323846;
constexpr double EPSILON = 10e-10;

using interpolatedValueType = std::pair<std::vector<double>, std::vector<Eigen::Vector3d>>;

class TrajectoryTemplate
{
private:
	interpolatedValueType calcInterpolPoints(const Eigen::Vector3d& t1, const Eigen::Vector3d& t2, double alpha, size_t nums = 100);
public:
	virtual int getCondition() = 0;
	virtual void points(CoordinateSystem coordinateSystem) = 0;
	virtual double length() = 0;
	virtual double getTortuosity() = 0;
	virtual double getAHD() = 0;
	virtual void getInitPoint(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) = 0;
	virtual void getTarget1Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) = 0;
	virtual void getTarget3Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) = 0;
	std::vector<Eigen::Vector3d> pointsCartesian;
	std::vector<Eigen::Vector4d> pointsMD;

	Eigen::Vector3d pointInitial,pointT1, pointT3;
	Eigen::Vector4d pointInitialMD, pointMDT1, pointMDT3;
	std::vector<Eigen::Vector3d> calcInterpolCartesianPoints(const Eigen::Vector3d& p1,
		const Eigen::Vector3d& t1,
		const Eigen::Vector3d& t2,
		double R, double alpha, size_t nums = 100);

	std::vector<Eigen::Vector4d> calcInterpolMDPoints(const Eigen::Vector3d& p1,
		const Eigen::Vector3d& t1,
		const Eigen::Vector3d& t2,
		double R, double alpha, size_t nums = 100);
};	


double allLength(std::vector<TrajectoryTemplate*>& Well);

std::vector<Eigen::Vector3d> allPointsCartesian(std::vector<TrajectoryTemplate*>& Well);

std::vector<Eigen::Vector4d> allPointsMD(std::vector<TrajectoryTemplate* >& Well);

int solve(std::vector<TrajectoryTemplate*>& Well);
