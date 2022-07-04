#include "TrajectoryTemplate.h"

Eigen::Vector3d TrajectoryTemplate::calcTangentVector(double azimuth, double inclination) {

	double x = sin(inclination * PI / 180.0) * cos(azimuth * PI / 180.0);
	double y = sin(inclination * PI / 180.0) * sin(azimuth * PI / 180.0);
	double z = cos(inclination * PI / 180.0);

	return Eigen::Vector3d{ fabs(x) > EPSILON ? x : 0.0,
							 fabs(y) > EPSILON ? y : 0.0,
							 fabs(z) > EPSILON ? z : 0.0
	};
}
// 
interpolatedValueType TrajectoryTemplate::calcInterpolPoints(const Eigen::Vector3d& t1, const Eigen::Vector3d& t2, double alpha, size_t nums) {

	std::vector<double> alphaSegment;
	std::vector<Eigen::Vector3d> tInter(nums); // здесь можно просто Eigen::Vector3d tInter;
	double eps = 1e-3;

	tInter[0] = t1; 

	for (size_t i = 0; i < nums; ++i) {
		double step = double(i) / (nums - 1);
		alphaSegment.push_back(alpha * step);
	}

	for (size_t idx = 1; idx < nums; ++idx) {
		tInter[idx] = t1 * sin(alpha - alphaSegment[idx]) / sin(alpha) + t2 * sin(alphaSegment[idx]) / sin(alpha);
		tInter[idx].normalize();
	}

	for (auto& vec3d : tInter) {
		for (auto& value : vec3d) {
			if (fabs(value) < eps) {
				value = 0.0;
			}
		}
	}

	return { alphaSegment, tInter };
}

std::vector<Eigen::Vector3d> TrajectoryTemplate::calcInterpolCartesianPoints(const Eigen::Vector3d& p1,
	const Eigen::Vector3d& t1,
	const Eigen::Vector3d& t2,
	double R, double alpha, size_t nums) {

	std::vector<Eigen::Vector3d> pointsCartesian(nums);

	pointsCartesian[0] = p1;

	interpolatedValueType interValues = calcInterpolPoints(t1, t2, alpha, nums);

	const auto& [alphaInter, tInter] = interValues;

	for (size_t idx = 1; idx < nums; ++idx) {
		pointsCartesian[idx] = pointsCartesian[0] + R * tan(alphaInter[idx] / 2.0) * (tInter[idx] + tInter[0]);
	}

	return pointsCartesian;
}

std::vector<Eigen::Vector4d> TrajectoryTemplate::calcInterpolMDPoints(const Eigen::Vector3d& p1,
	const Eigen::Vector3d& t1,
	const Eigen::Vector3d& t2,
	double R, double alpha, size_t nums) {

	std::vector<Eigen::Vector4d> pointsMD(nums);

	interpolatedValueType interValues = calcInterpolPoints(t1, t2, alpha, nums);

	const auto& [alphaInter, tInter] = interValues;

	for (size_t idx = 0; idx < nums; ++idx) {
		pointsMD[idx][0] = R * alphaInter[idx];
		pointsMD[idx][1] = tInter[idx][0];
		pointsMD[idx][2] = tInter[idx][1];
		pointsMD[idx][3] = tInter[idx][2];
	}

	return pointsMD;
}