#include "TrajectoryTemplate.h"

Eigen::Vector3d calcTangentVector(double azimuth, double inclination) {

	double x = sin(inclination * PI / 180.0) * cos(azimuth * PI / 180.0);
	double y = sin(inclination * PI / 180.0) * sin(azimuth * PI / 180.0);
	double z = cos(inclination * PI / 180.0);
	Eigen::Vector3d tangent{ fabs(x) > EPSILON ? x : 0.0,
							 fabs(y) > EPSILON ? y : 0.0,
							 fabs(z) > EPSILON ? z : 0.0
	};
	tangent.normalize();
	return tangent;
};


interpolatedValueType TrajectoryTemplate::calcInterpolPoints(const Eigen::Vector3d& t1, const Eigen::Vector3d& t2, double alpha, size_t nums) {

	std::vector<double> alphaSegment;
	std::vector<Eigen::Vector3d> tInter(nums); 
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

// AllPoint AllLength Solve
double allLength(std::vector<TrajectoryTemplate*>& Well) {
	double trajectoryLength = 0.;
	for (size_t idx = 0; idx < Well.size(); ++idx) {
		trajectoryLength += Well[idx]->length();
	}
	return trajectoryLength;
}

std::vector<Eigen::Vector3d> allPointsCartesian(std::vector<TrajectoryTemplate*>& Well) {
	std::vector<Eigen::Vector3d> pointsCartesian;
	for (size_t idx = 0; idx < Well.size(); ++idx) {
		Well[idx]->points(CoordinateSystem::CARTESIAN);
		pointsCartesian.insert(pointsCartesian.end(), Well[idx]->pointsCartesian.begin(), Well[idx]->pointsCartesian.end());
	}
	return pointsCartesian;
}

std::vector<Eigen::Vector4d> allPointsMD(std::vector<TrajectoryTemplate* >& Well) {
	std::vector<Eigen::Vector4d> pointsMD;
	double lastMD = 0.;
	for (size_t idx = 0; idx < Well.size(); ++idx) {
		Well[idx]->points(CoordinateSystem::MD);
		std::vector<Eigen::Vector4d> tmp = Well[idx]->pointsMD;
		for (size_t i = 0; i < tmp.size(); ++i) {
			pointsMD.push_back({ (lastMD + tmp[i][0]), tmp[i][1], tmp[i][2], tmp[i][3] });
		}
		lastMD = pointsMD.back()[0];
	}
	return pointsMD;
}

int solve(std::vector<TrajectoryTemplate*>& Well) {
	int numErrors = 0;
	for (size_t idx = 0; idx < Well.size(); ++idx) {
		try {
			Well[idx]->fit();
		}
		catch (const std::runtime_error& error) {
			numErrors += 1;
		}
	}
	return numErrors;
}

std::pair<double, double> CartesianToSpherical(Eigen::Vector3d t) {
	t.normalize();
	double theta = acos(t[2]), phi;
	if (theta < EPSILON)
	{
		return { 0,0 };
	}
	theta = theta * 180 / PI;
	if (abs(t[0]) <= EPSILON) 
	{
		phi = t[1] > 0 ? PI / 2 : 3 * PI / 2;
		return { theta,180 / PI * phi };
	}
	
	phi = atan(t[1] / t[0]);
	if (phi > EPSILON) 
	{
		if (t[0] > 0)
			return { theta,phi * 180 / PI };
		else
		{
			phi += PI;
			return { theta,180 / PI * phi };
		}
	}
	if(phi < -EPSILON)
	{
		if (t[0] > 0)
			return { theta,360 + phi * 180 / PI };
		else
		{ 
			phi += PI;
			return{ theta, 180 / PI * phi };
		}
	}
	else
	{
		if (t[0] > 0)
			return { theta,0 };
		else
			return { theta,180 };
	}
}