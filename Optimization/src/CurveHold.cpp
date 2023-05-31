#include "CurveHold.h"

CurveHold::CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3, double teta, double phi, double R, size_t nums) {
	this->p1 = p1;
	this->p3 = p3;
	this->teta = teta;
	this->phi = phi;
	this->t1 = calcTangentVector(phi, teta);
	this->R = R;
	this->nums = nums;
	this->condition = CurveHold::fit();
};

CurveHold::CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3, const Eigen::Vector3d& t1, double R, size_t nums) {
	this->p1 = p1;
	this->p3 = p3;
	this->t1 = t1;
	this->R = R;
	this->nums = nums;
	this->condition = CurveHold::fit();
}

double CurveHold::length() {
	double arc = alpha * R;
	return arc + betta;
}

double CurveHold::getTortuosity()
{
	return alpha;
}

double CurveHold::getAlpha() {
	return alpha;
}

Eigen::Vector3d CurveHold::getPointInterpol() {
	return p1+R*tan(alpha/2)*(t1+t2);
}

Eigen::Vector3d CurveHold::getTangent2() {
	return t2;
}

int CurveHold::fit() {
	
	Eigen::Vector3d b1 = (p3 - p1).cross(t1);
	Eigen::Vector3d n1 = t1.cross(b1);
	n1.normalize();
	Eigen::Vector3d r = p1 + R * n1;
	double norm = (p3 - r).norm();
	if (b1.norm() < EPSILON) {
		t2 = t1;
		alpha = 0;
		betta = (p3 - p1).norm();
	}
	else
	{
		if ((norm - R) < -EPSILON)
		{
			return -1;
		}

		double psi = (p3 - p1).norm();
		double etta = t1.dot(p3 - p1);
		double ksi = sqrt(psi * psi - etta * etta);

		if (etta < EPSILON && ksi - 2 * R < EPSILON) {
			return -1;
		}

		double bettaSqr = psi * psi - 2 * R * ksi;
		betta = 0.0;

		if (bettaSqr > 0.) {
			betta = sqrt(bettaSqr);
		}

		if (fabs(2 * R - ksi) < EPSILON) {
			alpha = 2 * atan(0.5 * ksi / etta);
		}
		else {
			alpha = 2 * atan((etta - betta) / (2 * R - ksi));
		}

		t2 = (p3 - p1 - R * tan(alpha / 2.0) * t1) / (betta + R * tan(alpha / 2.0));
		t2.normalize();

		for (auto& it : t2) {
			if (fabs(it) < EPSILON) {
				it = 0.0;
			}
		}
	}
	return 0;
}

Eigen::Vector3d CurveHold::FunctionPoint(double md) // md [0,1]
{
	double L = length();
	double s = md * L / alpha / R;
	if (s - 1 > EPSILON)
	{
		if (betta < EPSILON)
			return p3;
		return p3 - (1 - md) * L  * t2;
	}
	Eigen::Vector3d tstar = FunctionTangent(md);
	return p1 + R * tan(s*alpha/2) * (t1 + tstar);
}

Eigen::Vector3d CurveHold::FunctionTangent(double md) // md [0,1]
{
	double L = length(), s = md * L / alpha / R;
	if (s - 1 > EPSILON)
	{
		return t2;
	}
	return t1* sin((1 - s) * alpha) / sin(alpha) + t2 * sin(s * alpha) / sin(alpha);

}

int CurveHold::getCondition()
{
	return condition;
}

void CurveHold::getInitPoint(CoordinateSystem coordinateSystem) {
	if (coordinateSystem == CoordinateSystem::CARTESIAN)
		pointInitial = this->p1;
	else
		pointInitialMD = { 0,t1[0],t1[1],t1[2] };

}
void CurveHold::getTarget3Point(CoordinateSystem coordinateSystem) {
	if (coordinateSystem == CoordinateSystem::CARTESIAN) {
		pointT3 = this->p3;
	}
	else {
		pointMDT3 = { length(),t2[0],t2[1],t2[2] };
	}
}
void CurveHold::getTarget1Point(CoordinateSystem coordinateSystem) {
	if (coordinateSystem == CoordinateSystem::CARTESIAN) {
		pointT1 = this->p3;
	}
	else {
		pointMDT1 = { length(),t2[0],t2[1],t2[2] };
	}
}
void CurveHold::points(CoordinateSystem coordinateSystem) {
	double h = length() / nums;
	double arc = alpha * R;
	Eigen::Vector3d pArcEnd = p3 - betta * t2;

	size_t num = std::max(10, int(arc / h));
	size_t stepHold = std::max(2, int(betta / h));
	if (coordinateSystem == CoordinateSystem::CARTESIAN) {
		std::vector<Eigen::Vector3d> curvePoints;
		if (abs(alpha) < EPSILON) {
			curvePoints.push_back(p1);
		}
		else {
			curvePoints = calcInterpolCartesianPoints(p1, t1, t2, R, alpha, num);
		}
		std::vector<Eigen::Vector3d> holdPoints(stepHold + 1);
		for (size_t idx = 0; idx < stepHold + 1; ++idx) {
			holdPoints[idx] = pArcEnd + betta * t2 * idx / stepHold;
		}
			pointsCartesian = curvePoints;
		std::copy(holdPoints.begin(), holdPoints.end(), std::back_inserter(pointsCartesian)); // vstack
	}
	else {
		std::vector<Eigen::Vector4d> curvePoints;
		if (abs(alpha) < EPSILON) {
			curvePoints.push_back({ 0,t1[0],t1[1],t1[2] });
		}
		else {
			curvePoints = calcInterpolMDPoints(p1, t1, t2, R, alpha, num);
		}
		std::vector<Eigen::Vector4d> holdPoints(stepHold + 1);
		for (size_t idx = 0; idx < stepHold + 1; ++idx) {
			holdPoints[idx] = { (arc + betta * idx / stepHold) , t2[0], t2[1], t2[2] };
		}
		pointsMD = curvePoints;
		if (fabs(betta) > EPSILON) {
			std::copy(holdPoints.begin(), holdPoints.end(), std::back_inserter(pointsMD));
		}
	}
}