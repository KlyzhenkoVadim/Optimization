#include "Curve.h"

Curve::Curve(Eigen::Vector3d& pi,double R, double inc1, double azi1, double inc2, double azi2,size_t nums) {
	this->pi = pi;
	this->azi1 = azi1;
	this->inc1 = inc1;
	this->azi2 = azi2;
	this->inc2 = inc2;
	this->t1 = calcTangentVector(azi1, inc1);
	this->t2 = calcTangentVector(azi2, inc2);
	this->nums = nums;
	this->R = R;
}

void Curve::fit() {
	alpha = acos(t1.dot(t2) / t1.norm() / t2.norm());
	pf = pi + R * tan(alpha / 2) * (t1 + t2);
}

void Curve::points(CoordinateSystem coordsys) {
	if (coordsys == CoordinateSystem::CARTESIAN) {
		pointsCartesian = calcInterpolCartesianPoints(pi, t1, t2, R, alpha, nums);
	}
	else {
		pointsMD = calcInterpolMDPoints(pi, t1, t2, R, alpha, nums);
	}
}

double Curve::length() {
	return alpha * R;
}

void Curve::getInitPoint() {
	pointInitial = this->pi;
}

void Curve::getTarget1Point() {
	pointT1 = pf;
}

void Curve::getTarget3Point() {
	pointT3 = pf;
}