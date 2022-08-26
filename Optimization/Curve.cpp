#include "Curve.h"

Curve::Curve(const Eigen::Vector3d& pi, double inc1, double azi1, double inc2, double azi2,double RTVD, TypeCurve type, size_t nums) {
	this->pi = pi;
	this->azi1 = azi1;
	this->inc1 = inc1;
	this->azi2 = azi2;
	this->inc2 = inc2;
	this->t1 = calcTangentVector(azi1, inc1);
	this->t2 = calcTangentVector(azi2, inc2);
	this->nums = nums;
	double dotprod = t1.dot(t2) / t1.norm() / t2.norm();
	dotprod = fabs(dotprod - 1) < EPSILON ? 1 : dotprod;
	alpha = acos(dotprod);
	if (type == TypeCurve::DLS) {
		this->R = RTVD;
		pf = pi + R * tan(alpha / 2) * (t1 + t2);
	}
	else if (type == TypeCurve::TVD) {
		if (abs(t1[2] + t2[2]) < EPSILON) {
			throw(std::runtime_error("points lie on HorizontalPlane"));
		}
		if (abs(alpha) < EPSILON) { //  если alpha == 0 то дуга нулевая, и финальная точка = начальной. вне зависимости от TVD
			this->R = 1 / EPSILON;
			pf = pi;
		}
		else {
			this->R = (RTVD - pi[2]) / (tan(alpha / 2) * (t1[2] + t2[2]));
			pf = pi + R * tan(alpha / 2) * (t1 + t2);
		}
	}
}

void Curve::fit() {}

void Curve::points(CoordinateSystem coordsys) {
	if (coordsys == CoordinateSystem::CARTESIAN) {
		if (abs(alpha) < EPSILON)
			pointsCartesian.push_back(pi);
		else {
			pointsCartesian = calcInterpolCartesianPoints(pi, t1, t2, R, alpha, nums);
		}
	}
	else {
		if (abs(alpha) < EPSILON)
			pointsMD.push_back({ 0,t1[0],t1[1],t1[2] });
		else {
			pointsMD = calcInterpolMDPoints(pi, t1, t2, R, alpha, nums);
		}
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