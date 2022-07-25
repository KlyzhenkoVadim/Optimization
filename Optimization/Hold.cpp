#include "Hold.h"

Hold::Hold(const Eigen::Vector3d& pi,const Eigen::Vector3d& pf, size_t nums) {
	this->pi = pi;
	this->pf = pf;
	this->nums = nums;
	this->delta = pf - pi;
}

double Hold::length() {
	return delta.norm();
}

void Hold::getInitPoint() {
	pointInitial = this->pi;
}

void Hold::getTarget3Point() {
	pointT3 = this->pf;
}
void Hold::getTarget1Point() {
	pointT1 = this->pf;
}

void Hold::fit() {}

void Hold::points(CoordinateSystem coordinateSystem) {
	if (coordinateSystem == CoordinateSystem::CARTESIAN) {
		for (size_t i = 0; i < nums; ++i) {
			pointsCartesian.push_back(pi + delta * i / (nums - 1));
		}
	}
	else {
		Eigen::Vector3d t = delta;
		t.normalize();
		for (size_t id = 0; id < nums; ++id) {
			pointsMD.push_back({ (delta.norm() * id / (nums - 1)), t[0], t[1], t[2] });
		}
	}
}