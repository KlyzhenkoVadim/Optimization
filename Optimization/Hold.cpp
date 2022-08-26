#include "Hold.h"

Hold::Hold(const Eigen::Vector3d& pi,const Eigen::Vector3d& pf, size_t nums) {
	this->pi = pi;
	this->pf = pf;
	this->nums = nums;
	this->direction = pf - pi;
	this->Length = direction.norm();
	direction.normalize();


}

Hold::Hold(const Eigen::Vector3d& pi, double inc, double azi, double L, typeHold type, size_t nums) {
	this->pi = pi;
	this->nums = nums;
	this->direction = calcTangentVector(azi, inc);
	if (type == typeHold::md) {
		this->Length = L;
	}
	else if (type == typeHold::TVD) {
		if (abs(direction[2]) < EPSILON) {
			throw(std::runtime_error("HoldTVD is on HorizontalPlane"));
		}
		this->Length = (L - pi[2]) / direction[2];	
	}
	pf = pi + Length * direction;
}

double Hold::length() {
	return Length;
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
		if (Length < EPSILON)
			pointsCartesian.push_back(pi);
		else {
			for (size_t i = 0; i < nums; ++i) {
				pointsCartesian.push_back(pi + Length * direction * i / (nums - 1));
			}
		}
	}

	else {
		if (direction.norm() > EPSILON) {
			for (size_t id = 0; id < nums; ++id) {
				pointsMD.push_back({ (Length * id / (nums - 1)), direction[0], direction[1], direction[2] });
			}
		}
	}
}