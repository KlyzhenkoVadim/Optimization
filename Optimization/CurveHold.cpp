#include "CurveHold.h"

CurveHold::CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3, double teta, double phi, double R, size_t nums) {
	this->p1 = p1;
	this->p3 = p3;
	this->teta = teta;
	this->phi = phi;
	this->t1 = calcTangentVector(phi, teta);
	this->R = R;
	this->nums = nums;
	
}

//void CurveHold::fit() {
//	Eigen::Vector3d moduleP { fabs(p1[0] - p3[0]), fabs(p1[1] - p3[1]), fabs(p1[2] - p3[2]) };
//	double psiSqr = moduleP.dot(moduleP);
//	double etta = t1.dot(p3 - p1);
//	double ksi = sqrt(psiSqr - etta * etta);
//
//	if (ksi < EPSILON) { // ksi = 0, значит t1 || p3-p1, следовательно прямая
//		alpha = 0;
//		points.push_back(p1);
//		points.push_back(p3);
//	}
//
//	double betta = psiSqr - 2 * R * ksi > 0.0 ? sqrt(psiSqr - 2 * R * ksi) : 0.0;
//
//	double denominator = 2 * R - ksi;
//
//	if (denominator < EPSILON) {
//		if (etta <= EPSILON) {
//			alpha = PI;
//		}
//		else {
//			alpha = 2 * atan(ksi / (2 * etta));
//		}
//	}
//	else {
//		alpha = 2 * atan((etta - betta) / denominator);
//		if (alpha < EPSILON) {
//			alpha += 2 * PI;
//		}
//	}
//	
//	if (fabs(alpha - PI) < EPSILON) {
//		t2 = -t1;
//	}
//	else {
//		t2 = (p3 - p1 - R * tan(alpha / 2.0) * t1) / (betta + R * tan(alpha / 2.0));
//	}
//
//	for (auto& it : t2) {
//		if (it < EPSILON) {
//			it = 0.0;
//		}
//	}
//}

double CurveHold::length() {
	double arc = alpha * R;
	return arc + betta;
}

double CurveHold::getAlpha() {
	return alpha;
}

void CurveHold::fit() {
	//Eigen::Vector3d moduleP{ fabs(p1[0] - p3[0]), fabs(p1[1] - p3[1]), fabs(p1[2] - p3[2]) };
	Eigen::Vector3d b1 = -t1.cross(p3 - p1);
	Eigen::Vector3d n1 = t1.cross(b1);
	
	if (n1.dot(n1) > EPSILON) {
		n1.normalize();
	}

	Eigen::Vector3d r = p1 + R * n1;
	double norm = (p3 - r).norm();

	if (norm - R < EPSILON) {
		throw std::runtime_error("Target point lies inside  the sphere.");
	}

	double psi = (p3-p1).norm();
	double etta = t1.dot(p3 - p1);
	double ksi = sqrt(psi*psi - etta * etta);

	if (etta <= EPSILON and fabs(ksi - 2 * R) <= EPSILON) {
		throw std::runtime_error("error: alpha >= pi");
	}

	double bettaSqr = psi * psi - 2 * R * ksi;
	betta = 0.0;

	if (bettaSqr > 0.0) {
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

Eigen::Vector3d CurveHold::getInitPoint() {
	return this->p1;
}

Eigen::Vector3d CurveHold::getTargetPoint() {
	return this->p3;
}

void CurveHold::points(CoordinateSystem coordinateSystem) {
	double h = length() / nums;
	double arc = alpha * R;
	Eigen::Vector3d pArcEnd = p3 - betta * t2;

	size_t num = std::max(10, int(arc / h));
	size_t stepHold = std::max(2, int(betta / h));

	if (coordinateSystem == CoordinateSystem::CARTESIAN) {
		std::vector<Eigen::Vector3d> curvePoints = calcInterpolCartesianPoints(p1, t1, t2, R, alpha, num);
		std::vector<Eigen::Vector3d> holdPoints(stepHold + 1);

		for (size_t idx = 0; idx < stepHold + 1; ++idx) {
			holdPoints[idx] = pArcEnd + betta * t2 * idx / stepHold;
		}

		pointsCartesian = curvePoints;

		std::copy(holdPoints.begin(), holdPoints.end(), std::back_inserter(pointsCartesian)); // vstack
	}
	else {
		std::vector<Eigen::Vector4d> curvePoints = calcInterpolMDPoints(p1, t1, t2, R, alpha, num);
		std::vector<Eigen::Vector4d> holdPoints(stepHold + 1);

		for (size_t idx = 0; idx < stepHold + 1; ++idx) {
			holdPoints[idx] = { (arc + betta * idx / stepHold) , t2[0], t2[1], t2[2] };
		}

		pointsMd = curvePoints;

		if (fabs(betta) > EPSILON) {
			std::copy(holdPoints.begin(), holdPoints.end(), std::back_inserter(pointsMd));
		}
	}
}