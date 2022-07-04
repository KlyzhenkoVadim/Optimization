#include <algorithm>
#include "CurveHoldCurveHold.h"

CurveHoldCurveHold::CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1, double phi1, double R1, double R2, const Eigen::Vector3d& p4,
	double tetta4, double phi4, double betta, double eps, size_t nums) {

	this->p1 = p1;
	this->p4 = p4;
	this->t1 = calcTangentVector(tetta1, phi1);
	this->phi1 = phi1;
	this->phi4 = phi4;
	this->tetta1 = tetta1;
	this->tetta4 = tetta4;
	this->t4 = calcTangentVector(tetta4, phi4);
	this->R1 = R1;
	this->R2 = R2;
	this->eps = eps;
	this->nums = nums;
	this->betta = betta;
	this->pInter = p4 - t4 * betta;
}

void CurveHoldCurveHold::fit() {
	//ѕроверка на пересечение окружностей
	Eigen::Vector3d b1 = (pInter - p1).cross(t1);
	Eigen::Vector3d n1 = t1.cross(b1);
	n1.normalize();

	Eigen::Vector3d b4 = (p1 - pInter).cross(-t4);
	Eigen::Vector3d n4 = -t4.cross(b4);
	n4.normalize();

	r1 = p1 + R1 * n1;
	r4 = pInter + R2 * n4;

	auto isSegmentIntersect = [](const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, const Eigen::Vector3d& p4) {

		auto direction = [](const Eigen::Vector3d& pi, const Eigen::Vector3d& pj, const Eigen::Vector3d& pk) {
			Eigen::Vector3d crossProd = (pk - pi).cross(pj - pi);
			return crossProd[2];
		};

		auto onSegment = [](const Eigen::Vector3d& pi, const Eigen::Vector3d& pj, const Eigen::Vector3d& pk) {
			if (((std::min(pi[0], pj[0]) <= pk[0]) and (std::max(pi[0], pj[0]) >= pk[0])) and
				((std::min(pi[1], pj[1]) <= pk[1]) and (std::max(pi[1], pj[1]) >= pk[1])) and
				((std::min(pi[2], pj[2]) <= pk[2]) and (std::max(pi[2], pj[2]) >= pk[2]))) {
				return true;
			}
			return false;
		};

		int d1 = direction(p3, p4, p1);
		int d2 = direction(p3, p4, p2);
		int d3 = direction(p1, p2, p3);
		int d4 = direction(p1, p2, p4);

		if (((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0))) {
			return true;
		}
		else if (d1 == 0 and onSegment(p3, p4, p1) == true) {
			return true;
		}
		else if (d2 == 0 and onSegment(p3, p4, p2) == true) {
			return true;
		}
		else if (d3 == 0 and onSegment(p1, p2, p3) == true) {
			return true;
		}
		else if (d4 == 0 and onSegment(p1, p2, p4) == true) {
			return true;
		}
		else {
			return false;
		}
	};

	if ((pInter - r1).norm() < R1 or (p1 - r4).norm() < R2) {
		throw std::runtime_error("Error: Cannot reach p1/p4 from p4/p1 using Curve-hold!");
	}

	/*
	if (((r1 - r4).norm() < R1 + R2) and true == isSegmentIntersect(r1, r4, p1, pInter)) {
		throw std::runtime_error("Intersection of circles");
	}
	*/

	std::vector<double> alphaFirst(2);
	std::vector<double> alphaSecond(2);

	auto getaAlphaCurrCurvate = [&](bool flag, const auto& p) { // flag = 0 - first curveHold
		if (false == flag) {
			CurveHold firstCurveHold(p1, p, tetta1, phi1, R1);
			firstCurveHold.fit();
			return firstCurveHold.getAlpha();
		}
		else {
			CurveHold secondCurveHold(pInter, p, 180.0 - tetta4, 180.0 + phi1, R2);
			secondCurveHold.fit();
			return secondCurveHold.getAlpha();
		}
	};

	alphaFirst[0] = getaAlphaCurrCurvate(0, pInter);
	Eigen::Vector3d p1j = p1 + R1 * tan(alphaFirst[0] / 2) * t1;
	alphaSecond[0] = getaAlphaCurrCurvate(1, p1j);
	Eigen::Vector3d pInterj = pInter - R2 * tan(alphaSecond[0] / 2) * t4;

	const size_t maxIter = 10e3;
	size_t iter = 0;

	while (iter < maxIter) {
		alphaFirst[1] = getaAlphaCurrCurvate(0, pInterj);
		p1j = p1 + R1 * tan(alphaFirst[1] / 2) * t1;
		alphaSecond[1] = getaAlphaCurrCurvate(1, p1j);
		pInterj = pInter - R2 * tan(alphaSecond[1] / 2) * t4;

		if (sqrt((alphaFirst[1] - alphaFirst[0]) * (alphaFirst[1] - alphaFirst[0]) + (alphaSecond[1] - alphaSecond[0]) * (alphaSecond[1] - alphaSecond[0])) < eps) {
			break;
		}

		std::reverse(alphaFirst.begin(), alphaFirst.end());
		std::reverse(alphaSecond.begin(), alphaSecond.end());
	}

	if (iter == 999){
		throw std::runtime_error("Error: Iterations do not converge!");
	}

	t = (pInterj - p1j);
	t.normalize();
	alpha1 = alphaFirst[1];
	alpha2 = alphaSecond[1];
	p1Inter = p1 + R1 * tan(alpha1 / 2) * (t1 + t);
	p4Inter = pInter - R2 * tan(alpha2 / 2) * (t4 + t);
	holdLength = (p4Inter - p1Inter).norm();

	Eigen::Vector3d tmp = (p4Inter - p1Inter) / holdLength;

	for (size_t idx = 0; idx < t.size(); ++idx) {
		if (fabs(tmp[idx] - t[idx]) > 1e-4) {
			throw std::runtime_error("Error: Cannot reach the Target!");
		}
	}
}

Eigen::Vector3d CurveHoldCurveHold::getInitPoint() {
	return this->p1;
}

Eigen::Vector3d CurveHoldCurveHold::getTargetPoint() {
	return this->p4;
}

void CurveHoldCurveHold::points(CoordinateSystem coordinateSystem) {
	double arc1 = R1 * alpha1;
	double arc2 = R2 * alpha2;
	double h = length() / nums;

	if (coordinateSystem == CoordinateSystem::CARTESIAN) {
		std::vector<Eigen::Vector3d> pointsArc1 = calcInterpolCartesianPoints(p1, t1, t, R1, alpha1, std::max(10, int(arc1 / h)));
		double nHold1 = std::max(2, int(holdLength / h));
		std::vector<Eigen::Vector3d> pointsHold1(nHold1 + 1);

		for (size_t idx = 0; idx < nHold1 + 1; ++idx) {
			pointsHold1[idx] = p1Inter + idx * holdLength / nHold1 * t;
		}

		std::vector<Eigen::Vector3d> pointsArc2 = calcInterpolCartesianPoints(p4Inter, t, t4, R2, alpha2, std::max(10, int(arc2 / h)));


		double nHold2 = std::max(2, int(betta / h));

		std::vector<Eigen::Vector3d> pointsHold2(nHold2 + 1);

		for (size_t idx = 0; idx < nHold2 + 1; ++idx) {
			pointsHold2[idx] = pInter + idx * betta / nHold2 * t4;
		}
		/*
		pointsCartesian.pointsArc1 = pointsArc1;
		pointsCartesian.pointsHold1 = pointsHold1;
		pointsCartesian.pointsArc2 = pointsArc2;
		pointsCartesian.pointsHold2 = pointsHold2;*/
		pointsCartesian = pointsArc1;
			if (holdLength > EPSILON) {
				std::copy(pointsHold1.begin(), pointsHold1.end(), std::back_inserter(pointsCartesian));
			}
			std::copy(pointsArc2.begin(), pointsArc2.end(), std::back_inserter(pointsCartesian));
			if (betta > EPSILON) {
				std::copy(pointsHold2.begin(), pointsHold2.end(), std::back_inserter(pointsCartesian));
			}
	}
	else {
		std::vector<Eigen::Vector4d> pointsArc1 = calcInterpolMDPoints(p1, t1, t, R1, alpha1, std::max(10, int(arc1 / h)));
		double nHold1 = std::max(2, int(holdLength / h));
		std::vector<Eigen::Vector4d> pointsHold1(nHold1 + 1);

		std::vector<Eigen::Vector4d> pointsArc2 = calcInterpolMDPoints(p4Inter, t, t4, R2, alpha2, std::max(10, int(arc2 / h)));
		
		double nHold2 = std::max(2, int(betta / h));

		std::vector<Eigen::Vector4d> pointsHold2(nHold2 + 1);

		for (size_t idx = 0; idx < nHold1 + 1; ++idx) {
			pointsHold1[idx] = { arc1 + idx * holdLength / nHold1, t[0], t[1], t[2] };
		}
		for (size_t idx = 0; idx < pointsArc2.size(); ++idx) {
			pointsArc2[idx][0] += arc1 + holdLength;
		}

		for (size_t idx = 0; idx < nHold2 + 1; ++idx) {
			pointsHold2[idx] = { arc1 + idx * holdLength / nHold1, t4[0], t4[1], t4[2] };
		}
		/*
		pointsMD.pointsArc1 = pointsArc1;
		pointsMD.pointsHold1 = pointsHold1;
		pointsMD.pointsArc2 = pointsArc2;
		pointsMD.pointsHold2 = pointsHold2;*/
		pointsMd = pointsArc1;
		if (holdLength > EPSILON) {
			std::copy(pointsHold1.begin(), pointsHold1.end(), std::back_inserter(pointsMd));
		}
		std::copy(pointsArc2.begin(), pointsArc2.end(), std::back_inserter(pointsMd));
		if (betta > EPSILON) {
			std::copy(pointsHold2.begin(), pointsHold2.end(), std::back_inserter(pointsMd));
		}
	}
}


double CurveHoldCurveHold::length() {
	double arc1 = R1 * alpha1;
	double arc2 = R2 * alpha2;
	return arc1 + arc2 + holdLength + betta;
}