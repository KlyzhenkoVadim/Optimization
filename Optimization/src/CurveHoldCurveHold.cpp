#include <algorithm>
#include "CurveHoldCurveHold.h"

CurveHoldCurveHold::CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1, double phi1, double R1, double R2, const Eigen::Vector3d& p4,
	double tetta4, double phi4, double betta, double eps, size_t nums) {

	this->p1 = p1;
	this->p4 = p4;
	this->t1 = calcTangentVector(phi1, tetta1);
	this->phi1 = phi1;
	this->phi4 = phi4;
	this->tetta1 = tetta1;
	this->tetta4 = tetta4;
	this->t4 = calcTangentVector(phi4, tetta4);
	this->R1 = R1;
	this->R2 = R2;
	this->eps = eps;
	this->nums = nums;
	this->betta = betta;
	this->pInter = p4 - t4 * betta;
	this->condition = CurveHoldCurveHold::fit();
};

CurveHoldCurveHold::CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1, double phi1, double R1, double R2, const Eigen::Vector3d& pT1,
	const Eigen::Vector3d& pT3, double eps, size_t nums) {

	this->p1 = p1;
	this->p4 = pT3;
	this->pInter = pT1;
	this->phi1 = phi1;
	this->tetta1 = tetta1;
	this->tetta4 = 0.;
	this->phi4 = 0.;
	this->t1 = calcTangentVector(phi1, tetta1);
	Eigen::Vector3d tmpVec = pT3 - pT1;
	this->betta = tmpVec.norm();
	tmpVec.normalize();
	this->t4 = tmpVec;
	this->R1 = R1;
	this->R2 = R2;
	this->nums = nums;
	this->eps = eps;
	this->condition = CurveHoldCurveHold::fit();
};

int CurveHoldCurveHold::fit() 
{
	//ѕроверка на пересечение окружностей
	Eigen::Vector3d b1 = (pInter - p1).cross(t1);
	Eigen::Vector3d n1 = t1.cross(b1);
	n1.normalize();

	Eigen::Vector3d b4 = (p1 - pInter).cross(-t4);
	Eigen::Vector3d n4 = -t4.cross(b4);
	n4.normalize();
	r1 = p1 + R1 * n1;
	r4 = pInter + R2 * n4;

	if (b1.norm() < EPSILON and b4.norm() < EPSILON) {
		t = t1;
		holdLength = (pInter - p1).norm();
		alpha1 = 0;
		alpha2 = 0;
		p1Inter = p1;
		p4Inter = p4;
	}
	
	else
	{

		if ((pInter - r1).norm() < R1 or (p1 - r4).norm() < R2) {
			return -1;
		}

		std::vector<double> alphaFirst(2);
		std::vector<double> alphaSecond(2);

		auto getaAlphaCurrCurvate = [&](bool flag, const auto& p) { // flag = 0 - first curveHold
			if (false == flag) {
				CurveHold firstCurveHold(p1, p, t1, R1);
				if (firstCurveHold.getCondition() == 0)
					return firstCurveHold.getAlpha();
				else
					return -1.;
			}
			else {
				CurveHold secondCurveHold(pInter, p, -t4, R2);
				if (secondCurveHold.getCondition() == 0)
					return secondCurveHold.getAlpha();
				else
					return -1.;
			}
		};

		alphaFirst[0] = getaAlphaCurrCurvate(0, pInter);
		if (alphaFirst[0] < -EPSILON)
			return -1;
		Eigen::Vector3d p1j = p1 + R1 * tan(alphaFirst[0] / 2) * t1;
		alphaSecond[0] = getaAlphaCurrCurvate(1, p1j);
		if (alphaSecond[0] < -EPSILON)
			return -1;
		Eigen::Vector3d pInterj = pInter - R2 * tan(alphaSecond[0] / 2) * t4;

		const size_t maxIter = 1e3;
		size_t iter = 0;

		while (iter < maxIter) {
			alphaFirst[1] = getaAlphaCurrCurvate(0, pInterj);
			if (alphaFirst[1] < -EPSILON)
				return -1;
			p1j = p1 + R1 * tan(alphaFirst[1] / 2) * t1;
			alphaSecond[1] = getaAlphaCurrCurvate(1, p1j);
			if (alphaSecond[1] < -EPSILON)
				return -1;
			pInterj = pInter - R2 * tan(alphaSecond[1] / 2) * t4;

			if (sqrt((alphaFirst[1] - alphaFirst[0]) * (alphaFirst[1] - alphaFirst[0]) + (alphaSecond[1] - alphaSecond[0]) * (alphaSecond[1] - alphaSecond[0])) < eps) {
				break;
			}

			std::reverse(alphaFirst.begin(), alphaFirst.end());
			std::reverse(alphaSecond.begin(), alphaSecond.end());
			iter += 1;  // !!!
		}

		if (iter == 999) {
			return -1;
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
				return -1;
			}
		}
	}
	return 0;
}

int CurveHoldCurveHold::getCondition()
{
	return condition;
}

void CurveHoldCurveHold::getInitPoint(CoordinateSystem coordinateSystem) {
	if (coordinateSystem == CoordinateSystem::CARTESIAN)
		pointInitial = this->p1;
	else
		pointInitialMD = { 0,t1[0],t1[1],t1[2] };
}

void CurveHoldCurveHold::getTarget3Point(CoordinateSystem coordinateSystem) {
	if (coordinateSystem == CoordinateSystem::CARTESIAN)
		pointT3 = this->p4;
	else
		pointMDT3 = { length(),t4[0],t4[1],t4[2] };
}

void CurveHoldCurveHold::getTarget1Point(CoordinateSystem coordinateSystem) {
	if (coordinateSystem == CoordinateSystem::CARTESIAN)
		pointT1 = this->pInter;
	else
		pointMDT1 = { length() - betta,t4[0],t4[1],t4[2] };
}

Eigen::Vector3d CurveHoldCurveHold::FunctionPoint(double md) // md [0,1]
{
	// md*L [0,arc1]
	// md*L (arc1,arc1+HoldLength]
	// md*L (arc1+HoldLength,arc1+HoldLength+arc2]
	// md*L (length-betta,length]
	double L = length();
	double arc1 = alpha1 * R1;
	double arc2 = alpha2 * R2; 
	if (md * L < arc1)
	{
		return p1 + R1 * tan(md * alpha1 / 2) * (t1 + FunctionTangent(md));
	}
	if (0 < md * L - arc1 < holdLength)
	{
		return p1Inter + (md * L - arc1) * t;
	}
	if (0 < md * L - arc1 - holdLength < arc2)
	{
		double s = (md * L - arc1 - holdLength) / arc2;
		return pInter - R2 * tan(s * alpha2 / 2) * (t4 + FunctionTangent(md));
	}
	else
		return p4 - (1 - md) * L * t4;
}

Eigen::Vector3d CurveHoldCurveHold::FunctionTangent(double md) // md [0,1]
{
	double L = length();
	double arc1 = alpha1 * R1, arc2 = alpha2 * R2;
	if (md * L < arc1)
	{
		return t1 * sin((1 - md) * alpha1) / sin(alpha1) + t * sin(md * alpha1) / sin(alpha1);
	}
	if (0 < md * L - arc1 < holdLength)
	{
		return t;
	}
	if (0 < md * L - arc1 - holdLength < arc2)
	{
		return t * sin((1 - md) * alpha2) / sin(alpha2) + t4 * sin(md * alpha2) / sin(alpha2);
	}
	else
		return t4;
}


void CurveHoldCurveHold::points(CoordinateSystem coordinateSystem) {
	double arc1 = R1 * alpha1;
	double arc2 = R2 * alpha2;
	double h = length() / nums;
	int arc1Nums = std::max(10, int(arc1 / h));
	int arc2Nums = std::max(10, int(arc2 / h));
	int nHold1 = std::max(2, int(holdLength / h));
	int nHold2 = std::max(2, int(betta / h));
	if (coordinateSystem == CoordinateSystem::CARTESIAN) {
		std::vector<Eigen::Vector3d> pointsArc1;
		if (abs(alpha1) < EPSILON) {
			pointsArc1.push_back(p1);
		}
		else {
			pointsArc1 = calcInterpolCartesianPoints(p1, t1, t, R1, alpha1, arc1Nums);
		}
		
		std::vector<Eigen::Vector3d> pointsHold1(nHold1 + 1);
		for (size_t idx = 0; idx < nHold1 + 1; ++idx) {
			pointsHold1[idx] = p1Inter + idx * holdLength / nHold1 * t;
		}
		std::vector<Eigen::Vector3d> pointsArc2;
		if (abs(alpha2) < EPSILON) {
			pointsArc2.push_back(p4);
		}
		else {
			pointsArc2 = calcInterpolCartesianPoints(p4Inter, t, t4, R2, alpha2, arc2Nums);
		}
	
		std::vector<Eigen::Vector3d> pointsHold2(nHold2 + 1);
		for (size_t idx = 0; idx < nHold2 + 1; ++idx) {
			pointsHold2[idx] = pInter + idx * betta / nHold2 * t4;
		}
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
		std::vector<Eigen::Vector4d> pointsArc1;
		if (abs(alpha1) < EPSILON) {
			pointsArc1.push_back({ 0,t1[0],t1[1],t1[2] });
		}
		else {
			pointsArc1 = calcInterpolMDPoints(p1, t1, t, R1, alpha1,arc1Nums);
		}
		std::vector<Eigen::Vector4d> pointsHold1(nHold1 + 1);

		std::vector<Eigen::Vector4d> pointsArc2;
		if (abs(alpha2) < EPSILON) {
			pointsArc2.push_back({ 0,t4[0],t4[1],t4[2] });
		}
		else {
			pointsArc2 = calcInterpolMDPoints(p4Inter, t, t4, R2, alpha2, arc2Nums);
		}

		std::vector<Eigen::Vector4d> pointsHold2(nHold2 + 1);

		for (size_t idx = 0; idx < nHold1 + 1; ++idx) {
			pointsHold1[idx] = { arc1 + idx * holdLength / nHold1, t[0], t[1], t[2] };
		}
		for (size_t idx = 0; idx < pointsArc2.size(); ++idx) {
			pointsArc2[idx][0] += arc1 + holdLength;
		}

		for (size_t idx = 0; idx < nHold2 + 1; ++idx) {
			pointsHold2[idx] = { arc1 +arc2 + holdLength + idx * betta/ nHold2, t4[0], t4[1], t4[2] };
		}
		pointsMD = pointsArc1;
		if (holdLength > EPSILON) {
			std::copy(pointsHold1.begin(), pointsHold1.end(), std::back_inserter(pointsMD));
		}
		std::copy(pointsArc2.begin(), pointsArc2.end(), std::back_inserter(pointsMD));
		if (betta > EPSILON) {
			std::copy(pointsHold2.begin(), pointsHold2.end(), std::back_inserter(pointsMD));
		}
	}
}

double CurveHoldCurveHold::length() {
	double arc1 = alpha1 < EPSILON ? 0 : R1 * alpha1;
	double arc2 = alpha2 < EPSILON ? 0 : R2 * alpha2;
	return arc1 + arc2 + holdLength + betta;
}

double CurveHoldCurveHold::getTortuosity()
{
	return alpha1 + alpha2;
}