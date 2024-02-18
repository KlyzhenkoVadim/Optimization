#include "Penalties.h"

double PenaltyDLS(std::vector<TrajectoryTemplate*>& Well, double penalty) {
	double alpha, dotprod, dls;
	int s = solve(Well);
	for (size_t i = 0; i < Well.size() - 1; ++i) {
		Well[i]->getInitPoint(CoordinateSystem::MD);
		Well[i]->getTarget1Point(CoordinateSystem::MD);
        Eigen::Vector3d tprev{Well[i]->pointInitialMD[1],
                                Well[i]->pointInitialMD[2],
                                Well[i]->pointInitialMD[3]};
        Eigen::Vector3d tnext{Well[i]->pointMDT1[1],
                                Well[i]->pointMDT1[2],
                                Well[i]->pointMDT1[3]};
		double length = Well[i]->pointMDT1[0] - Well[i]->pointInitialMD[0];
		dotprod = tprev.dot(tnext) / tprev.norm() / tnext.norm();
		if (length < EPSILON)
			continue;
		dotprod = abs(dotprod - 1) < EPSILON ? 1 : dotprod;
		alpha = 180 / PI * acos(dotprod);
		dls = 10 * alpha / length;
		if (dls > 1.5) {
			return penalty;
		}
	}
	return 0;
};

double PenaltyIncWell(const std::vector<Eigen::Vector4d>& pMD, double incMin,
                      double incMax, double penalty) {
	double value = 0, currInc = -1;
	size_t size = pMD.size();
	for (size_t i = 0; i < size; ++i)
	{
		auto p = pMD[i];
		currInc = cartesianToSpherical(Eigen::Vector3d{ p[1],p[2],p[3] }).first;
		if (currInc - incMin > -EPSILON && currInc - incMax < EPSILON)
			continue;
		value += penalty / size;
	}
	return value;
}

double PenaltyIncTarget(std::vector<TrajectoryTemplate*>& Well, double penalty) {
	int s = solve(Well);
	Well.back()->getTarget1Point(CoordinateSystem::MD);
    Eigen::Vector3d tangent = {Well.back()->pointMDT1[1],
                                Well.back()->pointMDT1[2],
                                Well.back()->pointMDT1[3]};
	std::pair<double, double> incAzi = cartesianToSpherical(tangent);
	if (incAzi.first - 40 > EPSILON) {
		return penalty;
	}
	return 0;
};

double PenaltyAlphaTarget(std::vector<TrajectoryTemplate*>& Well, double penalty) {
	int s = solve(Well);
	Well.back()->getInitPoint(CoordinateSystem::MD);
	Well.back()->getTarget1Point(CoordinateSystem::MD);
    Eigen::Vector3d t1 = {Well.back()->pointInitialMD[1],
                            Well.back()->pointInitialMD[2],
                            Well.back()->pointInitialMD[3]};
    Eigen::Vector3d t2 = {Well.back()->pointMDT1[1],
                            Well.back()->pointMDT1[2],
                            Well.back()->pointMDT1[3]};
	double dotprod = t1.dot(t2) / t1.norm() / t2.norm();
	dotprod = fabs(dotprod - 1) < EPSILON ? 1 : dotprod;
	double alpha = acos(dotprod);
	if (alpha - PI / 2 > EPSILON) {
		return penalty;
	}
	return 0;
}

int findTVD(const std::vector<Eigen::Vector3d>& pC, double TVD) {
	int n = pC.size();
	int l = 0, r = n, mid = (l + r) / 2;
	while (!(l >= r)) {
		mid = (l + r) / 2;
		if (abs(pC[mid][2] - TVD) < EPSILON)
			return mid;
		if (pC[mid][2] - TVD > EPSILON) {
			r = mid;
		}
		else {
			l = mid + 1;
		}
	}
	return -(l + 1);
}

double PenaltyConstraint(const std::vector<Constraint>& cs,
                         const std::vector<Eigen::Vector3d>& pC,
                         const std::vector<Eigen::Vector4d>& pMD,
                         double penalty) {
	int id = 0;
	double pen = 0;
	std::pair<double, double>  incAzi = { 0,0 };
	for (size_t i = 0; i < cs.size(); ++i) {
		id = findTVD(pC, cs[i].lMin.TVD);
		if (id < 0) {
			id = -id - 2;
		}
		Eigen::Vector3d tmpVec = { pMD[id][1],pMD[id][2],pMD[id][3] };
		incAzi = cartesianToSpherical(tmpVec);
        if (!(incAzi.first - cs[i].lMin.theta > -EPSILON &&
                incAzi.first - cs[i].lMax.theta < EPSILON)) {
                pen += penalty / cs.size() / 2;
        }
        if ((incAzi.first > EPSILON) &&
            !((incAzi.second - cs[i].lMin.phi > -EPSILON &&
                incAzi.second - cs[i].lMax.phi < EPSILON) ||
                (incAzi.second - (cs[i].lMin.phi + 180) > -EPSILON &&
                incAzi.second - (cs[i].lMax.phi + 180) < EPSILON))) {
                pen += penalty / cs.size() / 2;
        }
	}
	return pen;
}

double PenaltyLength(double length, double MaxLength, double penalty) {
	double answer = length > MaxLength ? penalty : 0;
	return answer;
};

double PenaltyAHDNSEW(const std::vector<Eigen::Vector3d>& pC, double EWMAX,
                      double NSMAX, double penalty) {
	double value = 0, EW = 0, NS = 0;
	for (size_t i = 1; i < pC.size(); ++i)
	{
		NS += fabs(pC[i][0] - pC[i - 1][0]);
		EW += fabs(pC[i][1] - pC[i - 1][1]);
	}
	if (EW - EWMAX > EPSILON || NS - NSMAX > EPSILON)
	{
		value = penalty;
	}
	return value;
}