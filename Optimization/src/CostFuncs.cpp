#include "CostFuncs.h"

double sigmoid(double x, double penalty, double alpha, double x0) {
	return penalty / (1 + exp(alpha * (x - x0)));
}

int signum(double x) {
	if (x > 0)
		return 1;
	if (x < 0)
		return -1;
	return 0;
}

double sepFactor(const std::vector<Eigen::Vector3d>& pCartesianW1,
	const std::vector<Eigen::Vector4d>& pMDW1,
	const std::vector<Eigen::Vector3d>& pCartesianW2,
	const std::vector<Eigen::Vector4d>& pMDW2,
	double TVDstart ,bool actFunc, double penalty) 
{
	size_t n = pCartesianW1.size(), m = pCartesianW2.size();
	std::vector<std::vector<double>> distanceMatrix(n,std::vector<double>(m,1e3));
	std::vector<std::vector<double>> scalarProdMatrix(n, std::vector<double>(m,0));
	double sgm = 9e-3;
	sgm = sgm * 2.; // TODO: DELETE
	double eps = 1e-6;
	bool flg = true;
	double sepFactor = 1e3;
	size_t StartW1 = n - 1;
	size_t StartW2 = m - 1;
	size_t k;
	Eigen::Vector3d diffVector, tangentW1;

	for (size_t i = 0; i < n; ++i) {
		scalarProdMatrix[i][0] = -1.;
	}	
	
	for (size_t i = 0; i < n; ++i) {
		if (pCartesianW1[i][2] > TVDstart) {
			StartW1 = i;
			break;
		}
	}
	
	for (size_t j = 0; j < m; ++j) {
		if (pCartesianW2[j][2] > TVDstart){
			StartW2 = j;
			break;
		}
	}

	k = StartW2 + 1; // k = 1;
	for (size_t idxN = StartW1; idxN < n; ++idxN) {
		for (size_t idxM = k - 1; idxM < m; ++idxM) {
			diffVector = pCartesianW2[idxM] - pCartesianW1[idxN];
			tangentW1 = { pMDW1[idxN][1], pMDW1[idxN][2], pMDW1[idxN][3] };
			double mdsum = pMDW1[idxN][0] + pMDW2[idxM][0];
			if (diffVector.norm() < eps) {
				scalarProdMatrix[idxN][idxM] = 0;
				distanceMatrix[idxN][idxM] = 0.;
			}
			else {
				scalarProdMatrix[idxN][idxM] = round(tangentW1.dot(diffVector) / tangentW1.norm() / diffVector.norm() * 100000) / 100000;
				if (idxM > 1 && (signum(scalarProdMatrix[idxN][idxM]) != signum(scalarProdMatrix[idxN][idxM - 1]))) {
					if (mdsum > eps)
						distanceMatrix[idxN][idxM] = diffVector.norm() / (sgm * mdsum);
					if (flg) {
						k = idxM;
						flg = false;
					}
				}
			}
		}
		flg = true;
	}
		for (size_t i = 0; i < n; ++i){
			for (size_t j = 0; j < m; ++j) {
				sepFactor = std::min(distanceMatrix[i][j], sepFactor);
			}
		}

		if (actFunc) {
			return sigmoid(sepFactor, penalty, 10., 1.5);//return sigmoid(sepFactor, penalty, 100, 1.5);
		}
		return sepFactor;
	}

double AHD(std::vector<double>& pX, std::vector<double>& pY) {
	double AHD = 0;
	for (size_t i = 1; i < pX.size(); ++i) {
		double dx = pX[i] - pX[i - 1];
		double dy = pY[i] - pY[i - 1];
		AHD += sqrt(dx*dx + dy*dy);
	}
	return AHD;
}

double dls(Eigen::Vector3d& tangent1, Eigen::Vector3d& tangent2) {
	if (tangent1.norm() < EPSILON || tangent2.norm() < EPSILON) {
		return 0.;
	}
	double dotProd= tangent1.dot(tangent2) / tangent1.norm() / tangent2.norm();
	if (abs(dotProd) > 1)
		return 0;
	return (180. / PI) * acos(dotProd);
}

double Tortuosity(std::vector<TrajectoryTemplate*>& well) {
	double tortuos = 0.;
	for (size_t i = 0; i < well.size(); ++i) 
	{
		tortuos += well[i]->getTortuosity();
	}
	return 180/PI*tortuos;
}

double DDI(std::vector<TrajectoryTemplate*>& well,const std::vector<Eigen::Vector3d>& pCartesian, bool actFunc, double penalty) 
{
	double tortuos = Tortuosity(well); 
	double DDI;
	std::vector<double> pX, pY;
	size_t size = pCartesian.size();
	pX.reserve(size);
	pY.reserve(size);
	for (size_t idx = 0; idx < size; ++idx) {
		auto p = pCartesian[idx];
		pX.push_back(p.x());
		pY.push_back(p.y());
	}
	double ahd = AHD(pX,pY);
	double MD = allLength(well);
	double TVD = pCartesian.back()[2] - pCartesian[0][2] ;
	DDI = log10((1. / 0.305) * ahd * MD * tortuos / TVD);
	if (actFunc) {
		return sigmoid(DDI, penalty, -2.5, 6.25); // ???
	}
	return DDI;
}	
/*
std::vector<Eigen::Vector3d> p1{
	{0,0,1},{0,0,2}
};
std::vector<Eigen::Vector3d> p2{
	{0,0,0,},{1,0,0}
};
assert(ERD(p1) < EPSILON);
assert(std::abs(ERD(p2)-1.0) < EPSILON);
*/
double ERD(const std::vector<Eigen::Vector3d>& points)
{
	size_t size = points.size();
	std::vector<double> x, y;
	x.reserve(size);
	y.reserve(size);
	double max_z = 0.;
	for (const auto& p : points)
	{
		x.emplace_back(p.x());
		y.emplace_back(p.y());
		if (p.z() > max_z)
		{
			max_z = p.z();
		}
	}
	double ahd = AHD(x, y);
	double tvd = max_z - points.front().z();
	double res = tvd < std::numeric_limits<double>::epsilon() ? EPSILON : ahd / tvd;
	return res;
}

double orderScore1(std::vector<TrajectoryTemplate*>& mainWell, std::vector<std::vector<Eigen::Vector3d>>& pCTrajectories, 
	std::vector<std::vector<Eigen::Vector4d>>& pMDTrajectories,double SepFactorShift, double penalty) {
	double mainLength = 0., mainDDI = 0., rSepFactor = 0.;
	int condition = solve(mainWell);
	if (condition != 0) {
		return penalty *(-condition) / mainWell.size(); // penalty * percent of incorrect templates.
	}
	mainLength = allLength(mainWell); 
	std::vector<Eigen::Vector3d> mainPCartesian = allPointsCartesian(mainWell);
	std::vector<Eigen::Vector4d> mainPMD = allPointsMD(mainWell);
	double IdealLength;
	mainWell[0]->getInitPoint();
	mainWell.back()->getTarget1Point();
	mainWell.back()->getTarget3Point();
	IdealLength = (mainWell.back()->pointT1 - mainWell[0]->pointInitial).norm() + (mainWell.back()->pointT3 - mainWell.back()->pointT1).norm();
	if (std::abs(IdealLength) < std::numeric_limits<double>::epsilon())
		return 0.;
	mainDDI = DDI(mainWell,mainPCartesian);
	size_t size = pCTrajectories.size();

	for (size_t idx = 0; idx < size; ++idx) {
		auto p = pCTrajectories[idx];
		auto incl = pMDTrajectories[idx];
		rSepFactor += std::max(
			sepFactor(mainPCartesian,mainPMD,p,
			incl, SepFactorShift),
			sepFactor(p,incl,mainPCartesian,
				mainPMD,SepFactorShift));
	}

	double mainERD = ERD(mainPCartesian);

	return mainLength / IdealLength + rSepFactor + /*mainERD*/mainDDI;
}

double scoreSolver(std::vector<TrajectoryTemplate*>& tmp, const WellTrajectoryConstraints& cs, double penalty)
{
		double mainLength, mainDDI, IdealLength;;
		int condition = solve(tmp);
		if (condition != 0)
		{
			return penalty *(-condition)/tmp.size(); // penalty * percent of incorrect templates.
		}
		mainLength = allLength(tmp);
		std::vector<Eigen::Vector3d> mainPCartesian = allPointsCartesian(tmp);
		std::vector<Eigen::Vector4d> mainPMD = allPointsMD(tmp);
		tmp[0]->getInitPoint();
		tmp.back()->getTarget1Point();
		tmp.back()->getTarget3Point();
		IdealLength = (tmp.back()->pointT1 - tmp[0]->pointInitial).norm() + (tmp.back()->pointT3 - tmp.back()->pointT1).norm();
		if (IdealLength == 0)
			return 0.;
		mainDDI = DDI(tmp,mainPCartesian);
		double lenPen = PenaltyLength(mainLength,cs.maxMD), ahdPen = PenaltyAHDNSEW(mainPCartesian, cs.maxDistEastWest, cs.maxDistNorthSouth, 100);

		return mainLength / IdealLength + mainDDI + ahdPen + lenPen;
}