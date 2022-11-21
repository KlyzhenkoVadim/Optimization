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

double sepFactor(std::vector<Eigen::Vector3d>& pCartesianW1,
	std::vector<Eigen::Vector4d>& pMDW1,
	std::vector<Eigen::Vector3d>& pCartesianW2,std::vector<Eigen::Vector4d>& pMDW2,
	double TVDstart ,bool actFunc, double penalty) {
	size_t n = pCartesianW1.size(), m = pCartesianW2.size();
	std::vector<std::vector<double>> distanceMatrix(n,std::vector<double>(m,1e3));
	std::vector<std::vector<double>> scalarProdMatrix(n, std::vector<double>(m,0));
	double sgm = 9e-3;
	double eps = 1e-6;
	bool flg = true;
	double sepFactor = 1e3;
	size_t StartW1 = n - 1;
	size_t StartW2 = m - 1;
	size_t k;
	Eigen::Vector3d diffVector, tangentW1;// tangentW2;

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
			//tangentW2 = { pMDW2[idxM][1],pMDW2[idxM][2],pMDW2[idxM][3] };
			if (diffVector.norm() < eps) {
				scalarProdMatrix[idxN][idxM] = 0;
				distanceMatrix[idxN][idxM] = 0.;
			}
			else {
				scalarProdMatrix[idxN][idxM] = round(tangentW1.dot(diffVector) / tangentW1.norm() / diffVector.norm() * 100000) / 100000;
				if (idxM > 1 and (signum(scalarProdMatrix[idxN][idxM]) != signum(scalarProdMatrix[idxN][idxM - 1]))) {
					if (pMDW1[idxN][0] + pMDW2[idxM][0] > eps)
						distanceMatrix[idxN][idxM] = diffVector.norm() / (sgm * (pMDW1[idxN][0] + pMDW2[idxM][0]));
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
			return sigmoid(sepFactor, penalty, 100, 1.5);
		}
		return sepFactor;
	}

double AHD(std::vector<double>& pX, std::vector<double>& pY) {
	double AHD = 0;
	for (size_t i = 1; i < pX.size(); ++i) {
		
		AHD += sqrt((pX[i] - pX[i - 1])*(pX[i] - pX[i - 1]) + (pY[i] - pY[i - 1])*(pY[i] - pY[i - 1]));
	}
	return AHD;
}

double dls(Eigen::Vector3d& tangent1, Eigen::Vector3d& tangent2) {
	if (tangent1.norm() < EPSILON or tangent2.norm() < EPSILON) {
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

double DDI(std::vector<TrajectoryTemplate*>& well,const std::vector<Eigen::Vector3d>& pCartesian, bool actFunc, double penalty) {
	std::vector<Eigen::Vector4d> pmd = allPointsMD(well);
	double tortuos = Tortuosity(well); 
	double DDI;
	std::vector<double> pX, pY;
	for (size_t idx = 0; idx < pCartesian.size(); ++idx) {
		pX.push_back(pCartesian[idx][0]);
		pY.push_back(pCartesian[idx][1]);
	}
	double ahd = AHD(pX,pY);
	double MD = allLength(well);
	double TVD = pCartesian.back()[2] - pCartesian[0][2] ;
	DDI = log10((1. / 0.305) * ahd * MD * tortuos / TVD);
	if (actFunc) {
		return sigmoid(DDI, penalty, -2.5, 6.25);
	}
	return DDI;
}	

double orderScore1(std::vector<TrajectoryTemplate*>& mainWell, std::vector<std::vector<Eigen::Vector3d>>& pCTrajectories, 
	std::vector<std::vector<Eigen::Vector4d>>& pMDTrajectories,double SepFactorShift, double penalty) {
	double mainLength = 0., mainDDI = 0., rSepFactor = 0.;
	int condition = solve(mainWell);
	if (condition != 0) {
		return 10*penalty *(-condition) / mainWell.size(); // penalty * percent of incorrect templates.
	}
	mainLength = allLength(mainWell); 
	std::vector<Eigen::Vector3d> mainPCartesian = allPointsCartesian(mainWell);
	std::vector<Eigen::Vector4d> mainPMD = allPointsMD(mainWell);
	double IdealLength;
	mainWell[0]->getInitPoint();
	mainWell.back()->getTarget1Point();
	mainWell.back()->getTarget3Point();
	IdealLength = (mainWell.back()->pointT1 - mainWell[0]->pointInitial).norm() + (mainWell.back()->pointT3 - mainWell.back()->pointT1).norm();
	if (IdealLength == 0)
		return 0.;
	mainDDI = DDI(mainWell,mainPCartesian);
	if (pCTrajectories.size() == 0) {
		return mainLength/IdealLength + mainDDI;
	}
	for (size_t idx = 0; idx < pCTrajectories.size(); ++idx) {
		rSepFactor += std::max(sepFactor(mainPCartesian, mainPMD, pCTrajectories[idx],
			pMDTrajectories[idx], SepFactorShift),
			sepFactor(pCTrajectories[idx], pMDTrajectories[idx],
				mainPCartesian, mainPMD,SepFactorShift));
	}
	return mainLength / IdealLength + rSepFactor + mainDDI;
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