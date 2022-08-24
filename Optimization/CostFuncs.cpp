#include "CostFuncs.h"
#include <cmath>
#include <iostream>
double sigmoid(double x, double penalty, double alpha, double x0) {
	return penalty / (1 + exp(alpha * (x - x0)));
}

double stair(double x, double penalty, double x0) {
	if (x - x0 <= EPSILON)
		return penalty;
	return 0.;
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
	std::vector<Eigen::Vector3d>& pCartesianW2,std::vector<Eigen::Vector4d>& pMDW2, bool actFunc, double penalty) {
	size_t n = pCartesianW1.size(), m = pCartesianW2.size();
	Eigen::MatrixXd distanceMatrix(n, m);//std::vector<std::vector<double>> distanceMatrix;
	Eigen::MatrixXd scalarProdMatrix(n, m); //std::vector<std::vector<double>> scalarProdMatrix;
	double sgm = 2e-3;
	double eps = 1e-6;
	bool flg = true;
	for (size_t i = 0; i < n; ++i) {
		scalarProdMatrix(i,0) = -1.;
		for (size_t j = 0; j < m; ++j) {
			distanceMatrix(i,j) = 1e3;
		}
	}	
	size_t k = 1;
	for (size_t idxN = 0; idxN < n; ++idxN) {
		for (size_t idxM = k - 1; idxM < m; ++idxM) {
			Eigen::Vector3d diffVector = pCartesianW2[idxM] - pCartesianW1[idxN];
			Eigen::Vector3d tangentW1 = { pMDW1[idxN][1], pMDW1[idxN][2], pMDW1[idxN][3] };
			if (diffVector.norm() < eps) {
				scalarProdMatrix(idxN, idxM) = 0;
				distanceMatrix(idxN, idxM) = 0.;
			}
			else {
				scalarProdMatrix(idxN, idxM) = round(tangentW1.dot(diffVector) / tangentW1.norm() / diffVector.norm() * 100000) / 100000;
				if (idxM > 1 and (signum(scalarProdMatrix(idxN, idxM)) != signum(scalarProdMatrix(idxN, idxM - 1)))) {
					if (pMDW1[idxN][0] + pMDW2[idxM][0] > eps)
						distanceMatrix(idxN, idxM) = diffVector.norm() / (sgm * (pMDW1[idxN][0] + pMDW2[idxM][0]));
					if (flg) {
						k = idxM;
						flg = false;
					}
				}
			}
		}
		flg = true;
	}
		double sepFactor = 1e3;
		for (size_t i = 0; i < n; ++i){
			for (size_t j = 0; j < m; ++j) {
				sepFactor = std::min(distanceMatrix(i, j), sepFactor);
			}
		}

		if (actFunc) {
			return sigmoid(sepFactor, penalty, 100, 1.5);
		}
		return sepFactor;
	}

double AHD(std::vector<double>& pX, std::vector<double>& pY) {
	double AHD = 0;
	assert("length pX and pY must be equal", pX.size() != pY.size());
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
	if (dotProd + 1 < EPSILON or dotProd - 1 > EPSILON) {
		return 0.;
	}
	if (fabs(dotProd - 1.) < EPSILON)
		dotProd = 1.;
	if (fabs(dotProd + 1) < EPSILON)
		dotProd = -1.;
	return (180. / PI) * acos(dotProd);
}

double Tortuosity(std::vector<Eigen::Vector4d>& pointsMD) {
	double tortuos = 0.;
	for (size_t i = 1; i < pointsMD.size(); ++i) {
		Eigen::Vector3d tPresent = { pointsMD[i][1],pointsMD[i][2],pointsMD[i][3]};
		Eigen::Vector3d tPrev = { pointsMD[i-1][1],pointsMD[i-1][2],pointsMD[i-1][3]};
		tortuos += dls(tPresent, tPrev);
	}
	return tortuos;
}

double DDI(std::vector<Eigen::Vector3d>& pCartesian,std::vector<Eigen::Vector4d>& pMD, bool actFunc, double penalty) {
	double toruos = Tortuosity(pMD);
	double DDI;
	std::vector<double> pX, pY;
	for (size_t idx = 0; idx < pCartesian.size(); ++idx) {
		pX.push_back(pCartesian[idx][0]);
		pY.push_back(pCartesian[idx][1]);
	}
	double ahd = AHD(pX,pY);
	double MD = pMD.back()[0] - pMD[0][0];
	double TVD = pCartesian.back()[2] - pCartesian[0][2] ;
	DDI = log10((1. / 0.305) * ahd * MD * toruos / TVD);
	if (actFunc) {
		return sigmoid(DDI, penalty, -2.5, 6.4);
	}
	return DDI;
}	

double OneWellScore(std::vector<TrajectoryTemplate*>& mainWell, double penalty) {
	double mainLength, mainDDI,rSepFactor = 0.;
	int condition = solve(mainWell);
	if (condition != 0) {
		return penalty * condition/mainWell.size(); // penalty * percent of incorrect templates.
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
	mainDDI = DDI(mainPCartesian, mainPMD);

	return mainLength/IdealLength+mainDDI;

}

double orderScore1(std::vector<TrajectoryTemplate*>& mainWell, std::vector<std::vector<Eigen::Vector3d>>& pCTrajectories, 
	std::vector<std::vector<Eigen::Vector4d>>& pMDTrajectories, double penalty) {
	double mainLength = 0., mainDDI = 0., rSepFactor = 0.;
	int condition = solve(mainWell);
	if (condition != 0) {
		return 10*penalty * condition / mainWell.size(); // penalty * percent of incorrect templates.
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
	mainDDI = DDI(mainPCartesian, mainPMD);
	if (pCTrajectories.size() == 0) {
		return mainLength/IdealLength + mainDDI;
	}
	for (size_t idx = 0; idx < pCTrajectories.size(); ++idx) {
		rSepFactor += std::max(sepFactor(mainPCartesian, mainPMD, pCTrajectories[idx], pMDTrajectories[idx]),sepFactor(pCTrajectories[idx], pMDTrajectories[idx], mainPCartesian, mainPMD));
	}
	return mainLength/IdealLength + rSepFactor + mainDDI;
}
