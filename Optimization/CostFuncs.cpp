#include "CostFuncs.h"
#include <cmath>

double sigmoid(double x, double penalty, double alpha, double x0) {
	return penalty / (1 + exp(alpha * (x - x0)));
}

double signum(double x) {
	if (x > 0.) {
		return +1;
	}
	if (x < 0.) {
		return -1;
	}
	if (x < 1e-5) {
		return 0;
	}
}

double sepFactor(std::vector<Eigen::Vector3d>& pCartesianW1,
	std::vector<Eigen::Vector4d>& pMDW1,
	std::vector<Eigen::Vector3d>& pCartesianW2, bool actFunc = true, double penalty = 100) {
	size_t n = pCartesianW1.size(), m = pCartesianW2.size();
	std::vector<std::vector<double>> distanceMatrix;
	std::vector<std::vector<double>> scalarProdMatrix;
	double sgm = 4 * 5e-3;
	double eps = 1e-6;
	for (size_t i = 0; i < n; ++i) {
		scalarProdMatrix[i][0] = -1.;
		for (size_t j = 0; j < m; ++j) {
			distanceMatrix[i][j] = 1e6;
		}
	}	
	size_t k = 2;
	for (size_t idxN = 1; idxN < n; ++idxN) {
		for (size_t idxM = k - 1; idxM < m; ++idxM) {
			Eigen::Vector3d diffVector = pCartesianW2[idxM] - pCartesianW1[idxN];
			Eigen::Vector3d tangentW1 = {pMDW1[idxN][1], pMDW1[idxN][2], pMDW1[idxN][3]};
			if (tangentW1.norm() < eps or diffVector.norm() < eps) {
				scalarProdMatrix[idxN][idxM] = 0.;
			}
			else {
				scalarProdMatrix[idxN][idxM] = round(tangentW1.dot(diffVector) / tangentW1.norm() / diffVector.norm()*100000)/100000;
			}
			if (signum(scalarProdMatrix[idxN][idxM]) != signum(scalarProdMatrix[idxN][idxM - 1])) {
				distanceMatrix[idxN][idxM] = diffVector.norm() / sgm / pMDW1[idxN][0];
			}
		}
		double sepFactor = 1e3;
		for (size_t i = 0; i < n; i++){
			for (size_t j = 0; j < m; ++j) {
				if (sepFactor > distanceMatrix[i][j]) {
					sepFactor = distanceMatrix[i][j];
				}
				
			}
		}
		if (actFunc) {
			return sigmoid(sepFactor, penalty, 10, 1.5);
		}
		return sepFactor;
	}
	
	// ...
}

double AHD(std::vector<double>& pX, std::vector<double>& pY) {
	double AHD = 0;
	if (pX.size() != pY.size()) {
		throw(std::exception("length pX and pY must be equal"));
	}
	for (size_t i = 1; i < pX.size(); ++i) {
		
		AHD += sqrt((pX[i] - pX[i - 1])*(pX[i] - pX[i - 1]) + (pY[i] - pY[i - 1])*(pY[i] - pY[i - 1]));
	}
	return AHD;
}

double dls(Eigen::Vector3d& tangent1, Eigen::Vector3d& tangent2) {
	if (tangent1.norm() < EPSILON or tangent2.norm() < EPSILON) {
		return 0.;
	}
	double angle = tangent1.dot(tangent2) / tangent1.norm() / tangent2.norm();
	if (angle + 1 < EPSILON or angle - 1 > EPSILON) {
		return 0.;
	}
	return (180 / PI) * acos(angle);
}

double Tortuosity(std::vector<Eigen::Vector4d>& pointsMD) {
	double tortuos = 0;
	for (size_t i = 1; i < pointsMD.size(); ++i) {
		Eigen::Vector3d tPresent = { pointsMD[i][1],pointsMD[i][2],pointsMD[i][3]};
		Eigen::Vector3d tPrev = { pointsMD[i-1][1],pointsMD[i-1][2],pointsMD[i-1][3]};
		tortuos += dls(tPresent, tPrev);
	}
	return tortuos;
}

double DDI(std::vector<Eigen::Vector3d>& pCartesian,std::vector<Eigen::Vector4d>& pMD, bool actFunc = true, double penalty = 10.) {
	double toruos = Tortuosity(pMD);
	double DDI;
	std::vector<double> pX, pY;
	for (size_t idx = 0; idx < pCartesian.size(); ++idx) {
		pX[idx] = pCartesian[idx][0];
		pY[idx] = pCartesian[idx][1];
	}
	double ahd = AHD(pX,pY);
	double MD = pMD.back()[0] - pMD[0][0];
	double TVD = pCartesian.back()[2] - pCartesian[0][2] ;
	DDI = log10((1. / 0.305) * ahd * MD * toruos / TVD);
	if (actFunc) {
		return sigmoid(DDI, penalty, -30, 6.4);
	}
}	


double orderScore(std::vector<TrajectoryTemplate*>& mainWell, std::vector<std::vector<TrajectoryTemplate*>>& Trajectories, double penalty = 1000) {
	double mainLength, mainDDI,rSepFactor = 0;
	try {
		mainLength = allLength(mainWell); // здесь нужен fit,  a в роли аргументов должен быть динам. массив...
	}
	catch (std::runtime_error err) {
		return penalty;
	}
	Eigen::Vector3d tmpVec = mainWell.back()->pointsCartesian.back() - mainWell[0]->pointsCartesian[0];
	std::vector<Eigen::Vector3d> mainPCartesian = allPointsCartesian(mainWell);
	std::vector<Eigen::Vector4d> mainPMD = allPointsMD(mainWell);
	mainLength /= tmpVec.norm();
	mainDDI = DDI(mainPCartesian, mainPMD);
	if (Trajectories.size() == 0) {
		return mainLength + mainDDI;
	}
	for (size_t idx = 0; idx < Trajectories.size(); ++idx) {
		std::vector<Eigen::Vector3d> pCartesianTrajectory = allPointsCartesian(Trajectories[idx]);
		rSepFactor += sepFactor(mainPCartesian, mainPMD, pCartesianTrajectory);
	}
	return mainLength + mainDDI + rSepFactor;

}

