#include "PSO.h"
#include <fstream>
void writePSO(double func, size_t iter, std::string filename);

PSOvalueType PSO(std::function<double(const Eigen::VectorXd&)> func, const std::vector<double>& minValues, const std::vector<double>& maxValues,
	size_t numAgents, size_t dimension, size_t numIterations, const std::vector<double>& inertia, double socCoef, double indCoef){
	assert ("Number of bound conditions is not equal", minValues.size() != maxValues.size());
	bool flag = false;
	for (size_t i = 0; i < minValues.size();++i) {
		assert("MinValue more than MaxValues", minValues[i] - maxValues[i] > 1e-10);
	}
	std::vector<double> vMax;
	for (size_t idx = 0; idx < dimension; ++idx) {
		vMax.push_back(fabs(maxValues[idx] - minValues[idx]) / 2);
	}
	Eigen::VectorXd tmp(dimension);
	std::vector<Eigen::VectorXd> Xposition(numAgents,tmp);//(numAgents,dimension)
	std::vector<Eigen::VectorXd> Velocities(numAgents, tmp);//(numAgents,dimension);
	Eigen::VectorXd gBestPos(dimension);
	std::vector<Eigen::VectorXd> pBestPos;//(numAgents,dimension);
	
	//std::random_device rd;
	std::mt19937 gen(time(NULL));
	std::uniform_real_distribution<double> distCoefs(0., 1.);
	std::vector<std::uniform_real_distribution<double>> distrPos;
	std::vector<std::uniform_real_distribution<double>> distrVel;
	for (size_t i = 0; i < dimension; ++i) {
		distrPos.push_back(std::uniform_real_distribution<double>(minValues[i], maxValues[i]));
		distrVel.push_back(std::uniform_real_distribution<double>(-vMax[i], vMax[i]));
	}
	double minFunc = 1e3;
	size_t idxgBest;
		for (size_t i = 0; i < numAgents; ++i) {
			for (size_t j = 0; j < dimension; ++j) {
				Xposition[i](j) = distrPos[j](gen);
				//Velocities[i](j) = distrVel[j](gen);
				Velocities[i][j] = 0.;
				
			}
		}

		pBestPos = Xposition;
		idxgBest = 0;
		minFunc = func(Xposition[0]);
		for (size_t k = 1; k < numAgents; ++k) {
			double tmp = func(Xposition[k]);
			if (tmp - minFunc < 1e-10) {
				minFunc = tmp;
				idxgBest = k;
			}
		}
	gBestPos = Xposition[idxgBest];

	for (size_t iteration = 0; iteration < numIterations; ++iteration) {
		// r1 r2
		for (size_t idAg = 0; idAg < numAgents; ++idAg) {
			double r1 = distCoefs(gen), r2 = distCoefs(gen);
			for (size_t idim = 0; idim < dimension; ++idim) {
				Velocities[idAg](idim) = (inertia[iteration] * Velocities[idAg](idim) +
					r1 * indCoef * (pBestPos[idAg](idim) - Xposition[idAg](idim)) +
					r2 * socCoef * (gBestPos[idim] - Xposition[idAg](idim)));
				if (Velocities[idAg](idim) - vMax[idim] > 1e-10)
					Velocities[idAg](idim) = vMax[idim];
				else if (Velocities[idAg](idim) + vMax[idim] < -1e-10)
					Velocities[idAg](idim) = -vMax[idim];
				if ((Xposition[idAg](idim) + Velocities[idAg](idim) - maxValues[idim] <= 1e-10) and
					(Xposition[idAg](idim) + Velocities[idAg](idim) - minValues[idim] >= 1e-10))
					Xposition[idAg](idim) += Velocities[idAg](idim);
				else if((Xposition[idAg](idim) - Velocities[idAg](idim) - maxValues[idim] <= 1e-10) and
					(Xposition[idAg](idim) - Velocities[idAg](idim) - minValues[idim] >= 1e-10))
					Xposition[idAg](idim) -= Velocities[idAg](idim);
			}
		}
		for (size_t idxAg = 0; idxAg < numAgents; ++idxAg) {
			double tmpfuncMean = func(Xposition[idxAg]);
			if (tmpfuncMean - func(pBestPos[idxAg]) < 1e-10)
				pBestPos[idxAg] = Xposition[idxAg];
			if (tmpfuncMean - minFunc < 1e-10) {
				gBestPos = Xposition[idxAg];
				minFunc = func(gBestPos);
			}
		}
		//writePSO(minFunc, iteration,"dataPSO.csv");
	}
	return {gBestPos,minFunc};

}

void writePSO(double func, size_t iter, std::string filename) {
	std::fstream file;
	file.open(filename, std::ios::out | std::ios::app);
	file << iter << "," << func << std::endl;
	file.close();
}