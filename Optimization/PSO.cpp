#include "PSO.h"

PSOvalueType PSO(std::function<double(Eigen::VectorXd)> func, std::vector<double>& minValues, std::vector<double>& maxValues,
	size_t numAgents, size_t dimension, double socCoef, double indCoef,size_t numIterations , std::vector<double>& inertia) {
	
	assert ("Number of bound conditions is not equal", minValues.size() != maxValues.size());
	
	bool flag = false;
	for (size_t i = 0; i < minValues.size();++i) {
		assert("MinValue more than MaxValues", minValues[i] - maxValues[i] > 1e-10);
	}
	std::vector<double> vMax;
	for (size_t idx = 0; idx < dimension; ++idx) {
		vMax.push_back(fabs(maxValues[idx] - minValues[idx]) / 2);
	}
	Eigen::MatrixXd Xposition(numAgents,dimension);
	Eigen::MatrixXd Velocities(numAgents,dimension);
	Eigen::VectorXd gBestPos(dimension);
	Eigen::MatrixXd pBestPos(numAgents,dimension);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution distCoefs(0., 1.);

	for (size_t i = 0; i < dimension; ++i) {
		std::uniform_real_distribution distPos(minValues[i], maxValues[i]);
		std::uniform_real_distribution distVel(-vMax[i], vMax[i]);
		for (size_t j = 0; j < numAgents; ++j) {
			Xposition(j,i) = distPos(gen);
			Velocities(j,i) = distVel(gen);
		}
	}

	pBestPos = Xposition;
	size_t idxgBest = 0;
	double minFunc = func(Xposition.row(0));
	for (size_t k = 1; k < numAgents; ++k) {
		double tmp = func(Xposition.row(k));
		if (tmp - minFunc < 1e-10){
			minFunc = tmp;
			idxgBest = k;
		}
	}
	gBestPos = Xposition.row(idxgBest);

	for (size_t iteration = 0; iteration < numIterations; ++iteration) {
		// r1 r2 
		for (size_t idAg = 0; idAg < numAgents; ++idAg) {
			double r1 = distCoefs(gen), r2 = distCoefs(gen);
			for (size_t idim = 0; idim < dimension; ++idim) {
				Velocities(idAg, idim) = (inertia[iteration] * Velocities(idAg, idim) +
					r1 * indCoef * (pBestPos(idAg, idim) - Xposition(idAg, idim)) +
					r2 * socCoef * (gBestPos[idim] - Xposition(idAg, idim)));
				if (Velocities(idAg, idim) - vMax[idim] > 1e-10)
					Velocities(idAg, idim) = vMax[idim];
				else if (Velocities(idAg, idim) + vMax[idim] < -1e-10)
					Velocities(idAg, idim) = -vMax[idim];
				if ((Xposition(idAg, idim) + Velocities(idAg, idim) - maxValues[idim] <= 1e-10) and
					(Xposition(idAg, idim) + Velocities(idAg, idim) - minValues[idim] >= 1e-10))
					Xposition(idAg, idim) += Velocities(idAg, idim);
				else if((Xposition(idAg, idim) - Velocities(idAg, idim) - maxValues[idim] <= 1e-10) and
					(Xposition(idAg, idim) - Velocities(idAg, idim) - minValues[idim] >= 1e-10))
					Xposition(idAg, idim) -= Velocities(idAg, idim);
			}
		}
		for (size_t idxAg = 0; idxAg < numAgents; ++ idxAg){
			double tmpfuncMean = func(Xposition.row(idxAg)); 
			if (tmpfuncMean - func(pBestPos.row(idxAg)) < 1e-10)
				pBestPos.row(idxAg) = Xposition.row(idxAg);
			if (tmpfuncMean - func(gBestPos) < 1e-10)
				gBestPos = Xposition.row(idxAg);
		}
	}
	return { gBestPos, func(gBestPos) };

}