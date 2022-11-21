#include "PSO.h"
#include <fstream>
void writePSO(double func, size_t iter, std::string filename);

void writegBestPos(const Eigen::VectorXd& x, std::string filename)
{
	std::fstream file;
	file.open(filename, std::ios::out | std::ios::app);
	for (size_t i = 0; i < x.size();++i)
	{
		file << x[i]<<",";
	}
	file << std::endl;
	file.close();
}

void writeCurrbestCost(double c, std::string filename)
{
	std::fstream file;
	file.open(filename, std::ios::out | std::ios::app);
	file << c << std::endl;
	file.close();
}

PSOvalueType PSO(std::function<double(const Eigen::VectorXd&)> func, const std::vector<double>& minValues, const std::vector<double>& maxValues,
	size_t numAgents, size_t dimension, size_t numIterations, const std::vector<double>& inertia, double socCoef, double indCoef, const std::string& numtarget){
	bool flag = false;
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
	std::uniform_real_distribution<double> dist;
	double minFunc = 1e3;
	size_t idxgBest;
		for (size_t i = 0; i < numAgents; ++i) {
			for (size_t j = 0; j < dimension; ++j) {
				Xposition[i](j) = minValues[j] + (maxValues[j] - minValues[j]) * dist(gen);
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
	//writegBestPos(gBestPos,"bestPosEvo"+numtarget+".txt");
	//writeCurrbestCost(minFunc, "bestCurrEvo"+numtarget+".txt");
	double r1, r2;
	for (size_t iteration = 0; iteration < numIterations; ++iteration) {
		// r1 r2
		for (size_t idAg = 0; idAg < numAgents; ++idAg) {
			r1 = dist(gen);
			r2 = dist(gen);
			for (size_t idim = 0; idim < dimension; ++idim) {
				Velocities[idAg](idim) = (inertia[iteration] * Velocities[idAg](idim) +
					r1 * indCoef * (pBestPos[idAg](idim) - Xposition[idAg](idim)) +
					r2 * socCoef * (gBestPos[idim] - Xposition[idAg](idim)));
				if (Velocities[idAg](idim) - vMax[idim] > 1e-10)
					Velocities[idAg](idim) = vMax[idim];
				else if (Velocities[idAg](idim) + vMax[idim] < -1e-10)
					Velocities[idAg](idim) = -vMax[idim];
				if ((Xposition[idAg](idim) + Velocities[idAg](idim) - maxValues[idim] < 1e-10) and
					(Xposition[idAg](idim) + Velocities[idAg](idim) - minValues[idim] > 1e-10))
					Xposition[idAg](idim) += Velocities[idAg](idim);
				else if((Xposition[idAg](idim) - Velocities[idAg](idim) - maxValues[idim] < 1e-10) and
					(Xposition[idAg](idim) - Velocities[idAg](idim) - minValues[idim] > 1e-10))
					Xposition[idAg](idim) -= Velocities[idAg](idim);
			}
		}
		double tmpf = minFunc;
		for (size_t idxAg = 0; idxAg < numAgents; ++idxAg) {
			double tmpfuncMean = func(Xposition[idxAg]);
			if (tmpfuncMean - func(pBestPos[idxAg]) < 1e-10)
				pBestPos[idxAg] = Xposition[idxAg];
			if (tmpfuncMean - minFunc < 1e-10) {
				gBestPos = Xposition[idxAg];
				minFunc = func(gBestPos);
				//writegBestPos(gBestPos, "bestPosEvo" + numtarget + ".txt");
				//writeCurrbestCost(minFunc, "bestCurrEvo" + numtarget + ".txt");
			}
		}
	}
	return {gBestPos,minFunc};
}

void writePSO(double func, size_t iter, std::string filename) {
	std::fstream file;
	file.open(filename, std::ios::out | std::ios::app);
	file << iter << "," << func << std::endl;
	file.close();
}