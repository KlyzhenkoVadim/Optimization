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

PSOvalueType PSO(std::function<double(const std::vector<double>&)> func, const std::vector<double>& minValues, const std::vector<double>& maxValues,
	size_t numAgents, size_t dimension, size_t numIterations, const std::vector<double>& inertia, double socCoef, double indCoef, bool history)
{
	bool flag = false;
	std::vector<double> vMax;
	vMax.reserve(dimension);
	for (size_t idx = 0; idx < dimension; ++idx) {
		vMax.push_back(fabs(maxValues[idx] - minValues[idx]) / 2);
	}
	std::vector<std::vector<double>> Xposition;
	std::vector<std::vector<double>> Velocities;
	std::vector<double> gBestPos;
	std::vector<std::vector<double>> pBestPos;

	Xposition.resize(numAgents);
	Velocities.resize(numAgents);
	gBestPos.resize(dimension);
	pBestPos.resize(numAgents);

	//std::random_device rd;
	std::mt19937 gen(time(NULL));
	std::uniform_real_distribution<double> dist;
	double minFunc = 1e3;
	size_t idxMin;
	for (size_t i = 0; i < numAgents; ++i) {
		for (size_t j = 0; j < dimension; ++j) {
			Xposition[i].push_back(minValues[j] + (maxValues[j] - minValues[j]) * dist(gen));
			Velocities[i].push_back(0.);
		}
	}

	pBestPos = Xposition;
	idxMin = 0;
	minFunc = func(Xposition[0]);
	for (size_t k = 1; k < numAgents; ++k) {
		double tmp = func(Xposition[k]);
		if (tmp - minFunc < pso::EPSILON) {
			minFunc = tmp;
			idxMin = k;
		}
	}
	gBestPos = Xposition[idxMin];
	PSOvalueType result;
	if (history)
	{
		result.costHist.reserve(numIterations + 1);
		result.costHist.push_back(minFunc);
	}
	//writegBestPos(gBestPos,"bestPosEvo"+numtarget+".txt");
	//writeCurrbestCost(minFunc, numtarget+".txt");
	double r1, r2;
	double R = 0.5;
	for (size_t iteration = 0; iteration < numIterations; ++iteration)
	{
		// r1 r2
		for (size_t idAg = 0; idAg < numAgents; ++idAg)
		{
			r1 = dist(gen);
			r2 = dist(gen);
			for (size_t idim = 0; idim < dimension; ++idim) {
				auto& v = Velocities[idAg][idim];
				auto& x = Xposition[idAg][idim];
				
				if (idAg == idxMin)
				{
					v = inertia[iteration] * v + R * (1. - 2. * dist(gen));
				}
				else {
					v = (inertia[iteration] * v +
						r1 * indCoef * (pBestPos[idAg][idim] - x) +
						r2 * socCoef * (gBestPos[idim] - x));
				}
				
				if (v - vMax[idim] > pso::EPSILON)
					v = vMax[idim];
				else if (v + vMax[idim] < -pso::EPSILON)
					v = -vMax[idim];
				if ((x + v - maxValues[idim] < pso::EPSILON) &&
					(x + v - minValues[idim] > pso::EPSILON))
					x += v;
				else if ((x - v - maxValues[idim] < pso::EPSILON) &&
					(x - v - minValues[idim] > pso::EPSILON))
					x -= v;
			}
		}
		double tmpf = minFunc;
		for (size_t idxAg = 0; idxAg < numAgents; ++idxAg) {
			double tmpfuncMean = func(Xposition[idxAg]);
			if (tmpfuncMean - func(pBestPos[idxAg]) < pso::EPSILON)
				pBestPos[idxAg] = Xposition[idxAg];
			if (tmpfuncMean - minFunc < pso::EPSILON) {
				idxMin = idxAg;
				gBestPos = Xposition[idxAg];
				minFunc = func(gBestPos);
				//writegBestPos(gBestPos, "bestPosEvo" + numtarget + ".txt");
				//writeCurrbestCost(minFunc, numtarget + ".txt");
			}
		}
		//writeCurrbestCost(minFunc, numtarget + ".txt");
		if (history)
		{
			result.costHist.push_back(minFunc);
		}		
	}
	result.cost = minFunc;
	result.optVec = gBestPos;
	
	return result;
}

void writePSO(double func, size_t iter, std::string filename) {
	std::fstream file;
	file.open(filename, std::ios::out | std::ios::app);
	file << iter << "," << func << std::endl;
	file.close();
}