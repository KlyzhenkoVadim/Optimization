#include "Eigen/Dense"
#include <iostream>
#include <vector>
#include "TrajectoryTemplate.h"
#include "CurveHold.h"
#include "CurveHoldCurveHold.h"
#include "Hold.h"
#include <fstream>
#include "CostFuncs.h"
#include "PSO.h"
#include <cmath>
#include "API.h"

void writeDataCartesian(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename);
void writeDataMD(std::vector<Eigen::Vector4d>& pointsMD, std::string filename);
double Rastrigin(Eigen::VectorXd x) {
	size_t n = x.size();
	double result = 10. * n;
	for (size_t i = 0; i < n; ++i)
		result += x(i) * x(i) - 10 * cos(2 * PI * x(i));
	return result;

}
double Sphere(Eigen::VectorXd x) {
	size_t n = x.size();
	double result = 0.;
	for (size_t idx = 0; idx < n; ++idx)
		result += x(idx) * x(idx);
	return result;
}
void getOptData(PSOvalueType op);
void writeDataOpt(std::vector<size_t> order, size_t wellNum, PSOvalueType Opt, std::vector<std::vector<TrajectoryTemplate*>>& trajs);
void writeDataOptSep(std::vector<std::vector<Eigen::Vector3d>>& pCTrajs, std::vector<std::vector<Eigen::Vector4d>>& pMDTrajs,size_t& num);
std::pair<double,double> sqrtVar(std::vector<double> v) {
	double mean=0, var=0;
	size_t n = v.size();
	for (size_t i = 0; i < n; ++i) {
		mean += v[i] / n;
	}
	for (size_t i = 0; i < n; ++i) {
		var += (v[i] - mean) * (v[i] - mean)/n;
	}
	return { mean,sqrt(var) };
}
void normTransform(Eigen::VectorXd& x) {
	x[0] *= 1000;
	x[1]*= 1000;
	x[2]*= 1000;
	x[3] *= 180. / PI;
	x[4] *= 180. / PI;
	x[5]*= 1000;
}

int main()
{
	double m = 1800 / PI;
	// x = x,y,z,inc,azi,beta
	std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd& x)> Well1 = [&](Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> well1;
		Eigen::Vector3d pIHold = { 0,0,0 };
		Eigen::Vector3d pTHold = { 0,0,x(5) };
		Eigen::Vector3d pT3Chch1 = { x(0),x(1),x(2) };
		Eigen::Vector3d Tangent = calcTangentVector(x[4], x[3]);
		Eigen::Vector3d pT1Chch1 = pT3Chch1 - 110. * Tangent;
		Eigen::Vector3d pT1Chch2 = { 400.0,8.0,3000.0 };
		Eigen::Vector3d pT3Chch2 = { 1500. ,8.0 ,3010.0 };
		well1.push_back(new Hold(pIHold, pTHold));
		well1.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pT1Chch1, pT3Chch1));
		well1.push_back(new CurveHoldCurveHold(pT3Chch1, x[3], x[4], m, m, pT1Chch2, pT3Chch2));
		return well1;
	};
	std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd& x)> Well2 = [&](Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> well2;
		Eigen::Vector3d pIHold = { 0,15,0 };
		Eigen::Vector3d pTHold = { 0,15,x[5] };
		Eigen::Vector3d pT3Chch1 = { x[0],x[1],x[2] };
		Eigen::Vector3d Tangent = calcTangentVector(x[4], x[3]);
		Eigen::Vector3d pT1Chch1 = pT3Chch1 - 110. * Tangent;
		Eigen::Vector3d pT1Chch2 = { -700.,8.0,3000. };
		Eigen::Vector3d pT3Chch2 = { -1500.,8,3010. };
		well2.push_back(new Hold(pIHold, pTHold));
		well2.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pT1Chch1, pT3Chch1));
		well2.push_back(new CurveHoldCurveHold(pT3Chch1, x[3], x[4], m, m, pT1Chch2, pT3Chch2));
		return well2;
	};
	std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd& x)> Well3 = [&](Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> well3;
		Eigen::Vector3d pIHold = { 0,20,0 };
		Eigen::Vector3d pTHold = { 0,20,x[5] };
		Eigen::Vector3d pT3Chch1 = { x[0],x[1],x[2] };
		Eigen::Vector3d Tangent = calcTangentVector(x[4], x[3]);
		Eigen::Vector3d pT1Chch1 = pT3Chch1 - 110. * Tangent;
		Eigen::Vector3d pT1Chch2{ -650.,420.,3000. };
		Eigen::Vector3d pT3Chch2 = { -1500.,420.,3010. };
		well3.push_back(new Hold(pIHold, pTHold));
		well3.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pT1Chch1, pT3Chch1));
		well3.push_back(new CurveHoldCurveHold(pT3Chch1, x[3], x[4], m, m, pT1Chch2, pT3Chch2));
		return well3;
	};

	std::vector < std::function<std::vector<TrajectoryTemplate*>(Eigen::VectorXd&)>> mainWell = { Well1,Well2,Well3 };
	std::vector < std::vector<TrajectoryTemplate*>> trajectories;
	std::vector<std::vector<Eigen::Vector3d>> pCTrajectories, tmpC;
	std::vector<std::vector<Eigen::Vector4d>> pMDTrajectories, tmpMD;
	std::vector<double> inert(300, .9);
	std::vector<double> minValues = { -1000., -1000., 100,0., -180., 100. };//{ -1000., -1000., 2840, 60., -180., 100. };
	std::vector<double> maxValues = { 1000., 1000., 3000., 90., 180., 1000. };//{ 1000., 1000., 2955., 70., 180., 1000. };
	//std::vector<double> minValues = { -1., -1., 2.840, 60./180.*PI, -PI, .1 };
	//std::vector<double> maxValues = { 1., 1., 2.955, 70./180.*PI,PI, 1. };
	double lambda = 1e-2;
	std::vector<Eigen::VectorXd> args;
	std::vector<double> lens;
	size_t  N = 3, collisionNum = 0;
	bool opt_ON = true;
	std::vector<std::vector<size_t>> order = { {1,2,3},{1,3,2},{2,1,3},{2,3,1},{3,1,2},{3,2,1} };
	if (opt_ON) {
		for (int i = 0; i < 5; ++i) {
			for (size_t j = 0; j < 3; ++j) {
				std::function<double(Eigen::VectorXd)> scoreOne = [&](Eigen::VectorXd x) {
					std::vector<TrajectoryTemplate*> tmp = mainWell[j](x);
					return OneWellScore(tmp);
				};
				PSOvalueType optOne = PSO(scoreOne, minValues, maxValues, 30, 6, inert, 0.3, 0.5, 150);
				getOptData(optOne);
				args.push_back(optOne.first);
				if (j == 0) {
				trajectories.push_back(mainWell[0](args[0])); // order[i][0] - 1;
				int tmpCond = solve(trajectories.back());
				if (tmpCond > 0) {
					std::cout << "Unbuilt Trajectory" << std::endl;
				}
				writeDataOpt(order[0], order[0][0], optOne, trajectories);
				pCTrajectories.push_back(allPointsCartesian(trajectories.back()));
				pMDTrajectories.push_back(allPointsMD(trajectories.back()));
			}
				
			}
			
			for (size_t id = 1; id < 3; ++id) {
				size_t idx = id;
				//idx = order[0][id] - 1;
				std::function<double(Eigen::VectorXd)> score = [&](Eigen::VectorXd x) {
					double regularize = lambda * (x - args[idx]).lpNorm<2>();
					std::vector<TrajectoryTemplate*> tmp = mainWell[idx](x);
					return orderScore1(tmp, pCTrajectories, pMDTrajectories) + regularize; };
				PSOvalueType opt = PSO(score, minValues, maxValues, 30, 6, inert, 0.3, 0.5, 150);
				trajectories.push_back(mainWell[idx](opt.first));
				int cond = solve(trajectories.back());
				getOptData(opt);
				if (cond > 0) {
					std::cout << "Unbuild Trajectory\n";
					break;
				}
				writeDataOpt(order[0], order[0][id], opt, trajectories);
				pCTrajectories.push_back(allPointsCartesian(trajectories.back()));
				pMDTrajectories.push_back(allPointsMD(trajectories.back()));
			}
			pCTrajectories.clear();
			pMDTrajectories.clear();
			trajectories.clear();
		}
	}

	Point2d Init{0,0};
	GeoPoint Targets{ 400.0,8.0,3000.0 , 1500. ,8.0 ,3010.0 };
	GeoPoint Targets2{ -700.,8.0,3000. , -1500.,8,3010. };
	Solver S;
	//S.setData(Init, Targets);
	//S.Optimize();
	//double Length = S.getTrajectoryLength();
	//getOptData(S.getPSOdata());
	//std::cout << "Length of first well is: " << Length << std::endl;
	/*std::function<double(Eigen::VectorXd)> scoreOne = [&](Eigen::VectorXd x) {
		std::vector<TrajectoryTemplate*> tt = Well1(x);
		return OneWellScore(tt);
	};
	size_t numsIters = 100;
	std::vector<double> mins{ -1000,-1000,100,0,-180,100 }, max{ 1000,1000,3000,90,180,1000 },inertIa(numsIters,0.9);
	PSOvalueType opt1 = PSO(scoreOne, mins, max, 30, 6, inertIa, 0.3, 0.5, numsIters);
	getOptData(opt1);
	std::vector<TrajectoryTemplate*> tmpWells = Well1(opt1.first);
	int c = solve(tmpWells);
	std::cout << allLength(tmpWells);*/
	return 0;
}

void getOptData(PSOvalueType op) {
	std::cout << "gBestCost: " << op.second << std::endl;
	std::cout << "gBestPos: ";
	for (auto x : op.first)
		 std::cout<< x << ", ";
	std::cout << std::endl;
}

void writeDataCartesian(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename) {
	std::ofstream output;
	output.open(filename);
	for (size_t i = 0; i < pointsCartesian.size(); ++i) {
		output << pointsCartesian[i][0] << "," << pointsCartesian[i][1] << "," << pointsCartesian[i][2] << std::endl;
	}
	output.close();
}

void writeDataMD(std::vector<Eigen::Vector4d>& pointsMD, std::string filename) {
	std::ofstream output;
	output.open(filename);
	for (size_t i = 0; i < pointsMD.size(); ++i) {
		output << pointsMD[i][0] << "," << pointsMD[i][1] << "," << pointsMD[i][2] <<"," << pointsMD[i][3] << std::endl;
	}
	output.close();
}

void writeDataOpt(std::vector<size_t> order,size_t wellNum, PSOvalueType Opt, std::vector<std::vector<TrajectoryTemplate*>>& trajs) {
	std::fstream file;
	file.open("C:/Users/klyzhenko.vs/0CET/optDataCpp1.csv", std::ios::out | std::ios::app); // C:/Users/HP/WellboreOpt/optDataCpp2.csv
	file << "[";
	for (size_t i = 0; i < order.size(); ++i) {
		file << order[i];
		if(i!=order.size()-1)
			file << " ";
	}
	file << "]," << wellNum << ",[";
	for (size_t i = 0; i < Opt.first.size(); ++i) {
		file << Opt.first[i];
		if (i != Opt.first.size() - 1)
			file << " ";
	}
	file << "]," << Opt.second << ",";
	file << allLength(trajs.back()) << ",";
	std::vector<Eigen::Vector3d> pC = allPointsCartesian(trajs.back());
	std::vector<Eigen::Vector4d> pMD = allPointsMD(trajs.back());
	file << DDI(pC, pMD, false)<<",";
	size_t n = trajs.size();
	if (n == 1)
		file << "-\n";
	if (n == 2) {
		std::vector<Eigen::Vector3d>pC1 = allPointsCartesian(trajs[0]);
		std::vector<Eigen::Vector4d>pMD1 = allPointsMD(trajs[0]);
		file << sepFactor(pC, pMD, pC1, pMD1, false) << std::endl;
	}
	if (n == 3) {
		std::vector<Eigen::Vector3d>pC1 = allPointsCartesian(trajs[0]), pC2 = allPointsCartesian(trajs[1]);
		std::vector<Eigen::Vector4d>pMD1 = allPointsMD(trajs[0]), pMD2 = allPointsMD(trajs[1]);
		double sf = std::min(sepFactor(pC, pMD, pC1, pMD1, false), sepFactor(pC, pMD, pC2, pMD2, false));
		file << sf << std::endl;
	}
	file.close();
}

void writeDataOptSep(std::vector<std::vector<Eigen::Vector3d>>& pCTrajs, std::vector<std::vector<Eigen::Vector4d>>& pMDTrajs,size_t & num){
	double sf12 = std::min(sepFactor(pCTrajs[0], pMDTrajs[0], pCTrajs[1], pMDTrajs[1], false), sepFactor(pCTrajs[1], pMDTrajs[1], pCTrajs[0], pMDTrajs[0], false));
	double sf13 = std::min(sepFactor(pCTrajs[0], pMDTrajs[0], pCTrajs[2], pMDTrajs[2], false), sepFactor(pCTrajs[2], pMDTrajs[2], pCTrajs[0], pMDTrajs[0], false));
	double sf23 = std::min(sepFactor(pCTrajs[1], pMDTrajs[1], pCTrajs[2], pMDTrajs[2], false), sepFactor(pCTrajs[2], pMDTrajs[2], pCTrajs[1], pMDTrajs[1], false));
	if (sf12 - 1.5 < EPSILON or sf13 - 1.5 < EPSILON or sf23 - 1.5 < EPSILON) {
		num += 1;
	}
	std::cout << "SepF 1-2: " << sf12 << std::endl;
	std::cout << "SepF 1-3: " << sf13 << std::endl;
	std::cout << "SepF 2-3: " << sf23<< std::endl;
}