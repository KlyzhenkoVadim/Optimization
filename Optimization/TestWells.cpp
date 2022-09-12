#include "TestWells.h"

void getOptData(PSOvalueType op) {
	std::cout << "gBestCost: " << op.second << std::endl;
	std::cout << "gBestPos: ";
	std::cout<< std::setprecision(7);
	for (auto x : op.first)
		std::cout << x << ", ";
	std::cout << std::endl;
}

void writeDataCartesian(std::vector<Eigen::Vector3d>& pointsCartesian, std::string filename) {
	std::ofstream output;
	output.open(filename);
	output << std::setprecision(9);
	for (size_t i = 0; i < pointsCartesian.size(); ++i) {
		output << pointsCartesian[i][0] << "," << pointsCartesian[i][1] << "," << pointsCartesian[i][2] << std::endl;
	}
	output.close();
}

void writeDataMD(std::vector<Eigen::Vector4d>& pointsMD, std::string filename) {
	std::ofstream output;
	output.open(filename);
	for (size_t i = 0; i < pointsMD.size(); ++i) {
		output << pointsMD[i][0] << "," << pointsMD[i][1] << "," << pointsMD[i][2] << "," << pointsMD[i][3] << std::endl;
	}
	output.close();
}

void writeInclinometry(const std::vector<Eigen::Vector4d>& pMD, std::string filename) {
	std::ofstream output;
	output.open(filename);
	output << "md, inc, azi\n";
	for (size_t i = 0; i < pMD.size(); ++i) {
		std::pair<double, double> incAzi = CartesianToSpherical(Eigen::Vector3d{ pMD[i][1],pMD[i][2],pMD[i][3] });
		output << pMD[i][0] << "," << incAzi.first << "," << incAzi.second << std::endl;
	}
	output.close();
}

struct GeoP {
	Eigen::Vector3d pInitial, pT1, pT3;
};

void test() {
	//double m = 1800 / PI;
	//// x = x,y,z,inc,azi,beta
	//std::function<std::vector<TrajectoryTemplate*>(const Eigen::VectorXd& x, GeoP geo)> Well = [&](const Eigen::VectorXd& x, GeoP geo) {
	//	std::vector<TrajectoryTemplate*> well;
	//	Eigen::Vector3d pIHold = geo.pInitial;//{0,0,0}
	//	Eigen::Vector3d pTHold = { 0,0,x(5) };
	//	Eigen::Vector3d pT3Chch1 = { x(0),x(1),x(2) };
	//	Eigen::Vector3d Tangent = calcTangentVector(x[4], x[3]);
	//	Eigen::Vector3d pT1Chch1 = pT3Chch1 - 110. * Tangent;
	//	Eigen::Vector3d pT1Chch2 = geo.pT1;//{ 400.0,8.0,3000.0 };
	//	Eigen::Vector3d pT3Chch2 = geo.pT3;//{ 1500. ,8.0 ,3010.0 };
	//	well.push_back(new Hold(pIHold, pTHold));
	//	well.push_back(new CurveHoldCurveHold(pTHold, 0, 0, m, m, pT1Chch1, pT3Chch1));
	//	well.push_back(new CurveHoldCurveHold(pT3Chch1, x[3], x[4], m, m, pT1Chch2, pT3Chch2));
	//	return well;
	//};
	//std::vector < std::vector<TrajectoryTemplate*>> trajectories;
	//std::vector<std::vector<Eigen::Vector3d>> pCTrajectories, tmpC;
	//std::vector<std::vector<Eigen::Vector4d>> pMDTrajectories, tmpMD;
	//std::vector<GeoP> geoPointsWells{
	//	{{0,0,0,},{ 400.0,8.0,3000.0},{1500. ,8.0 ,3010.0} },
	//	{{ 0,15,0 }, { -700.,8.0,3000. },{ -1500.,8,3010. }},
	//	{{ 0,20,0 },{ -650.,420.,3000. },{ -1500.,420.,3010. }}
	//};
	//std::vector<double> inert(300, .9);
	//std::vector<double> minValues = { -1000., -1000., 100,0., -180., 100. };//{ -1000., -1000., 2840, 60., -180., 100. };
	//std::vector<double> maxValues = { 1000., 1000., 3000., 90., 180., 1000. };//{ 1000., 1000., 2955., 70., 180., 1000. };
	//double lambda = 1e-3;
	//std::vector<Eigen::VectorXd> args;
	//std::vector<double> lens;
	//size_t  N = 3, collisionNum = 0;

	//bool opt_ON = !true;
	//std::vector<std::vector<size_t>> order = { {1,2,3},{1,3,2},{2,1,3},{2,3,1},{3,1,2},{3,2,1} };
	//std::time_t start = std::time(NULL);
	//if (opt_ON) {
	//	for (int i = 0; i < 1; ++i) {
	//		{
	//			for (size_t j = 0; j < 1; ++j) {
	//				size_t idd = order[i][j] - 1;
	//				std::function<double(const Eigen::VectorXd&)> scoreOne = [&](const Eigen::VectorXd& x) {
	//					std::vector<TrajectoryTemplate*> tmp = Well(x, geoPointsWells[idd]);
	//					double oneScore = 0; // OneWellScore(tmp);
	//					for (auto x : tmp) {
	//						delete x;
	//					}
	//					return oneScore;
	//				};
	//				PSOvalueType optOne = PSO(scoreOne, minValues, maxValues, 40, 13, 150);
	//				getOptData(optOne);
	//				args.push_back(optOne.first);
	//				if (j == 0) {
	//					trajectories.push_back(Well(args[0],geoPointsWells[idd])); // order[i][0] - 1;
	//					int tmpCond = solve(trajectories.back());
	//					if (tmpCond > 0) {
	//						std::cout << "Unbuilt Trajectory" << std::endl;
	//					}
	//					//writeDataOpt(order[i], order[i][j], optOne, trajectories);
	//					pCTrajectories.push_back(allPointsCartesian(trajectories.back()));
	//					pMDTrajectories.push_back(allPointsMD(trajectories.back()));
	//					writeDataCartesian(pCTrajectories.back(), "output.txt");
	//				}
	//			}
	//		}
	//		for (size_t id = 1; id < 3; ++id) {
	//			//size_t idx = id;
	//			size_t idx = order[i][id] - 1;
	//			std::function<double(Eigen::VectorXd)> score = [&](Eigen::VectorXd x) {
	//				double regularize = lambda * (x - args[idx]).lpNorm<2>();
	//				std::vector<TrajectoryTemplate*> tmp = Well(x,geoPointsWells[idx]);
	//				return orderScore1(tmp, pCTrajectories, pMDTrajectories) + regularize; };
	//			PSOvalueType opt = PSO(score, minValues, maxValues, 30, 6, 150);
	//			trajectories.push_back(Well(opt.first, geoPointsWells[idx]));
	//			int cond = solve(trajectories.back());
	//			getOptData(opt);
	//			if (cond > 0) {
	//				std::cout << "Unbuild Trajectory\n";
	//				break;
	//			}
	//			writeDataOpt(order[i], order[i][id], opt, trajectories);
	//			pCTrajectories.push_back(allPointsCartesian(trajectories.back()));
	//			pMDTrajectories.push_back(allPointsMD(trajectories.back()));
	//		}
	//		pCTrajectories.clear();
	//		pMDTrajectories.clear();
	//		trajectories.clear();
	//		args.clear();
	//	}
	//}
}