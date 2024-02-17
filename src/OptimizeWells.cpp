#include "OptimizeWells.h"

#include "CostFuncs.h"
#include "PSO.h"
#include "CurveHold.h"
#include "CurveHoldCurveHold.h"
#include "Hold.h"
#include "Curve.h"
#include "TestWells.h"
#include "Penalties.h"


double DlsToRadius(double dls) // dls - grad/10m
{
	return dls < EPSILON ? 1. / EPSILON : 1800. / PI / dls;
}

double OneWellScore(std::vector<TrajectoryTemplate*>& mainWell, double penalty) {
	double mainLength, mainDDI;
	int condition = solve(mainWell);
	if (condition != 0) {
		return penalty * (-condition) / mainWell.size(); // penalty * percent of incorrect templates.
	}
	mainLength = allLength(mainWell);
	std::vector<Eigen::Vector3d> mainPCartesian = allPointsCartesian(mainWell);
	double IdealLength;
	mainWell[0]->getInitPoint();
	mainWell.back()->getTarget1Point();
	mainWell.back()->getTarget3Point();
    IdealLength =
            (mainWell.back()->pointT1 - mainWell[0]->pointInitial).norm() +
            (mainWell.back()->pointT3 - mainWell.back()->pointT1).norm();
    if (IdealLength < std::numeric_limits<double>::epsilon())
		return 0.;
	mainDDI = DDI(mainWell, mainPCartesian);
	return mainLength / IdealLength + mainDDI;
}

std::vector<TrajectoryTemplate*> well2CHCH(const std::vector<double>& x,
                                           const Eigen::Vector3d& pinit,
                                           const Eigen::Vector3d& pT1,
                                           const Eigen::Vector3d& pT3) {
	std::vector<TrajectoryTemplate*> well;
	well.reserve(3);
	double hold0 = 100. + 900. * x[0];

	double R1 = DlsToRadius(x[1]);
	double R2 = DlsToRadius(x[2]);
	double R3 = DlsToRadius(x[3]);
	double R4 = DlsToRadius(x[4]);

	double inc = 180. * x[5];
	double azi = 360. * x[6];

	double NS = pT1[0] + 1000 * (1 - 2 * x[7]);
	double EW = pT1[1] + 1000 * (1 - 2 * x[8]);
	
	double BETA = 110.;
	Eigen::Vector3d tmpT3 = { NS,EW,pT1[2] - 150 };
	well.push_back(new Hold(pinit, 0, 0, hold0, typeHold::TVD));
	well.back()->getTarget1Point();
    well.push_back(new CurveHoldCurveHold(well.back()->pointT1, 0, 0, R1, R2,
										  tmpT3, inc, azi, BETA, 0.001,50));
	well.back()->getTarget3Point();
    well.push_back(new CurveHoldCurveHold(well.back()->pointT3, inc, azi, R3,
                                          R4, pT1, pT3, 0.001, 50));
	return well;
}

void OptimizeHorizontal(const Eigen::Vector3d& pinit,
                        const Eigen::Vector3d& pT1,
                        const Eigen::Vector3d& pT3) {
    std::function<double(const std::vector<double>& x)> score =
            [&](const std::vector<double>& x) {
		std::vector<TrajectoryTemplate*> tt = well2CHCH(x, pinit, pT1, pT3);
		double score = OneWellScore(tt);
		for (auto& x : tt)
		{
			delete x;
		}
		return score;
	};
	std::vector<double> minValues{ 400, 0, 0, 0, 0, 35, 0, 0, 0 };
	std::vector<double> maxValues{ pT1[2], 1.5, 1.5, 1.5, 1.5, 55, 360, 1, 1 };
	PsoValueType opt;
	for (size_t i = 0; i < 500; ++i)
	{
		opt = PSO(score, minValues, maxValues, 5 * minValues.size(),500);
		if (opt.cost < 5)
		{
			break;
		}
	}
	getOptData(opt);
}

void OptimizeHorizontals(const std::vector<Eigen::Vector3d>& pinits,
                         const std::vector<Eigen::Vector3d>& targets1,
                         const std::vector<Eigen::Vector3d>& targets3) {
	const std::string aloneOptDir{"alone_optimization/"};
	const std::string colriskOptDir{"colrisk_optimization/"};
	const std::string l1OptDir{"l1_optimization/"};
	const std::string l2OptDir{"l2_optimization/"};

	Eigen::Vector3d pinit = pinits[0];
	Eigen::Vector3d target1 = targets1[0], target3 = targets3[0];
	double TVDShift = 0;
	if (!(pinits.size() == targets1.size() && targets1.size() == targets3.size()))
	{
		std::cout << "error";
		return;
	}
	const size_t trajNumber = targets1.size();
	const size_t paramNumber = 9;
	const size_t psoNparticles = 40;
	const size_t psoNiterations = 200;
	std::vector<double> psoInertiaCoef(2 * psoNiterations, 0.7298);
	const double psoCoef = 1.49618;
	const double lambda_reg = 1.2;
	const size_t p_reg = 2;

	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector<std::vector<Eigen::Vector4d>> pMDWells;

	pCWells.reserve(trajNumber);
	pMDWells.reserve(trajNumber);

	std::vector<std::vector<double>> regVectors;
	regVectors.reserve(trajNumber);
	std::vector<double> x_reg;
	x_reg.resize(paramNumber);

	auto regNorm = [](const std::vector<double>& x,
                          const std::vector<double>& x_reg, int p)
	{
		size_t size = x.size();
		assert(size == x_reg.size());
		double res = 0.;
		if (p == 0)
			return 0.;
		for (size_t i = 0; i < size; ++i)
		{
			if (p == 1)
				res += std::abs(x[i] - x_reg[i]);
			else if (p == 2)
				res += (x[i] - x_reg[i]) * (x[i] - x_reg[i]);
		}
		if (p == 2)
			res = sqrt(res);
		return res;
	};

	std::function<double(const std::vector<double>& x)> score =
            [&](const std::vector<double>& x) {
		std::vector<TrajectoryTemplate*> trajectory = well2CHCH(x, pinit, target1,
														target3);
		double score = orderScore1(trajectory, pCWells, pMDWells, TVDShift,100.0);
		for (auto& x : trajectory)
		{
			delete x;
		}
		return score + lambda_reg * regNorm(x, x_reg,p_reg);
	};

	auto onescore = [&](const std::vector<double>& x)
	{
		std::vector<TrajectoryTemplate*> tt = well2CHCH(x, pinit, target1,
														target3);
		double score = OneWellScore(tt);
		for (auto& y : tt)
		{
			delete y;
		}
		return score;
	};

	std::vector<PsoValueType> optResults;
	PsoValueType currentOptRes;
	std::vector<TrajectoryTemplate*> currentTraj;
	
	std::vector<double> minValues(paramNumber, 0.0);
	std::vector<double> maxValues(paramNumber, 1.0);
	//FOR REGULARIZATION
	for (size_t i = 0; i < trajNumber; ++i)
	{
		pinit = pinits[i];
		target1 = targets1[i];
		target3 = targets3[i];
		
		for (size_t j = 0; j < 500; ++j) 
		{
			currentOptRes = PSO(onescore, minValues, maxValues, psoNparticles,
						 2 * psoNiterations, psoInertiaCoef, psoCoef, psoCoef);
			if (currentOptRes.cost < 4.0)
			{
				regVectors.push_back(currentOptRes.optVec);
				break;
			}
		}
		
		
	}//FOR REGULARIZATION
	
	//std::vector<std::array<size_t, 3>> indexes{{0, 1, 2}, { 0, 2, 1 },
	//  { 1,0,2 }, { 1,2,0 }, { 2,0,1 }, { 2,1,0 }};
	std::vector<std::array<size_t, 3>> optimizationOrders{{0, 1, 2}, {1, 2, 0}};
	std::vector<std::array<size_t, 3>> idx123(50, {0, 1, 2 });
	nlohmann::json jstat;
	for (const auto& currentOptOrder : idx123)
	{
		std::string name;
		nlohmann::json jRes;
		for (size_t idx = 0; idx < trajNumber;++idx)
		{
			size_t i = currentOptOrder[idx];
			name += std::to_string(i + 1);
			
			pinit = pinits[i];
			target1 = targets1[i];
			target3 = targets3[i];
			x_reg = regVectors[i];

			for (size_t j = 0; j < 100; ++j)
			{
				currentOptRes = PSO(score, minValues, maxValues, psoNparticles,
                             psoNiterations, psoInertiaCoef, psoCoef, psoCoef, true);
				double cost1 = currentOptRes.cost - lambda_reg * regNorm(currentOptRes.optVec, x_reg, p_reg);
				if (cost1 < 5.0) // TODO: уточнить, какое значение должно быть.. was (5.0)
				{
					getOptData(currentOptRes);
					TVDShift = std::max(TVDShift, currentOptRes.optVec[0]);
					optResults.push_back(currentOptRes);
					currentTraj = well2CHCH(currentOptRes.optVec, pinit, target1, target3);
					int s = solve(currentTraj);
					const auto& points_c = allPointsCartesian(currentTraj);
					const auto& points_inc = allPointsMD(currentTraj);
					double length = allLength(currentTraj);
					double ddi = DDI(currentTraj, points_c, false);
					pCWells.push_back(points_c);
					pMDWells.push_back(points_inc);
					for (auto& y : currentTraj)
						delete y;
					std::cout << "trajectory #" + std::to_string(i + 1) +
								" optimized in " + std::to_string(j + 1) +
								" iterations." << std::endl;
					nlohmann::json j;
					j["order"] = currentOptOrder;
					j["name"] = "well-" + std::to_string(i + 1);
					j["pInit"] = pinit;
					j["pT1"] = target1;
					j["pT3"] = target3;

					auto& jOpt = j["opt-result"];
					jOpt["reg"] = p_reg;
					jOpt["l-reg"] = lambda_reg;
					jOpt["cost"] = cost1; // tmpopt.cost - l_reg * regNorm(tmpopt.optVec, x_reg, p_reg);
					jOpt["DDI"] = ddi;
					jOpt["length"] = length;
					jOpt["opt-vec"] = currentOptRes.optVec;
					jOpt["VHold"] = currentOptRes.optVec[0] * 900. + 100.;
					jOpt["dls1"] = currentOptRes.optVec[1];
					jOpt["dls2"] = currentOptRes.optVec[2];
					jOpt["dls3"] = currentOptRes.optVec[3];
					jOpt["dls4"] = currentOptRes.optVec[4];
					/*auto& jSf = jOpt["separation"];
					for (int k = 0; k < idx; ++k)
					{
						nlohmann::json jt;
						jt["name-ref"] = "well-" + std::to_string(i + 1);
						jt["name-offset"] = jRes[k]["name"].get<std::string>();
						jt["sf"] = sepFactor(points_c, points_inc, pCWells[k], pMDWells[k], 0, false);
						jSf.push_back(jt);

					}*/
					jOpt["inc"] = 180. * currentOptRes.optVec[5];
					jOpt["azi"] = 360. * currentOptRes.optVec[6];
					jOpt["NS"] = target1[0] + 1000. * (1. - 2. * currentOptRes.optVec[7]);
					jOpt["EW"] = target1[1] + 1000. * (1. - 2. * currentOptRes.optVec[8]);

					jRes.push_back(j);

					break;
				}
			}
		}

		pCWells.clear();
		pMDWells.clear();
		jstat.push_back(jRes);
	}
}

void OptimizeTogether(const std::vector<Eigen::Vector3d>& pinits,
                      const std::vector<Eigen::Vector3d>& targets1,
                      const std::vector<Eigen::Vector3d>& targets3) {
	Eigen::Vector3d pinit = pinits[0];
	Eigen::Vector3d target1 = targets1[0];
	Eigen::Vector3d target3 = targets3[0];
	double TVDShift = 0;
	if (!(pinits.size() == targets1.size() &&
		targets1.size() == targets3.size()))
	{
		std::cout << "error";
		return;
	}
	size_t size = targets1.size();
	size_t nParams = 9;
	size_t nParticles = 20;
	size_t nIterations = 300;

	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector<std::vector<Eigen::Vector4d>> pMDWells;

	pCWells.reserve(size);
	pMDWells.reserve(size);

	std::vector<std::vector<double>> arr_xreg;
	arr_xreg.reserve(size);
	std::vector<double> x_reg;
	x_reg.resize(nParams);

	double l_reg = 0.01;
    auto regNorm = [](const std::vector<double>& x,
                        const std::vector<double>& x_reg, int p) {
		size_t size = x.size();
		assert(size == x_reg.size());
		double res = 0.;
		for (size_t i = 0; i < size; ++i)
		{
			if (p == 1)
				res += std::abs(x[i] - x_reg[i]);
			else if (p == 2)
				res += (x[i] - x_reg[i]) * (x[i] - x_reg[i]);
		}
		if (p == 2)
			res = sqrt(res);
		return res;
	};

	auto onescore = [&](const std::vector<double>& x)
	{
		std::vector<TrajectoryTemplate*> tt =
          well2CHCH(x, pinit, target1, target3);
		int condition = solve(tt);
		if (!condition)
			return 1000.;
		double score = OneWellScore(tt,1000);
		for (auto& y : tt)
		{
			delete y;
		}
		return score;
	};

	auto scoretogether = [&](const std::vector<double>& x)
	{
		double result = 0.;
		for (size_t i = 0; i < size; ++i)
		{
			std::vector<double> xi;
			xi.resize(nParams);
			pinit = pinits[i];
			target1 = targets1[i];
			target3 = targets3[i];
			for (size_t j = 0; j < nParams; ++j)
				xi[j] = x[nParams * i + j];
			double score = onescore(xi);
			if (score > 999.)
				return 1000.;
			result += score;
		}
		return result;
	};

	PsoValueType opt;
	std::vector<double> minValues(nParams * size, 0.0);
	std::vector<double> maxValues(nParams * size, 1.0);
	for (size_t iter = 0; iter < 500; ++iter)
	{
		opt = PSO(scoretogether, minValues, maxValues, nParticles,
			nIterations, std::vector<double>(nIterations, 0.9),
			1.49, 1.49, true); // WAS 150 ITERATIONS
		if (opt.cost < 100.)
		{
			break;
		}
	}
	getOptData(opt);
	nlohmann::json jRes;
	nlohmann::json jj;
	jj["PSO-history"] = opt.costHist;
	jj["cost"] = opt.cost;
	jj["opt-vec"] = opt.optVec;
	jRes.push_back(jj);
	for (size_t i = 0; i < size; ++i)
	{
		std::vector<double> xi;
		xi.resize(nParams);
		for (size_t j = 0; j < nParams; ++j)
		{
			xi[j] = opt.optVec[nParams * i + j];
		}
		auto tt = well2CHCH(xi, pinits[i], targets1[i], targets3[i]);
		int condition = solve(tt);
		const auto& points_c = allPointsCartesian(tt);
		const auto& points_inc = allPointsMD(tt);
		double length = allLength(tt);
		double ddi = DDI(tt, points_c, false);
		for (auto& y : tt)
			delete y;
		nlohmann::json j;
		if(!condition)
		{
			j["name"] = "well-" + std::to_string(i + 1);
			j["pInit"] = pinit;
			j["pT1"] = target1;
			j["pT3"] = target3;
			j["reg"] = 0;
			j["DDI"] = ddi;
			j["length"] = length;
			j["opt-vec"] = xi;
			j["points"] = points_c;
		}
		jRes.push_back(j);
	}
	writeGG("../output/together_worisks0.json", jRes);
}

void testHorizontal(const Eigen::VectorXd& x,std::vector<TrajectoryTemplate*>& tt)
{
	int s = solve(tt);
	std::cout << "Condition of Well: " << s;
	std::vector<Eigen::Vector3d> pC = allPointsCartesian(tt);
	writeDataCartesian(pC, "output/Cartesian/pointsHorizontal.txt");
}

void writeGG(const std::string& filename, const nlohmann::json& jresults)
{
	std::ofstream ofs(filename);
	ofs << jresults << std::endl;
	ofs.close();
}

void aloneOpt()
{
	Eigen::Vector3d pinit{0, 0, 0};
	Eigen::Vector3d pt1{5., 0, 1000.};
	Eigen::Vector3d pt3{1000,0,1000};
	double mDls = 1.;
	double mHold = 400.;
	std::vector<double> minValues{100.0,0.0};
	std::vector<double> maxValues{300.0,1.0};
	size_t nAgents = 10;
	size_t nDim = 2;
	size_t nIterations = 150;
	double w = 0.7298;
	double c1 = 1.49618;
	auto horizontalProfile = [&](const std::vector<double>& x)
	{//x[0] - hold, x[1] - dls1 , x[2] - dls2
		double hold = x[0];
		double r1 = DlsToRadius(x[1]);
		double r2 = DlsToRadius(x[2]);
		std::vector<TrajectoryTemplate*> trajectory;
		trajectory.reserve(3);
		trajectory.push_back(new Hold(pinit, 0, 0, hold));
		trajectory.back()->getTarget3Point();
		trajectory.push_back(new CurveHoldCurveHold(trajectory.back()->pointT3,
							0, 0, r1, r2, pt1, pt3));
		return trajectory;
	};

	auto well = [&](const std::vector<double>& x)
	{
		double hold = x[0];
		double r = DlsToRadius(x[1]);
		std::vector<TrajectoryTemplate*> tt;
		tt.reserve(2);
		tt.push_back(new Hold(pinit, 0., 0., hold));
		tt.back()->getTarget3Point();
		tt.push_back(new CurveHold(tt.back()->pointT3, pt1, 0., 0., r));
		return tt;
	};

	auto score = [&](const std::vector<double>& x)
	{
		auto tt = well(x);
		double result = OneWellScore(tt,35);
		for (auto& t : tt)
		{
			delete t;
		}
		return result;
	};

	nlohmann::json jres;
	auto& jbeta = jres["beta"];
	auto& jdls = jres["dls"];
	auto& jcost = jres["cost"];

	for (size_t i = 0; i < 100; ++i)
	{
		double beta = 100. + 200.*double(i) / 99.;
		for (size_t j = 0; j < 100; ++j)
		{
			double dls = double(j) / 99.;
			double cost = score(std::vector<double>{beta, dls});
	
			jbeta.push_back(beta);
			jdls.push_back(dls);
			jcost.push_back(cost);
		}
	}
	
	nlohmann::json jt;
	jt["name"] = "example";
    std::vector<TrajectoryTemplate*> tExample =
        well(std::vector<double>{mHold, mDls});
	size_t condition = solve(tExample);
	if (condition)
		std::cout << "error-tExample";
	auto points = allPointsCartesian(tExample);
	double len_example = allLength(tExample);
	
	jt["points"] = points;
	writeGG("../output/cost_graph/2.json", jt);
}
