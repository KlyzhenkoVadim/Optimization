#include "OptimizeWells.h"

double AzimuthChooser(Constraint cs, double r) 
{
	// r \in U[0,1]
	if (r < 0.5) {
		return cs.lMin.phi + 2 * r * (cs.lMax.phi - cs.lMin.phi);
	}
	else {
		return 180 + cs.lMin.phi + 2 * (r - 0.5) * (cs.lMax.phi - cs.lMin.phi);
	}
};

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
	IdealLength = (mainWell.back()->pointT1 - mainWell[0]->pointInitial).norm() + (mainWell.back()->pointT3 - mainWell.back()->pointT1).norm();
	if (IdealLength == 0)
		return 0.;
	mainDDI = DDI(mainWell, mainPCartesian);
	return mainLength / IdealLength + mainDDI;
}

std::vector<TrajectoryTemplate*> well2CHCH(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3)
{
	std::vector<TrajectoryTemplate*> well;
	well.reserve(3);
	double hold0 = 400 + (pT3[2] - 400) * x[0];//x[0]
	double DLSMax = 1.5;
	double R1 = x[1] < EPSILON ? 1 / EPSILON : 1800 / PI / x[1];
	double R2 = x[2] < EPSILON ? 1 / EPSILON : 1800 / PI / x[2];
	double R3 = x[3] < EPSILON ? 1 / EPSILON : 1800 / PI / x[3];
	double R4 = x[4] < EPSILON ? 1 / EPSILON : 1800 / PI / x[4];
	//double inc = x[5], azi = x[6];
	double inc = 35. + 20*x[5], azi = 360.*x[6];
	Eigen::Vector3d v = pT3 - pT1;
	v.normalize();
	double NS = pT1[0] + 2000 * (1 - 2 * x[7]);
	double EW = pT1[1] + 2000 * (1 - 2 * x[8]);
	double BETA = 100;
	Eigen::Vector3d tmpT3 = { NS,EW,pT1[2] - 150 }; 
	
	well.push_back(new Hold(pinit, 0, 0, hold0, typeHold::TVD)); // !!!
	well.back()->getTarget1Point();
	well.push_back(new CurveHoldCurveHold(well.back()->pointT1, 0, 0, R1, R2, tmpT3,inc,azi,BETA));
	well.back()->getTarget3Point();
	well.push_back(new CurveHoldCurveHold(well.back()->pointT3, inc, azi, R3, R4, pT1, pT3));
	return well;
}

void OptimizeHorizontal(const Eigen::Vector3d& pinit, const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3)
{
	std::function<double(const Eigen::VectorXd& x)> score = [&](const Eigen::VectorXd& x)
	{
		std::vector<TrajectoryTemplate*> tt = well2CHCH(x, pinit, pT1, pT3);
		double score = OneWellScore(tt);
		for (auto x : tt)
		{
			delete x;
		}
		return score;
	};
	std::vector<double> minValues{ 400,0,0,0,0,35,0,0,0}, maxValues{ pT1[2],1.5,1.5,1.5,1.5,55,360,1,1};
	PSOvalueType opt;
	for (size_t i = 0; i < 500; ++i) // 500
	{
		opt = PSO(score, minValues, maxValues, 5 * minValues.size(), minValues.size(),500);
		if (opt.second < 5)
		{
			break;
		}
	}
	getOptData(opt);
	//std::vector<TrajectoryTemplate*> tt = well2CHCH(opt.first,pinit,pT1,pT3);
	//std::vector<Eigen::Vector3d> pC = allPointsCartesian(tt);
	//writeDataCartesian(pC, "../output/Cartesian/horizontal_test.txt");
}

void OptimizeHorizontals(const Eigen::Vector3d& pinit, const std::vector<Eigen::Vector3d>& targets1, const std::vector<Eigen::Vector3d>& targets3)
{
	Eigen::Vector3d target1 = targets1[0], target3 = targets3[0];
	double TVDShift = 0;
	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector<std::vector<Eigen::Vector4d>> pMDWells;
	//std::vector<Constraint> cs = { { {900, 30, 110}, { 900,90,170 }}, { {2000,0,20},{2000,50,80} },{ {2600,0,20},{2600,40,80} }};
	
	Eigen::VectorXd x_reg;
	double l_reg = 0.1;

	std::function<double(const Eigen::VectorXd& x)> score = [&](const Eigen::VectorXd& x)
	{
		std::vector<TrajectoryTemplate*> tt = well2CHCH(x,pinit, target1, target3); //  wellSolverHorizontal(x, pinit, target1, target3);
		double score = orderScore1(tt, pCWells, pMDWells, TVDShift);
		for (auto x : tt)
		{
			delete x;
		}
		return score;// +l_reg * (x - x_reg).lpNorm<1>();
	};

	std::vector<PSOvalueType> opts;
	PSOvalueType tmpopt;
	std::vector<TrajectoryTemplate*> tt;
	//std::vector<double> minValues{ 400,0,0,0,0,35,0,0,0 }, maxValues{target1[2],1.5,1.5,1.5,1.5,55,360,1,1};
	std::vector<double> minValues{ 0,0,0,0,0,0,0,0,0 }, maxValues{ 1,1.5,1.5,1.5,1.5,1,1,1,1 };
	
	std::vector<Eigen::VectorXd> arr_xreg{
		Eigen::VectorXd{{(863.9466 - 400) / 2500, 0.7543149, 1.017725, 1.323914, 1.083109, (51.04885 - 35) / 20, 281.5318 / 360, 0.5132608, 0.4064161}},
		Eigen::VectorXd{{(410.3396 - 400) / (2900 - 400), 1.001136, 1.04123, 1.472675, 1.451098, (41.65297 - 35) / 20, 276.3627 / 360., 0.5348235, 0.4225081}}
	};
	/*
	std::vector<Eigen::VectorXd> arr_xreg{
	Eigen::VectorXd{{863.9466, 0.7543149, 1.017725, 1.323914, 1.083109, 51.04885, 281.5318, 0.5132608, 0.4064161}},
	Eigen::VectorXd{{410.3396, 1.001136, 1.04123, 1.472675, 1.451098, 41.65297, 276.3627, 0.5348235, 0.4225081}}};
	};*/

	/*for (size_t i = 0; i < target1.size(); ++i)
	{
		tmpopt = PSO(onescore, minValues, maxValues, 5 * maxValues.size(), maxValues.size(), 150);
		arr_xreg.push_back(tmpopt.first);
	}*/

	for (size_t i = 0; i < targets1.size(); ++i)
	{
		target1 = targets1[i];
		target3 = targets3[i];
		x_reg = arr_xreg[i];

		for (size_t j = 0; j < 100; ++j)
		{
			tmpopt = PSO(score, minValues, maxValues, 5 * maxValues.size(), maxValues.size(), 150); // WAS 150 ITERATIONS
			if (tmpopt.second < 5)//5
			{
				getOptData(tmpopt);
				TVDShift = std::max(TVDShift, tmpopt.first[0]);
				opts.push_back(tmpopt);
				tt = well2CHCH(tmpopt.first,pinit,target1,target3);
				int s = solve(tt);
				pCWells.push_back(allPointsCartesian(tt));
				pMDWells.push_back(allPointsMD(tt));
				for (auto y : tt)
					delete y;
				break;
			}
		}
		writeDataCartesian(pCWells[i], "../output/Cartesian/pHorizontal" + std::to_string(i + 1) + ".txt");
	}
	double sum_cost = 0.;
	for (auto x : opts)
	{
		sum_cost += x.second;
		std::cout << x.second << ", ";
	}
	std::cout << sum_cost << "\n";	
}

void testHorizontal(const Eigen::VectorXd& x,std::vector<TrajectoryTemplate*>& tt)
{
	int s = solve(tt);
	std::cout << "Condition of Well: " << s;
	std::vector<Eigen::Vector3d> pC = allPointsCartesian(tt);
	writeDataCartesian(pC, "output/Cartesian/pointsHorizontal.txt");
}

