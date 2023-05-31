#include "Slotting.h"

DirectionalDrillingSolver::DirectionalDrillingSolver(){}

void DirectionalDrillingSolver::setData(std::string filename)
{
	std::ifstream ifs(filename);
	bool flag = ifs.is_open();
	nlohmann::json jf = nlohmann::json::parse(ifs);
	Eigen::Vector3d tmp;
	std::pair<double, double> pp;
	pinit = { jf["Platform"]["coord"][0], jf["Platform"]["coord"][1] ,0 };
	for (auto x : jf["Platform"]["Targets"].items())
	{
		names.push_back(x.key());
		tmp[0] = x.value()["T1"][0];
		tmp[1] = x.value()["T1"][1];
		tmp[2] = x.value()["T1"][2];
		targets.push_back(tmp);
	}
	bool isnan;
	for (auto y : jf["Constraints"]["Layers"])
	{
		isnan = false;
		for (auto z : y) {
			if (z.is_null())
			{
				isnan = true;
				break;
			}
		}
		if (!isnan)
		{
			constraints.cs.push_back(Constraint{ {(y["TVD"]),(y["inc"][0]),(y["azi"][0])},{(y["TVD"]),(y["inc"][1]),(y["azi"][1])} });
		}
	}
	if (!jf["Constraints"]["MinHoldLength"].is_null())
		constraints.wc.minDepthFirstHold = jf["Constraints"]["MinHoldLength"];
	if (!jf["Constraints"]["MaxDLS"].is_null())
		constraints.wc.maxDLS = jf["Constraints"]["MaxDLS"];
	if (!jf["Constraints"]["MaxIncTarget"].is_null())
		constraints.MaxIncTarget = jf["Constraints"]["MaxIncTarget"];
	if (!jf["Constraints"]["MaxMD"].is_null())
		constraints.wc.maxMD = jf["Constraints"]["MaxMD"];
	if (!jf["Constraints"]["MaxAHD"].is_null())
	{
		constraints.wc.maxDistEastWest = jf["Constraints"]["MaxAHD"];
		constraints.wc.maxDistNorthSouth = jf["Constraints"]["MaxAHD"];
	}
	minValues.push_back(constraints.wc.minDepthFirstHold);
	minValues.push_back(0.);
	minValues.push_back(0.);
	minValues.push_back(0);
	minValues.push_back(0);
	minValues.push_back(0);
	maxValues.push_back(targets[0][2]);
	maxValues.push_back(constraints.wc.maxDLS);
	maxValues.push_back(constraints.wc.maxDLS);
	maxValues.push_back(constraints.MaxIncTarget);
	maxValues.push_back(360);
	maxValues.push_back(targets[0][2]);	
}

std::vector<TrajectoryTemplate*> DirectionalDrillingSolver::wellfunction(const std::vector<double>& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target)
{
	std::vector<TrajectoryTemplate*> tmp;
	tmp.push_back(new Hold(pinit, 0, 0, x[0]));
	tmp.back()->getTarget1Point();
	double R1 = abs(x[1]) < EPSILON ? 1 / EPSILON : 1800 / PI / x[1];
	double R2 = abs(x[2]) < EPSILON ? 1 / EPSILON : 1800 / PI / x[2];
	tmp.push_back(new CurveHoldCurveHold(tmp.back()->pointT1, 0, 0, R1, R2, target, x[3], x[4], x[5]));
	return tmp;
}

void DirectionalDrillingSolver::Optimize()
{
	Eigen::Vector3d target = targets[0];
	double TVDShift = 0;
	bool flag = true;

	std::function<double(const std::vector<double>&)> score = [&](const std::vector<double>& x)
	{
		std::vector<TrajectoryTemplate*> tmp = wellfunction(x, pinit, target);
		int s = solve(tmp);
		std::vector<Eigen::Vector3d> pCtmp = allPointsCartesian(tmp);
		std::vector<Eigen::Vector4d>  pMDtmp = allPointsMD(tmp);
		double cost = orderScore1(tmp, pCWells, pMDWells, TVDShift);
		double penLen = abs(constraints.wc.maxMD - 1 / EPSILON) < EPSILON ? 0 : PenaltyLength(allLength(tmp), constraints.wc.maxMD, 20);
		double penConstraints = constraints.cs.size() == 0 ? 0 : PenaltyConstraint(constraints.cs, pCtmp, pMDtmp, 100);
		double penIncWell = PenaltyIncWell(pMDtmp, 0, 65, 100);
		for (auto& x : tmp) {
			delete x;
		}
		return cost + penConstraints + penLen + penIncWell;
	};

	size_t index = 0;
	PSOvalueType opt;
	
	std::ofstream of, ofc;
	std::cout << "Optimization started...\n\n";
	for (size_t i = 0; i < targets.size(); ++i)
	{
		std::cout << "Iteration: " << i << "\n";
		target = targets[i];
		for (size_t j = 0; j < 50; ++j)
		{
			opt = PSO(score, minValues, maxValues, 3 * minValues.size(), minValues.size(),
				150, std::vector<double>(150, 0.7298), 1.49618, 1.49618);
			if (opt.cost - 10 < EPSILON)
			{
				std::cout << std::endl;
				break;
			}
		}
		std::vector<TrajectoryTemplate*> well = wellfunction(opt.optVec, pinit, target);
		int cond = solve(well);
		opts.push_back(opt);
		pCWells.push_back(allPointsCartesian(well));
		pMDWells.push_back(allPointsMD(well));
		TVDShift = std::max(TVDShift, opt.optVec[0]);
	}
}

void DirectionalDrillingSolver::writeResults(std::string outputdirpath) 
{
	for (size_t i = 0; i < opts.size(); ++i) {
		getOptData(opts[i]);
	}
	std::string tmp = outputdirpath + "cartesian",
		tmpInc = outputdirpath + "inclinometry";
	for (size_t i = 0; i < pCWells.size(); ++i) {
		writeDataCartesian(pCWells[i], tmp + names[i] + ".txt");
		writeInclinometry(pMDWells[i], tmpInc + names[i] + ".csv");
	}
}

