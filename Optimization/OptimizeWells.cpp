#include "OptimizeWells.h"

double AzimuthChooser(Constraint cs, double r) {
	// r \in U[0,1]
	if (r < 0.5) {
		return cs.lMin.phi + 2 * r * (cs.lMax.phi - cs.lMin.phi);
	}
	else {
		return 180 + cs.lMin.phi + 2 * (r - 0.5) * (cs.lMax.phi - cs.lMin.phi);
	}
};

std::vector<TrajectoryTemplate*> wellCHCH(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target) {
	std::vector<TrajectoryTemplate*> tmp;
	tmp.push_back(new Hold(pinit, 0, 0, x[0]));
	tmp.back()->getTarget1Point();
	double R1 = abs(x[1]) < EPSILON ? 1 / EPSILON : 1800 / PI / x[1];
	double R2 = abs(x[2]) < EPSILON ? 1 / EPSILON : 1800 / PI / x[2];
	tmp.push_back(new CurveHoldCurveHold(tmp.back()->pointT1, 0,0, R1, R2, target, x[3], x[4],x[5]));
	return tmp;
}

double OneWellScore(std::vector<TrajectoryTemplate*>& mainWell, double penalty) {
	double mainLength, mainDDI;
	int condition = solve(mainWell);
	if (condition != 0) {
		return penalty * condition / mainWell.size(); // penalty * percent of incorrect templates.
	}
	mainLength = allLength(mainWell);
	std::vector<Eigen::Vector3d> mainPCartesian = allPointsCartesian(mainWell);
	std::vector<Eigen::Vector4d> mainPMD = allPointsMD(mainWell);
	double IdealLength;
	mainWell[0]->getInitPoint();
	mainWell.back()->getTarget1Point();
	mainWell.back()->getTarget3Point();
	IdealLength = (mainWell.back()->pointT1 - mainWell[0]->pointInitial).norm() + (mainWell.back()->pointT3 - mainWell.back()->pointT1).norm();
	if (IdealLength == 0)
		return 0.;
	mainDDI = DDI(mainPCartesian, mainPMD);

	return mainLength / IdealLength + mainDDI;

}

std::vector<TrajectoryTemplate*> Well(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs) 
{
	// Hold1 inc1 azi1 inc2 azi2 inc3 azi3 tvd1 tvd2 tvd3 R
	std::vector<TrajectoryTemplate*> tr;
	tr.push_back(new Hold(pinit, 0, 0, x[0])); // holdLength
	tr.back()->getTarget1Point();
	std::vector<Layer> layers{ {x[0],0,0}, { 900,x[1],AzimuthChooser(cs[0],x[2])},
		{2000,x[3],AzimuthChooser(cs[1],x[4])},{2600,x[5],AzimuthChooser(cs[2],x[6])} };
	double tvd0 = x[0] + (900 - x[0]) * x[7];
	std::vector<double> tvds{ tvd0,x[8],x[9] };// tvd (x[0],900)// tvd (900,2000)// tvd (2000,2600)
	double RCurveHold = x[10] < EPSILON ? 1 / EPSILON : 1800 / PI / x[10]; // x[10] = dls
	for (size_t idx = 1; idx < layers.size(); ++idx) {
		tr.push_back(new Curve(tr.back()->pointT1, layers[idx - 1].theta, layers[idx - 1].phi, layers[idx].theta, layers[idx].phi, tvds[idx - 1], TypeCurve::TVD));
		tr.back()->getTarget1Point();
		tr.push_back(new Hold(tr.back()->pointT1, layers[idx].theta, layers[idx].phi, layers[idx].TVD, typeHold::TVD));
		tr.back()->getTarget1Point();
	}
	tr.push_back(new CurveHold(tr.back()->pointT1, target, layers.back().theta, layers.back().phi, RCurveHold));	
	return tr;
};

void OptimizeWells(const Eigen::Vector3d& pinit, const std::vector<Eigen::Vector3d>& targets, const std::vector<Constraint>& cs,
					const std::vector<double>& minValues, const std::vector<double>& maxValues, std::vector<std::vector<Eigen::Vector3d>>& pCWells,
					std::vector<std::vector<Eigen::Vector4d>>& pMDWells, std::vector<PSOvalueType>& opts)
{
	Eigen::Vector3d target = targets[0];
	double TVDShift = 0;
	bool flag = true;
	std::function<double(const Eigen::VectorXd&)> scoreOne = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> tmp = Well(x, pinit, target, cs);
		int condition = solve(tmp);
		if (condition != 0) {
			return 1000. * condition / tmp.size(); // penalty * percent of incorrect templates.
		}
		double  length = allLength(tmp), dlspenalty = PenaltyDLS(tmp,500);
		double IdealLength = length;
		tmp[0]->getInitPoint();
		tmp.back()->getTarget1Point();
		tmp.back()->getTarget3Point();
		IdealLength = (tmp.back()->pointT1 - tmp[0]->pointInitial).norm() + (tmp.back()->pointT3 - tmp.back()->pointT1).norm();
		for (auto x : tmp) {
			delete x;
		}
		return length / IdealLength + PenaltyLength(length,5000,100) + dlspenalty;
	};

	std::function<double(const Eigen::VectorXd&)> CollRiskScore = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> currWell = Well(x, pinit, target, cs);
		double score = orderScore1(currWell, pCWells, pMDWells, TVDShift);
		double length = allLength(currWell), dlsPen = PenaltyDLS(currWell, 10), incPen = PenaltyIncTarget(currWell, 10);
		double alphaCH = PenaltyAlphaTarget(currWell, 10);
		for (auto x : currWell)
			delete x;
		return score + PenaltyLength(length,5000,10) + dlsPen + incPen;
	};

	size_t index = 0;

	while (flag) {
		flag = false;
		for (size_t idd = index; idd < targets.size(); ++idd) {
			target = targets[idd];
			PSOvalueType opt = PSO(CollRiskScore, minValues, maxValues, 15, 11, 150);
			if(opt.second > 10){
				flag = true;
				index = idd;
				break;
			}
			std::cout << index <<"\n";
			index = idd + 1;
			std::vector<TrajectoryTemplate*> well = Well(opt.first, pinit, target, cs);
			int cond = solve(well);
			opts.push_back(opt);
			pCWells.push_back(allPointsCartesian(well));
			pMDWells.push_back(allPointsMD(well));
			TVDShift = std::max(TVDShift, opt.first[0]);
		}
	}
}

void writeIteration(size_t id)
{
	std::string filename = "IterationInput.txt";
	std::ofstream output;
	output.open(filename, std::ios::out || std::ios::app);
	output << id << std::endl;
	output.close();
}

std::vector<TrajectoryTemplate*> wellSolverHorizontal(const Eigen::VectorXd& x, const Eigen::Vector3d& pInit, const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3)
{	
	double MAXDLS = 1.5;
	double R1 = x[1] < EPSILON ? 1 / EPSILON : 1800 / x[1] / PI;
	double R2 = x[2] < EPSILON ? 1 / EPSILON : 1800 / x[2] / PI;
	double TVD1 = x[0] + (pT3[2] - x[0]) * x[3];
	double TVD2 = TVD1 + (pT3[2] - TVD1) * x[4];
	double inc1 = 180/PI*asin((TVD1-x[0])*PI*MAXDLS/1800*x[5]), azi1 = x[6];
	std::vector<TrajectoryTemplate*> well;
	well.push_back(new Hold(pInit, 0, 0, x[0], typeHold::TVD));
	well.back()->getTarget1Point();
	well.push_back(new Curve(well.back()->pointT1, 0, 0, inc1, azi1, TVD1, TypeCurve::TVD));
	well.back()->getTarget1Point();
	well.push_back(new Hold(well.back()->pointT1, inc1, azi1, TVD2, typeHold::TVD));
	well.back()->getTarget1Point();
	well.push_back(new CurveHoldCurveHold(well.back()->pointT1,inc1,azi1, R1, R2, pT1, pT3));
	return well;
}

void OptimizeHorizontal(const Eigen::Vector3d& pinit, const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3)
{
	std::function<double(const Eigen::VectorXd& x)> score = [&](const Eigen::VectorXd& x)
	{
		std::vector<TrajectoryTemplate*> tt = wellSolverHorizontal(x, pinit, pT1, pT3);
		double score = OneWellScore(tt);
		for (auto x : tt)
		{
			delete x;
		}
		return score;
	};
	std::vector<double> minValues{ 40,0,0,0,0,0,0 }, maxValues{ pT3[2],1.5,1.5,1,1,50,360 };
	PSOvalueType opt;
	for (size_t i = 0; i < 500; ++i)
	{
		opt = PSO(score, minValues, maxValues, 5 * minValues.size(), minValues.size());
		if (opt.second < 10)
		{
			break;
		}
	}
	getOptData(opt);
	testHorizontal(opt.first, pinit, pT1, pT3);
}

void OptimizeHorizontals(const Eigen::Vector3d& pinit, const std::vector<Eigen::Vector3d>& targets1, const std::vector<Eigen::Vector3d>& targets3)
{
	Eigen::Vector3d target1 = targets1[0], target3 = targets3[0];
	double TVDShift = 0;
	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector<std::vector<Eigen::Vector4d>> pMDWells;
	std::vector<Constraint> cs = { { {900, 30, 110}, { 900,90,170 }}, { {2000,0,20},{2000,50,80} } };// , { {2600,0,20},{2600,40,80} }};

	std::function<double(const Eigen::VectorXd& x)> score = [&](const Eigen::VectorXd& x)
	{
		std::vector<TrajectoryTemplate*> tt = wellSolverHorizontal(x, pinit, target1, target3);
		double score = orderScore1(tt, pCWells, pMDWells, TVDShift);
		int s = solve(tt);
		auto pMD = allPointsMD(tt);
		auto pC = allPointsCartesian(tt);
		//double penConstr = PenaltyConstraint(cs, pC, pMD, 100);
		for (auto x : tt)
		{
			delete x;
		}
		return score;// penConstr;
	};
	std::vector<PSOvalueType> opts;
	PSOvalueType tmpopt;
	std::vector<TrajectoryTemplate*> tt;
	std::vector<double> minValues{ 400,0,0,0,0,0,0}, maxValues{2900,1.5,1.5,1,1,1,360};;
	for (size_t i = 0; i < targets1.size(); ++i)
	{
		target1 = targets1[i];
		target3 = targets3[i];
		for (size_t j = 0; j < 100; ++j)
		{
			tmpopt = PSO(score, minValues, maxValues, 5 * maxValues.size(), maxValues.size(), 150);
			if (tmpopt.second < 10)
			{
				getOptData(tmpopt);
				TVDShift = std::max(TVDShift, tmpopt.first[0]);
				opts.push_back(tmpopt);
				tt = wellSolverHorizontal(tmpopt.first, pinit, target1, target3);
				int s = solve(tt);
				pCWells.push_back(allPointsCartesian(tt));
				pMDWells.push_back(allPointsMD(tt));
				for (auto y : tt)
					delete y;
				break;
			}
		}
		writeDataCartesian(pCWells[i], "output/Cartesian/pHorizontal" + std::to_string(i + 1) + ".txt");
	}
	
}

void testHorizontal(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3)
{
	std::vector<TrajectoryTemplate*> tt = wellSolverHorizontal(x,pinit,pT1,pT3);
	int s = solve(tt);
	std::vector<Eigen::Vector3d> pC = allPointsCartesian(tt);
	writeDataCartesian(pC, "output/Cartesian/pointsHorizontal.txt");

}

void getXfromFile(Eigen::VectorXd& x, std::string word)
{
	size_t j = -1, lineNum = 0;
	std::string tmp;
		for (size_t i = 0; i < word.size(); ++i)
		{
			if (word[i] != ',')
			{
				tmp.push_back(word[i]);
			}
			else 
			{
				j += 1;
				x[j] = std::stod(tmp, nullptr);
				tmp.clear();
			}
		}
}

void MakePointsEvo(const Eigen::Vector3d& pinit,const Eigen::Vector3d& target,const std::string& filename,const std::string& numtarget)
{
	std::ifstream file;
	file.open(filename);
	Eigen::VectorXd x{{0,0,0,0,0,0}};
	std::string dir = "C:/Users/klyzhenko.vs/Desktop/optimization/Optimization/Optimization/output/to_GIF/"+numtarget+"/";
	std::vector<TrajectoryTemplate*> tmp;
	std::vector<Eigen::Vector3d> pC;
	std::string word;
	size_t lineNum = 0;
	while (std::getline(file,word))
	{
		getXfromFile(x, word);
		tmp = wellCHCH(x,pinit,target);
		int s = solve(tmp);
		if (s == 0)
		{
			pC = allPointsCartesian(tmp);
			std::string fileNum = std::to_string(lineNum);
			if (lineNum < 10)
				fileNum = "00" + fileNum;
			if (lineNum >= 10 and lineNum < 100)
				fileNum = "0" + fileNum;

			writeDataCartesian(pC, dir+fileNum+ ".txt");
		}
		for (auto w : tmp)
		{
			delete w;
		}
		pC.clear();
		lineNum += 1;
	}
	file.close();
}