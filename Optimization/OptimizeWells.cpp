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

void ShowOptData(PSOvalueType opt, std::vector<Constraint> cs) {
	std::cout << "Optimization result:\n_____________________________________\n";
	std::cout << "CostValue: " << opt.second << '\n';
	std::cout << "Hold0TVD: " << opt.first[0] << '\n';
	for (size_t i = 0; i < 3; ++i) {
		std::cout << "Curve" << i + 1 << " (inc" << i + 1 << ", azi" << i + 1 << "): " << opt.first[2 * i + 1] << ", " << AzimuthChooser(cs[i], opt.first[2 * (i + 1)]) << "\n";
		if (i == 0) {
			std::cout << "Hold" << i + 1 << "TVDstart: " << opt.first[0] + opt.first[i + 7] * (900 - opt.first[0]) << "\n";
		}
		else {
			std::cout << "Hold" << i + 1 << "TVDstart: " << opt.first[i + 7] << "\n";
		}
	}

	std::cout << "CurveHoldDLS: " << 1800 / PI / opt.first[10] << "\n";
}

double PenaltyLength(double length, double penalty) {
	double answer = length > 5000 ? penalty : 0;
	return answer;
};

double PenaltyDLS(std::vector<TrajectoryTemplate*>& Well, double penalty) {
	double alpha, dotprod, dls;
	int s = solve(Well);
	for (size_t i = 0; i < Well.size() - 1; ++i) {
		Well[i]->getInitPoint(CoordinateSystem::MD);
		Well[i]->getTarget1Point(CoordinateSystem::MD);
		Eigen::Vector3d tprev{ Well[i]->pointInitialMD[1],Well[i]->pointInitialMD[2],Well[i]->pointInitialMD[3] };
		Eigen::Vector3d tnext{ Well[i]->pointMDT1[1],Well[i]->pointMDT1[2],Well[i]->pointMDT1[3] };
		double length = Well[i]->pointMDT1[0] - Well[i]->pointInitialMD[0];
		dotprod = tprev.dot(tnext) / tprev.norm() / tnext.norm();
		if (length < EPSILON)
			continue;
		dotprod = abs(dotprod - 1) < EPSILON ? 1 : dotprod;
		alpha = 180 / PI * acos(dotprod);
		dls = 10 * alpha / length;
		if (dls > 1.5) {
			return penalty;
		}
	}
	return 0;
};

double PenaltyIncTarget(std::vector<TrajectoryTemplate*>& Well, double penalty) {
	int s = solve(Well);
	Well.back()->getTarget1Point(CoordinateSystem::MD);
	Eigen::Vector3d tangent = {Well.back()->pointMDT1[1],Well.back()->pointMDT1[2],Well.back()->pointMDT1[3] };
	std::pair<double, double> incAzi = CartesianToSpherical(tangent);
	if (incAzi.first - 40 > EPSILON) {
		return penalty;
	}
	return 0;
};

double PenaltyAlphaTarget(std::vector<TrajectoryTemplate*>& Well, double penalty) {
	int s = solve(Well);
	Well.back()->getInitPoint(CoordinateSystem::MD);
	Well.back()->getTarget1Point(CoordinateSystem::MD);
	Eigen::Vector3d t1 = { Well.back()->pointInitialMD[1],Well.back()->pointInitialMD[2],Well.back()->pointInitialMD[3] };
	Eigen::Vector3d t2 = { Well.back()->pointMDT1[1] ,Well.back()->pointMDT1[2] ,Well.back()->pointMDT1[3] };
	double dotprod = t1.dot(t2) / t1.norm() / t2.norm();
	dotprod = fabs(dotprod - 1) < EPSILON ? 1 : dotprod;
	double alpha = acos(dotprod);
	if (alpha - PI / 2 > EPSILON) {
		return penalty;
	}
	return 0;

}

int findTVD(const std::vector<Eigen::Vector3d>& pC, double TVD) {
	int n = pC.size();
	int l = 0, r = n, mid = (l + r) / 2;
	while (!(l >= r)) {
		mid = (l + r) / 2;
		if (abs(pC[mid][2] - TVD) < EPSILON)
			return mid;
		if (pC[mid][2] - TVD > EPSILON) {
			r = mid;
		}
		else {
			l = mid + 1;
		}
	}
	return -(l + 1);
}

double PenaltyConstraint(std::vector<Constraint> cs, std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD, double penalty) {
	int id = 0;
	double pen = 0;
	std::pair<double, double>  incAzi;
	for (size_t i = 0; i < cs.size(); ++i) {
		id = findTVD(pC, cs[i].lMin.TVD);
		if (id < 0) {
			id = -id - 1;
		}
		incAzi = CartesianToSpherical(Eigen::Vector3d{ pMD[id][1],pMD[id][2],pMD[id][3] });
		if (!(incAzi.first >= cs[i].lMin.theta and incAzi.first <= cs[i].lMax.theta)) {
			pen += penalty / cs.size() / 2;
		}
		if ((incAzi.first > EPSILON) and !((incAzi.second > cs[i].lMin.phi and incAzi.second < cs[i].lMax.phi) or (incAzi.second > cs[i].lMin.phi + 180 and incAzi.second < cs[i].lMax.phi + 180))) {
			pen += penalty / cs.size() / 2;
		}
	}
	return pen;
}

std::vector<TrajectoryTemplate*> wellCHCH(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target) {
	std::vector<TrajectoryTemplate*> tmp;
	tmp.push_back(new Hold(pinit, 0, 0, x[0]));
	tmp.back()->getTarget1Point();
	double R1 = abs(x[1]) < EPSILON ? 1 / EPSILON: 1800 / PI / x[1];
	double R2 = abs(x[2]) < EPSILON ? 1/EPSILON : 1800 / PI / x[2];
	
	tmp.push_back(new CurveHoldCurveHold(tmp.back()->pointT1, 0, 0, R1, R2, target, x[3], x[4]));
	return tmp;
}

std::vector<TrajectoryTemplate*> Well(const Eigen::VectorXd& x, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs) {
	// const Eigen::Vector3d & pInit, const Eigen::Vector3d & pTarget
	// Hold1 inc1 azi1 inc2 azi2 inc3 azi3 tvd1 tvd2 tvd3 R
	std::vector<TrajectoryTemplate*> tr;
	tr.push_back(new Hold(pinit, 0, 0, x[0])); // holdLength
	tr.back()->getTarget1Point();
	std::vector<Layer> layers{ {x[0],0,0}, { 900,x[1],AzimuthChooser(cs[0],x[2])},
		{2000,x[3],AzimuthChooser(cs[1],x[4])},{2600,x[5],AzimuthChooser(cs[2],x[6])} };
	double tvd0 = x[0] + (900 - x[0]) * x[7];
	std::vector<double> tvds{ tvd0,x[8],x[9] };
	// tvd (x[0],900)
	// tvd (900,2000)
	// tvd (2000,2600)
	double RCurveHold = x[10] < EPSILON ? 1 / EPSILON : 1800 / PI / x[10]; // x[10] = dls
	//  theta (30,90) phi (110-170)
	//	theta (0,50) phi (200-260)	 выбрал  было (20-80)
	//	theta (0,40) phi (200-260)	выбрал	 было (20-80)
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
		return length / IdealLength + PenaltyLength(length,100) + dlspenalty;
	};

	std::function<double(const Eigen::VectorXd&)> CollRiskScore = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> currWell = Well(x, pinit, target, cs);
		double score = orderScore1(currWell, pCWells, pMDWells, TVDShift);
		double length = allLength(currWell), dlsPen = PenaltyDLS(currWell, 500), incPen = PenaltyIncTarget(currWell, 100);
		double alphaCH = PenaltyAlphaTarget(currWell, 10);
		for (auto x : currWell)
			delete x;
		return score + PenaltyLength(length, 100) + dlsPen + incPen;
	};

	while (flag) {
		flag = false;
		for (size_t idd = 0; idd < targets.size(); ++idd) {
			target = targets[idd];
			PSOvalueType opt = PSO(CollRiskScore, minValues, maxValues, 30, 11, 150);

			if(opt.second > 3.8){
				flag = true;
				pCWells.clear();
				pMDWells.clear();
				opts.clear();
				break;
			}
			std::vector<TrajectoryTemplate*> well = Well(opt.first, pinit, target, cs);
			int cond = solve(well);
			opts.push_back(opt);
			pCWells.push_back(allPointsCartesian(well));
			pMDWells.push_back(allPointsMD(well));
			TVDShift = std::max(TVDShift, opt.first[0]);
		}
	}
}

void OptimizeCHCHWells(const Eigen::Vector3d& pinit, const std::vector<Eigen::Vector3d>& targets, const std::vector<Constraint>& cs,
	const std::vector<double>& minValues, const std::vector<double>& maxValues, std::vector<std::vector<Eigen::Vector3d>>& pCWells,
	std::vector<std::vector<Eigen::Vector4d>>& pMDWells, std::vector<PSOvalueType>& opts) 
{
	Eigen::Vector3d target = targets[0];
	double TVDShift = 0;
	bool flag = true;

	std::function<double(const Eigen::VectorXd&)> scoreCHCH = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> tmp = wellCHCH(x, pinit, target);
		int s = solve(tmp);
		std::vector<Eigen::Vector3d> pCtmp = allPointsCartesian(tmp);
		std::vector<Eigen::Vector4d>  pMDtmp = allPointsMD(tmp);
		double cost = orderScore1(tmp, pCWells, pMDWells), dlsPen = PenaltyDLS(tmp, 500), length = allLength(tmp);
		double penConstraints = PenaltyConstraint(cs, pCtmp, pMDtmp, 100);
		for (auto x : tmp) {
			delete x;
		}
		return cost + penConstraints + PenaltyLength(length,100);
	};

	while (flag) {
		flag = false;
		for (size_t idd = 0; idd < targets.size(); ++idd) {
			target = targets[idd];
			PSOvalueType opt = PSO(scoreCHCH, minValues, maxValues, 2, 5, 150);
			if (opt.second > 20) {
				flag = true;
				pCWells.clear();
				pMDWells.clear();
				opts.clear();
				break;
			}
			std::vector<TrajectoryTemplate*> well = wellCHCH(opt.first, pinit, target);
			int cond = solve(well);
			opts.push_back(opt);
			pCWells.push_back(allPointsCartesian(well));
			pMDWells.push_back(allPointsMD(well));
			TVDShift = std::max(TVDShift, opt.first[0]);
		}
	}
}

void FindBestWellHead(const std::vector<Eigen::Vector3d>& targets) 
{
	std::vector<double> minVs{ 5800236,680000 }, maxVs{ 5810000,685000 };
	std::function<double(const Eigen::VectorXd&)> headScore = [&](const Eigen::VectorXd& x) {
		Eigen::Vector3d InitialPoint = { x[0],x[1],0. };
		double azi, inc, cost = 0;
		for (size_t i = 0; i < 4; ++i) {
			inc = CartesianToSpherical(targets[i] - InitialPoint).first;
			azi = CartesianToSpherical(targets[i] - InitialPoint).second;
			if (azi - 110 > EPSILON and azi - 170 < EPSILON and inc - 0 > EPSILON and inc - 40 < EPSILON) {
				cost += (targets[i] - InitialPoint).norm();
				continue;
			}
			if (azi - 110 - 180 > EPSILON and azi - 170 - 180 < EPSILON and inc - 0 > EPSILON and inc - 40 < EPSILON) {
				cost += (targets[i] - InitialPoint).norm();
				continue;
			}
			cost += 1 / EPSILON / EPSILON;
		}
		return cost;
	};
	PSOvalueType optHead = PSO(headScore, minVs, maxVs, 50, 2, 500);
	getOptData(optHead);
}


void testFindTVD(const std::vector<Eigen::Vector3d>& pC, const std::vector<Eigen::Vector4d>& pMD, double TVD)
{
	int id = findTVD(pC, TVD);
	if (id < 0) {
		id = -id - 2;
	}
	std::cout << "Index of Points is: " << id << "\nTVD in dat Point: " << pC[id][2] << "\n";
	std::pair<double, double> incAzi = CartesianToSpherical(Eigen::Vector3d{pMD[id][1],pMD[id][2],pMD[id][3]});
	std::cout << "Inc,Azi in dat Point: " << incAzi.first << "," << incAzi.second << std::endl;
	int id2 = 0;
	for (size_t i = 0; i < pC.size(); ++i) {
		if (pC[i][2] > TVD) {
			id2 = i-1;
			break;
		}
	}
	std::cout << "Index Real of Point is: " << id2 << "\nTVD in dat Point: " << pC[id2][2] << "\n";
	incAzi = CartesianToSpherical(Eigen::Vector3d{ pMD[id2][1],pMD[id2][2],pMD[id2][3] });
	std::cout << "Inc,Azi in dat Point: " << incAzi.first << "," << incAzi.second << std::endl;

}

void testPenaltyConstraints(const std::vector<Constraint>& cs, std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD) 
{
	std::cout << "Constraints are:\n";
	for (auto x : cs) {
		std::cout << "TVD" << x.lMin.TVD << "\ninc: [" << x.lMin.theta << ", " << x.lMax.theta << "]\nazi: [" << x.lMin.phi << "," << x.lMax.phi << "]  OR ["
			<< 180 + x.lMin.phi << "," << 180 + x.lMax.phi << "]\n_______________________________________________________________\n";
	}
	std::cout << "Penalty is: " << PenaltyConstraint(cs, pC, pMD, 100);

}

void testWellCHCH(const Eigen::VectorXd& x, const Eigen::Vector3d&pinit, const Eigen::Vector3d& target,
					std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD)
{
	std::vector<TrajectoryTemplate*> currWell = wellCHCH(x, pinit, target);
	std::string cond = solve(currWell) == 0 ?"Good": "Bad";
	std::cout << "Condition of Well is: " << cond << "\n";
	std::cout << "Parameters:\n" << "Hold: " << x[0] << "\nDLS1: " << x[1] << "\nDLS2: " << x[2] << "\ninc,azi in Target: " << x[3] << "," << x[4] << std::endl;
	std::cout << "Length: " << allLength(currWell);
	pC = allPointsCartesian(currWell);
	pMD = allPointsMD(currWell);
	writeDataCartesian(pC, "WellCHCHtest.txt");
	writeInclinometry(pMD, "WellCHCHtestInc.txt");
}

void testCH(double dls, const Eigen::Vector3d& pinit, const Eigen::Vector3d& target)
{
	double R = dls < EPSILON ? 1 / EPSILON : 1800 / PI / dls;
	CurveHold ch(pinit, target, 0, 0, R);
	ch.fit();
	ch.points(CoordinateSystem::CARTESIAN);
	std::cout << "alpha is: " << ch.getAlpha() << "\n";
	std::cout << "tangent is:\n " << ch.getTangent2() << "\n";
	writeDataCartesian(ch.pointsCartesian, "WellCHtest.txt");
}

void testScoreCHCH(const Eigen::VectorXd& x,const Eigen::Vector3d& pinit, const Eigen::Vector3d& target, const std::vector<Constraint>& cs,
					const std::vector<double>& minValues, const std::vector<double>& maxValues,bool argKnow) 
{
	std::vector<std::vector<Eigen::Vector3d>> pCWells;
	std::vector<std::vector<Eigen::Vector4d>> pMDWells;

	std::function<double(const Eigen::VectorXd&)> scoreCHCH = [&](const Eigen::VectorXd& x) {
		std::vector<TrajectoryTemplate*> tmp = wellCHCH(x, pinit, target);
		int s = solve(tmp);
		double cost = orderScore1(tmp, pCWells, pMDWells), dlsPen = PenaltyDLS(tmp, 500), length = allLength(tmp);
		std::vector<Eigen::Vector3d> pCtmp = allPointsCartesian(tmp);
		std::vector<Eigen::Vector4d>  pMDtmp = allPointsMD(tmp);
		double penConstraints = PenaltyConstraint(cs, pCtmp, pMDtmp, 100);
		for (auto x : tmp) {
			delete x;
		}
		return cost + penConstraints + PenaltyLength(length, 100);
	};
	std::random_device rd;
	std::mt19937 gen(rd());
	std::vector<std::uniform_real_distribution<double>> distr;
	for (size_t i = 0; i < minValues.size(); ++i) {
		distr.push_back(std::uniform_real_distribution(minValues[i], maxValues[i]));
	}
	double tmpPen = 0;
	Eigen::VectorXd y{ {0,0,0,0,0} };
	if (!argKnow)
	{
		while (true) {
			for (size_t i = 0; i < minValues.size(); ++i)
			{
				y[i] = distr[i](gen);
				std::cout << y[i] << ",";
			}
			tmpPen = scoreCHCH(y);
		}
	}
	else
	{
		double sc = scoreCHCH(x);
	}
}