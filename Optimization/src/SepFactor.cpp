#include "./SepFactor.h"

double SeparationFactor(const std::vector<Eigen::Vector3d>& pC1,
	const std::vector<Eigen::Vector4d>& pMD1,
	const std::vector<Eigen::Vector3d>& pC2,
	const std::vector<Eigen::Vector4d>& pMD2,
	double TVDShift)
{
	double sf = std::numeric_limits<double>::max();
	double sgm = 9e-3;
	int w1_start = 1, w2_start = 1;
	if (TVDShift > EPSILON)
	{
		w1_start = 0;
		w2_start = 0;
		for (auto x : pC1)
		{
			if (x[2] > TVDShift)
				break;
			w1_start += 1;
		}
		for (auto x : pC2)
		{
			if(x[2] > TVDShift)
				break;
			w2_start += 1;
		}
	}

	for (size_t i = w1_start; i < pC1.size(); ++i)
	{
		for (size_t j = w2_start; j < pC2.size(); ++j)
		{
			sf = std::min(sf, (pC1[i] - pC2[j]).norm() / sgm/(pMD1[i][0] + pMD2[j][0]));
		}
	}
	return sf;
}


double SeparationFactorNewton(std::vector<TrajectoryTemplate*>& well1, std::vector<TrajectoryTemplate*>& well2)
{
	double sgm = 9e-3;
	double L1 = allLength(well1), L2 = allLength(well2);
	auto g = [&](double t, double s)
	{
		return (FunctionWellPoint(t, well1) - FunctionWellPoint(s, well2)).norm() / (sgm * (L1 * t + L2 * s));
	};
	// TODO: Доделать до конца, и сравнить насколько быстро сходится
}


void gabbage_collector(std::vector<TrajectoryTemplate*>& well)
{
	for (auto x : well)
		delete x;
	well.clear();
}

void set_points(std::vector<TrajectoryTemplate*>& well, std::vector<Eigen::Vector3d>& pC, std::vector<Eigen::Vector4d>& pMD)
{
	pC.clear();
	pMD.clear();
	pC = allPointsCartesian(well);
	pMD = allPointsMD(well);
	gabbage_collector(well);
}



void test_separationFactor()
{
	std::vector<TrajectoryTemplate*> well1 = { new CurveHold(Eigen::Vector3d{0,0,0},Eigen::Vector3d{-100,0,1000},0,0,300,100) },
		well2 = { new CurveHold(Eigen::Vector3d{1,0,0},Eigen::Vector3d{100,0,1000},0,0,300,100) };
	std::vector<Eigen::Vector3d> pC1, pC2;
	std::vector < Eigen::Vector4d> pMD1, pMD2;
	set_points(well1, pC1, pMD1);
	set_points(well2, pC2, pMD2);
	double sfactor = SeparationFactor(pC1, pMD1, pC2, pMD2);
	double sf12 = sepFactor(pC1, pMD1, pC2, pMD2, 0, false);
	double sf21 = sepFactor(pC2, pMD2, pC1, pMD1, 0, false);

	assert(fabs(sfactor - std::min(sf12,sf21)) < EPSILON);
	double sf22 = sepFactor(pC1, pMD1, pC1, pMD1, 0, false);
	assert(sf22 < EPSILON);
	double sf11 = SeparationFactor(pC1, pMD1, pC1, pMD1);
	assert(sf11 < EPSILON);
	well1.push_back(new Hold(Eigen::Vector3d{ -10,0,0 }, Eigen::Vector3d{ 10,0,100 },1000));
	well2.push_back(new Hold(Eigen::Vector3d{ 10,0,0 }, Eigen::Vector3d{ -10,0,100 }, 1000));
	set_points(well1, pC1, pMD1);
	set_points(well2, pC2, pMD2);
	sfactor = SeparationFactor(pC1, pMD1, pC2, pMD2);
	assert(sfactor - 1.5 < EPSILON);
	sf12 = sepFactor(pC1, pMD1, pC2, pMD2, 0, false);
	sf21 = sepFactor(pC2, pMD2, pC1, pMD1, 0, false);
	assert(std::min(sf12,sf21) - 1.5 < EPSILON);

	
}



void test_sepfactor_speed()
{
	std::vector<TrajectoryTemplate*> well1 = { new CurveHold(Eigen::Vector3d{0,0,0},Eigen::Vector3d{-100,0,1000},0,0,300,500) },
		well2 = { new CurveHold(Eigen::Vector3d{1,0,0},Eigen::Vector3d{100,0,1000},0,0,300,500) };
	std::vector<Eigen::Vector3d> pC1 = allPointsCartesian(well1), pC2 = allPointsCartesian(well2);
	std::vector < Eigen::Vector4d> pMD1 = allPointsMD(well1), pMD2 = allPointsMD(well2);

	double sfactor, sf12, sf21;
	double start = time(NULL), end = start;
	int numiter = 0;
	while (end - start - 30 < EPSILON)
	{
		sf12 = sepFactor(pC1, pMD1, pC2, pMD2, 0, false);
		sf21 = sepFactor(pC2, pMD2, pC1, pMD1, 0, false);

		end = time(NULL);
		numiter += 1;
	}
	std::cout << "number of iterations sepFactor: " << numiter << "Value: " << std::min(sf12,sf21) << std::endl;
	start = time(NULL);
	end = start;
	numiter = 0;
	while (end - start - 30 < EPSILON)
	{
		sfactor = SeparationFactor(pC1, pMD1, pC2, pMD2,0);
		end = time(NULL);
		numiter += 1;
	}
	std::cout << "number of iterations SeparationFactor:" << numiter << "Value: " << sfactor << std::endl;
}