#include "trajectory_test.h"

#include "TrajectoryTemplate.h"
#include "Hold.h"
#include "Curve.h"
#include "CurveHold.h"
#include "CurveHoldCurveHold.h"
#include "TestWells.h"

void printvec(const Eigen::Vector3d& v)
{
	std::cout << "[ ";
	for (auto x : v)
		std::cout << x << " ";
	std::cout << "]";
}

void well_deleter(std::vector<TrajectoryTemplate*>& well)
{
	for (auto x : well)
		delete x;
	well.clear();
}

void writeforcheck(std::vector<Eigen::Vector3d>& pFromFunc, std::vector<Eigen::Vector3d>& pC)
{
	writeDataCartesian(pFromFunc, "../output/Cartesian/fromfunction.txt");
	writeDataCartesian(pC, "../output/Cartesian/pC.txt");

}

double GetMinimalDistance(std::vector<Eigen::Vector3d>& p_from_func, std::vector<Eigen::Vector3d>& pC)
{
	double distance = 0;
	double dpoint = std::numeric_limits<double>::max();
	for (size_t i = 0; i < p_from_func.size(); ++i)
	{
		for (size_t j = 0; j < pC.size(); ++j)
		{
			dpoint = std::min((p_from_func[i] - pC[j]).norm(), dpoint);
		}
		distance += dpoint;
		dpoint = std::numeric_limits<double>::max();
	}
	return distance/p_from_func.size();
}

void SetPointsFromFunc(std::vector<Eigen::Vector3d>& points_from_func,std::vector<TrajectoryTemplate*>& well, size_t N)
{
	for(size_t i = 0; i < N; ++i )
		points_from_func.push_back(FunctionWellPoint(double(i) / (N - 1), well));
}

void checkAnswer(std::vector<TrajectoryTemplate*>& well)
{
	assert(solve(well) == 0);
	std::vector<Eigen::Vector3d> pC = allPointsCartesian(well), points_from_function;
	SetPointsFromFunc(points_from_function, well, pC.size());
	double dist = 1. / 2. * (GetMinimalDistance(pC, points_from_function) + GetMinimalDistance(points_from_function, pC));
	double n = allLength(well) / pC.size();
	assert(dist < n);
	std::cout << "good\n";
	well_deleter(well);
	pC.clear();
	points_from_function.clear();
}

void test_trajectory_pspeed()
{
	std::vector<Eigen::Vector3d> pfunc;
	std::vector<Eigen::Vector3d> pC;
	size_t nums = 500;
	std::vector<TrajectoryTemplate*> well = { new CurveHoldCurveHold(Eigen::Vector3d{ 0,0,0 }, 0, 0, 300, 300,
		Eigen::Vector3d{ 500,100,1500 }, Eigen::Vector3d{ 1000,1000,1500 }, 0.001, nums)};
	well_deleter(well);
	well = { new CurveHold(Eigen::Vector3d{0,0,0},Eigen::Vector3d{100,100,1000},0,0,300,500) };
	double start = time(NULL), end = start;
	int num = 0;
	while (end - start - 30 < EPSILON)
	{
		//pC = allPointsCartesian(well);
		SetPointsFromFunc(pfunc, well, nums);
		num += 1;
		end = time(NULL);
	}
	std::cout << num;
}

void test_trajectory_points()
{
	//Hold
	std::vector<TrajectoryTemplate*> well = { new Hold(Eigen::Vector3d{0,0,0},Eigen::Vector3d{0,0,100},100) };
	std::vector<Eigen::Vector3d> pC = allPointsCartesian(well);
	double s = 0;
	for (size_t i = 0; i < pC.size(); ++i)
	{
		assert((FunctionWellPoint(double(i) / (pC.size()-1), well) - pC[i]).norm() < EPSILON);
	}
	well_deleter(well);
	pC.clear();
	//Curve
	well = { new Curve(Eigen::Vector3d{0,0,0},0,0,90,0,300,TypeCurve::TVD,100) };
	pC = allPointsCartesian(well);
	for (size_t i = 0; i < pC.size(); ++i)
	{
		assert((FunctionWellPoint(double(i) / (pC.size() - 1), well) - pC[i]).norm() < EPSILON);
	}
	well_deleter(well);
	pC.clear();
	//CurveHold
	well = { new CurveHold(Eigen::Vector3d{0,0,0},Eigen::Vector3d{100,100,1000},0,0,300,500) };
	pC = allPointsCartesian(well);
	std::vector<Eigen::Vector3d> points_from_function;
	for (size_t i = 0; i < pC.size(); ++i)
	{
		points_from_function.push_back(FunctionWellPoint(double(i) / (pC.size()-1), well));
	}
	double dist = 1. / 2. * (GetMinimalDistance(points_from_function, pC) + GetMinimalDistance(pC, points_from_function));
	assert(dist < 1);
	well_deleter(well);
	pC.clear();
	points_from_function.clear();
	// CurveHoldCurveHold
	well = { new CurveHoldCurveHold(Eigen::Vector3d{0,0,0},0,0,300,300,Eigen::Vector3d{500,100,1500},Eigen::Vector3d{1000,1000,1500},0.001,500) };
	pC = allPointsCartesian(well);
	for (size_t i = 0; i < pC.size(); ++i)
	{
		points_from_function.push_back(FunctionWellPoint(double(i) / (pC.size() - 1), well));
	}
	dist = 1. / 2. * (GetMinimalDistance(points_from_function, pC) + GetMinimalDistance(pC, points_from_function));
	int n = well[0]->length() / pC.size();
	assert(dist < n );
	well_deleter(well);
	pC.clear();
	points_from_function.clear();

	// Hold - CurveHold
	well = { new Hold(Eigen::Vector3d{0,0,0},Eigen::Vector3d{0,0,400},20), new CurveHold(Eigen::Vector3d{0,0,400},Eigen::Vector3d{500,100,1500},0,0,300,100) };
	pC = allPointsCartesian(well);
	for (size_t i = 0; i < pC.size(); ++i)
	{
		points_from_function.push_back(FunctionWellPoint(double(i) / (pC.size() - 1), well));
	}
	dist = 1. / 2 * (GetMinimalDistance(points_from_function, pC) + GetMinimalDistance(pC, points_from_function));
	n = allLength(well) / pC.size();
	assert(dist < n);
	well_deleter(well);
	pC.clear();
	points_from_function.clear();
	// Hold - Curve
	well = { new Hold(Eigen::Vector3d{0,0,0},Eigen::Vector3d{0,0,400},20), new Curve(Eigen::Vector3d{0,0,400},0,0,90,10,1000,TypeCurve::TVD,100) };
	pC = allPointsCartesian(well);
	for (size_t i = 0; i < pC.size(); ++i)
	{
		points_from_function.push_back(FunctionWellPoint(double(i) / (pC.size() - 1), well));
	}
	n = allLength(well) / pC.size();
	dist = 1. / 2 * (GetMinimalDistance(pC, points_from_function) + GetMinimalDistance(points_from_function, pC));
	assert(dist < n);
	well_deleter(well);
	pC.clear();
	points_from_function.clear();
	// Hold-CurveHoldCurveHold
	well = { new Hold(Eigen::Vector3d{0,0,0},Eigen::Vector3d{0,0,400},20), 
		new CurveHoldCurveHold(Eigen::Vector3d{0,0,400},0,0,400,400,
			Eigen::Vector3d{100,500,1500},Eigen::Vector3d{1000,500,1500}) };
	pC = allPointsCartesian(well);
	SetPointsFromFunc(points_from_function, well, pC.size());
	dist = 1. / 2 * (GetMinimalDistance(pC, points_from_function) + GetMinimalDistance(points_from_function, pC));
	n = allLength(well) / pC.size();
	assert(dist < n);
	well_deleter(well);
	pC.clear();
	points_from_function.clear();
	// Curve - Hold
	well = { new Curve(Eigen::Vector3d{0,0,0},0,0,90,0,500,TypeCurve::TVD) };
	well.back()->getTarget1Point();
	well.push_back(new Hold(well.back()->pointT1, 90, 0, 1000, typeHold::md, 25));
	pC = allPointsCartesian(well);
	SetPointsFromFunc(points_from_function, well, pC.size());
	dist = 1. / 2. * (GetMinimalDistance(pC, points_from_function) + GetMinimalDistance(points_from_function, pC));
	n = allLength(well) / pC.size();
	assert(dist < n);
	well_deleter(well);
	pC.clear();
	points_from_function.clear();
	//Curve-CurveHold;
	well = { new Curve(Eigen::Vector3d{0,0,0},0,0,90,0,300) };
	well[0]->getTarget1Point();
	well.push_back(new CurveHold(well[0]->pointT1, Eigen::Vector3d{ 700,500,1000 }, 90, 0, 300));
	checkAnswer(well);
	// Curve-CurveHoldCurveHold
	well = { new Curve(Eigen::Vector3d{0,0,0},0,0,90,0,300,TypeCurve::DLS) };
	well[0]->getTarget1Point();
	well.push_back(new CurveHoldCurveHold(well[0]->pointT1, 90, 0, 400, 400, Eigen::Vector3d{ 700,700,1700 }, Eigen::Vector3d{ 1000,500,1700 }));
	checkAnswer(well);
	// CurveHold - Hold
	CurveHold ch(Eigen::Vector3d{ 0,0,0 }, Eigen::Vector3d{ 500,600,600 }, 0, 0, 300);
	auto [inc, azi] = CartesianToSpherical(ch.getTangent2());
	well = {new CurveHold(Eigen::Vector3d{ 0,0,0 }, Eigen::Vector3d{ 500,600,600 }, 0, 0, 300)
		, new Hold(Eigen::Vector3d{500,600,600},inc,azi,100) };
	checkAnswer(well);
	// CurveHold - Curve
	well = { new CurveHold(Eigen::Vector3d{ 0,0,0 }, Eigen::Vector3d{ 500,600,600 }, 0, 0, 300)
		,new Curve(Eigen::Vector3d{500,600,600},inc,azi,90,0,1000,TypeCurve::TVD) };
	checkAnswer(well);
	// CurveHold - CurveHoldCurveHold
	well = { new CurveHold(Eigen::Vector3d{ 0,0,0 }, Eigen::Vector3d{ 500,600,600 }, 0, 0, 300)
		,new CurveHoldCurveHold(Eigen::Vector3d{500,600,600},inc,azi,300,300,
			Eigen::Vector3d{1500,1000,1800}, Eigen::Vector3d{2500,1000,1800}) };
	checkAnswer(well);
	// CHCH - Hold
	well = { new CurveHoldCurveHold(Eigen::Vector3d{0,0,0},0,0,300,300,Eigen::Vector3d{100,600,1200},90,0,700),
			new Hold(Eigen::Vector3d{100,600,1200},90,0,1500) };
	checkAnswer(well);
	//CHCH - Curve
	well = { new CurveHoldCurveHold(Eigen::Vector3d{0,0,0},0,0,300,300,Eigen::Vector3d{100,600,1200},90,0,700),
			new Curve(Eigen::Vector3d{100,600,1200},90,0,0,0,1700,TypeCurve::TVD) };
	checkAnswer(well);
	//CHCH - CurveHold
	well = { new CurveHoldCurveHold(Eigen::Vector3d{ 0,0,0 }, 0, 0, 300, 300, Eigen::Vector3d{ 100,600,1200 }, 90, 0, 700),
			new CurveHold(Eigen::Vector3d{100,600,1200},Eigen::Vector3d{1000,1000,1800},90,0,300) };
	checkAnswer(well);
}
