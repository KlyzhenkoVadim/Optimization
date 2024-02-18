#include "OptimizeWells.h"
#include "trajectory_test.h"
#include "optimize_wells_test.h"
#include "Separation_factor.h"
#include <nlohmann/json.hpp>
#include <fstream>

typedef Eigen::Vector3d Vec3d;

void well2_example()
{
	Vec3d pinit{ 0,0,0 };
	Vec3d pt1{ 3000,0,2000 };
	Vec3d pt3{ 3800,0,2000 };
	std::vector<double> x{0.1, 1., 1.2, 1, 1, 0., 0, 1.5, 0.5};
	auto w = well2CHCH(x, pinit, pt1, pt3);
	size_t cond = solve(w);
	if (cond)
	{
		std::cout << "error";
		return;
	}
	auto points = allPointsCartesian(w);
	nlohmann::json j;
	j["points"] = points;
	j["VHold"] = x[0] * 900. + 100.;
	j["dls1"] = x[1];
	j["dls2"] = x[2];
	j["dls3"] = x[3];
	j["dls4"] = x[4];
	j["inc"] = 180. * x[5];
	j["azi"] = 360. * x[6];
	j["NS"] = pt1[0] + 1000. * (1. - 2. * x[7]);
	j["EW"] = pt1[1] + 1000. * (1. - 2. * x[8]);

	writeGG("../output/well2_example.json", j);

}

void testCase1()
{
	std::ifstream ifs("../input/optimization_input.json");
	nlohmann::json jInput = nlohmann::json::parse(ifs)["wells"];
	size_t size = jInput.size();
	auto j0 = jInput[0];
	std::pair<double, double> dls = {0., j0["DLS"].get<double>() };
    std::pair<double, double> firstHold =
        j0["VHold"].get<std::pair<double, double>>();
    std::pair<double, double> lastHold =
        j0["Final_hold"].get<std::pair<double, double>>();
	std::vector<Vec3d> pinits;
	std::vector<Vec3d> targets1;
	std::vector<Vec3d> targets3;
	
	pinits.resize(size);
	targets1.resize(size);
	targets3.resize(size);

	for (int i = 0; i < size; ++i)
	{
		const auto& j = jInput[i];
		pinits[i] << j["x0"][0].get<double>(), j["x0"][1].get<double>(),0.;
		targets1[i] << j["T1"][0].get<double>(),
			j["T1"][1].get<double>(), j["T1"][2].get<double>();
		targets3[i] << j["T3"][0].get<double>(),
			j["T3"][1].get<double>(), j["T3"][2].get<double>();
	}
	optimizeHorizontals(pinits, targets1, targets3);
	//OptimizeTogether(pinits, targets1, targets3);
}

int main()
{
	/*Vec3d pInit{ 5803198,684790,0 }, pT1{ 500,0,1000}, pT3{ 1000,50,1000 },target40R = { 5803236,682857,2900 }, 
		target4001 = { 5803529,682498,2900 }, target4003 = { 5803536,683257,2900 }, target4000 = { 5803409,683700,2900 };
	std::vector<Vec3d> targets3 = { target4001,target4003 }, targets1 = { target40R,target4000 };
	OptimizeHorizontals(pInit, targets1, targets3);*/
	testCase1();
	//aloneOpt();
	//well2_example();
	//checkCorrectnessOfUpperBound();

	return 0;
}
