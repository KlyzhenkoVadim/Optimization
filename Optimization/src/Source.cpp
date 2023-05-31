#include "OptimizeWells.h"
#include "trajectory_test.h"
#include "optimize_wells_test.h"
#include "Separation_factor.h"
#include "C:\Users\klyzhenko.vs\Desktop\ArchiveWProjects\optimization\Optimization\packages\nlohmann.json.3.11.2\build\native\include\nlohmann\json.hpp"
#include <fstream>

typedef Eigen::Vector3d Vec3d;

void simpleTemplates()
{

}

void testCase1()
{
	std::ifstream ifs("../input/optimization_input.json");
	nlohmann::json jInput = nlohmann::json::parse(ifs)["wells"];
	size_t size = jInput.size();
	auto j0 = jInput[0];
	std::pair<double, double> dls = {0., j0["DLS"].get<double>() };
	std::pair<double, double> firstHold = j0["VHold"].get<std::pair<double,double>>();
	std::pair<double, double> lastHold = j0["Final_hold"].get<std::pair<double, double>>();
	std::vector<Vec3d> pinits;
	std::vector<Vec3d> targets1;
	std::vector<Vec3d> targets3;
	
	pinits.resize(size);
	targets1.resize(size);
	targets3.resize(size);

	for (int i = 0; i < size; ++i)
	{
		const auto& j = jInput[i];
		pinits[i] << j["x0"][0].get<double>(), j["x0"][1].get<double>(), 0.;
		targets1[i] << j["T1"][0].get<double>(), j["T1"][1].get<double>(), j["T1"][2].get<double>();
		targets3[i] << j["T3"][0].get<double>(), j["T3"][1].get<double>(), j["T3"][2].get<double>();
	}
	OptimizeHorizontals(pinits, targets1, targets3);
	//OptimizeTogether(pinits, targets1, targets3);
}

int main()
{
	/*Vec3d pInit{ 5803198,684790,0 }, pT1{ 500,0,1000}, pT3{ 1000,50,1000 },target40R = { 5803236,682857,2900 }, 
		target4001 = { 5803529,682498,2900 }, target4003 = { 5803536,683257,2900 }, target4000 = { 5803409,683700,2900 };
	std::vector<Vec3d> targets3 = { target4001,target4003 }, targets1 = { target40R,target4000 };
	OptimizeHorizontals(pInit, targets1, targets3);*/
	testCase1();
	return 0;
}
