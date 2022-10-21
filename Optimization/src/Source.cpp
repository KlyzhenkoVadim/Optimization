#include "OptimizeWells.h"
#include "SepFactor.h"


int main()
{
	
	Eigen::Vector3d pInit{ 5803198,684790,0 }, pT1{ 500,0,1000}, pT3{ 1000,50,1000 },target40R = { 5803236,682857,2900 }, 
		target4001 = { 5803529,682498,2900 }, target4003 = { 5803536,683257,2900 }, target4000 = { 5803409,683700,2900 };
	std::vector<Eigen::Vector3d> targets3 = { target4001,target4003 }, targets1 = { target40R,target4000 };

	test_separationFactor();
	test_sepfactor_speed();

	return 0;
}
