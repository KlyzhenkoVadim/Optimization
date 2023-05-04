#include "OptimizeWells.h"
#include "trajectory_test.h"
#include "optimize_wells_test.h"
#include "Separation_factor.h"

typedef Eigen::Vector3d Vector3d;

int main()
{
	
	Eigen::Vector3d pInit{ 5803198,684790,0 }, pT1{ 500,0,1000}, pT3{ 1000,50,1000 },target40R = { 5803236,682857,2900 }, 
		target4001 = { 5803529,682498,2900 }, target4003 = { 5803536,683257,2900 }, target4000 = { 5803409,683700,2900 };
	std::vector<Eigen::Vector3d> targets3 = { target4001,target4003 }, targets1 = { target40R,target4000 };
	//test_WellTrajectorySolver();
	targets3 = { target4003,target4001 };
	targets1 = { target4000,target40R };

	/*gBestCost: 2.97501
	gBestPos: 410.3396, 1.001136, 1.04123, 1.472675, 1.451098, 41.65297, 276.3627, 0.5348235, 0.4225081,
	gBestCost: 2.891164
	gBestPos: 863.9466, 0.7543149, 1.017725, 1.323914, 1.083109, 51.04885, 281.5318, 0.5132608, 0.4064161,*/
	//OptimizeHorizontal(pInit, target40R, target4001);
	//OptimizeHorizontal(pInit, target4000, target4001);
	OptimizeHorizontals(pInit, targets1, targets3);
	return 0;
}
