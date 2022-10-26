#include "optimize_wells_test.h"

void brute_force_minimum(int dim, std::vector<double>& values, double* currMin, Eigen::VectorXd& xBest,
	std::function<double(const Eigen::VectorXd& x)> cost, const int& max_iter,const int& max_dim, const std::vector<double>& minValues,
	const std::vector<double>& maxValues)
{

	for (size_t i = 0; i < max_iter; ++i)
	{
		if (dim > 0)
		{
			values[dim - 1] = (minValues[dim - 1]) * (1 - double(i) / (max_iter - 1)) + maxValues[dim - 1] * (double(i) / (max_iter - 1));
			brute_force_minimum(dim - 1, values, currMin, xBest, cost,max_iter,max_dim,minValues,maxValues);
		}
		else
		{
			Eigen::VectorXd v(max_dim);
			for (size_t j = 0; j < max_dim; ++j)
				v[j] = values[j];
			double tmpcost = cost(v);
			if (tmpcost < *currMin)
			{
				*currMin = tmpcost;
				xBest = v;
			}
		}
	}
}

void horizontal_optimize_test()
{
	Eigen::Vector3d pinit{ 5803198,684790,0 }, target1{ 5803236,682857,2900 }, target3{ 5803529,682498,2900 };
	std::vector<double> minValues{ 400,0,0,0,0,35,0,0,0 }, maxValues{ target1[2],1.5,1.5,1.5,1.5,55,360,1,1 };
	
	auto score = [&](const Eigen::VectorXd& x)
	{
		std::vector<TrajectoryTemplate*> well = well2CHCH(x,pinit,target1,target3);
		double score_value = OneWellScore(well);
		return score_value;
	};
	int max_iter = 5, max_dim = 9;
	int dim = max_dim;
	std::vector<double> value = { 400,0,0,0,0,35,0,0,0 };
	double min = std::numeric_limits<double>::max();
	Eigen::VectorXd xBest(9);
	brute_force_minimum(dim, value,&min,xBest, score, max_iter, max_dim, minValues, maxValues);
	std::cout << min << "\n";
	for (auto x : xBest)
		std::cout << x << " ";
	std::cout << "\n";
	OptimizeHorizontal(pinit, target1, target3);
}

void ddi_test()
{
	std::vector<TrajectoryTemplate*> well_for_test = { new Hold(Eigen::Vector3d{0,0,0},0,0,100),
		new Curve(Eigen::Vector3d{ 0,0,0}, 0, 0, 90, 0, 600, TypeCurve::TVD,100),
		new CurveHold(Eigen::Vector3d{0,0,0},Eigen::Vector3d{100,500,1000},0,0,500,500),
		new CurveHoldCurveHold(Eigen::Vector3d{0,0,0},0,0,300,300,Eigen::Vector3d{500,600,1200},Eigen::Vector3d{1200,600,1200}) };
	std::vector<TrajectoryTemplate*> tmpwell = { well_for_test[0] };
	for (auto x : well_for_test)
	{
		tmpwell = { x };
		std::vector<Eigen::Vector4d> pmd = allPointsMD(tmpwell);
		std::cout << Tortuosity(tmpwell) << "," << Tortuosity_sum(pmd)<<std::endl;
		assert(abs(Tortuosity(tmpwell) - Tortuosity_sum(pmd)) < 1);
		delete x;
	}
	Eigen::VectorXd x{ {412.723, 0.9665206, 0.5035774, 1.452895, 1.071906, 41.73488, 280.6988, 0.5478006, 0.4144362} };
	Eigen::Vector3d pinit{ 5803198,684790,0 }, target1{ 5803236,682857,2900 }, target3{ 5803529,682498,2900 };
	std::vector<TrajectoryTemplate*> well = well2CHCH(x, pinit, target1, target3);
	std::vector<Eigen::Vector4d> pmd = allPointsMD(well);
	assert(abs(Tortuosity(well) - Tortuosity_sum(pmd)) < 1);

}
