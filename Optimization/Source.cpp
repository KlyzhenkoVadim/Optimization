#include "Eigen/Dense"
#include <iostream>
#include <vector>
#include "TrajectoryTemplate.h"

class Calculator {
public:
	virtual void calc() = 0;
};
class Calculater : public Calculator {
public:
	void calc() override {
		std::cout << "Calculater " << std::endl;
	}
};
class Calculatir : public Calculator {
public:
	void calc() override {
		std::cout << "Calculatir " << std::endl;

	}
};


int main()
{
	auto isSegmentIntersect = [](const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, const Eigen::Vector3d& p4) {

		auto direction = [](const Eigen::Vector3d& pi, const Eigen::Vector3d& pj, const Eigen::Vector3d& pk) {
			Eigen::Vector3d crossProd = (pk - pi).cross(pj - pi);
			return crossProd[2];
		};

		auto onSegment = [](const Eigen::Vector3d& pi, const Eigen::Vector3d& pj, const Eigen::Vector3d& pk) {
			if (((std::min(pi[0], pj[0]) <= pk[0]) and (std::max(pi[0], pj[0]) >= pk[0])) and
				((std::min(pi[1], pj[1]) <= pk[1]) and (std::max(pi[1], pj[1]) >= pk[1])) and
				((std::min(pi[2], pj[2]) <= pk[2]) and (std::max(pi[2], pj[2]) >= pk[2]))) {
				return true;
			}
			return false;
		};

		int d1 = direction(p3, p4, p1);
		int d2 = direction(p3, p4, p2);
		int d3 = direction(p1, p2, p3);
		int d4 = direction(p1, p2, p4);

		if (((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0))) {
			return true;
		}
		else if (d1 == 0 and onSegment(p3, p4, p1) == true) {
			return true;
		}
		else if (d2 == 0 and onSegment(p3, p4, p2) == true) {
			return true;
		}
		else if (d3 == 0 and onSegment(p1, p2, p3) == true) {
			return true;
		}
		else if (d4 == 0 and onSegment(p1, p2, p4) == true) {
			return true;
		}
		else {
			return false;
		}
	};

	//Eigen::MatrixXd A{ {1., 2.}, {3., 4.} };
	//Eigen::MatrixXd B{ {1., 2.}, {3., 4.} };
	//std::cout << A * B << '\n';

	Eigen::Vector3d p1{ 3,1,0 };
	Eigen::Vector3d p2{ 3,3,4 };
	Eigen::Vector3d p3{ 1,1,0 };
	Eigen::Vector3d p4{ 4,2,0 };

	bool flag = isSegmentIntersect(p1, p2, p3, p4);

	std::cout << (p1 - p2).cross(p3 - p4);


	//std::cout << -p1;

	Calculator* c1 = new Calculater();
	Calculator* c2 = new Calculatir();

	std::vector<Calculator*> calcs;

	calcs.push_back(c1);
	calcs.push_back(c2);

	calcs[0]->calc();
	calcs[1]->calc();



	return 0;
}

