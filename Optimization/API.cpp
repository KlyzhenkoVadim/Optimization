#include "API.h"

Eigen::Vector3d calcTangentVector(double azimuth, double inclination) {

	double x = sin(inclination * PI / 180.0) * cos(azimuth * PI / 180.0);
	double y = sin(inclination * PI / 180.0) * sin(azimuth * PI / 180.0);
	double z = cos(inclination * PI / 180.0);

	return Eigen::Vector3d{ fabs(x) > EPSILON ? x : 0.0,
							 fabs(y) > EPSILON ? y : 0.0,
							 fabs(z) > EPSILON ? z : 0.0
	};
};

someAPI::someAPI(const std::string type) {
	this->type = type;
};

void someAPI::setData() {
	if (type == "hold2chch") {

	}
};

void someAPI::setPSOdata() {};

void someAPI::getData() {};

void someAPI::getPSOdata() {};

void someAPI::Optimize() {};

void someAPI::TypeTrajectory() {
	if (type == "hold2chch") {
		typeWell = [&](Eigen::VectorXd& x) {
			std::vector<TrajectoryTemplate*> well1;
			Eigen::Vector3d pIHold = { x[0],x[1],x[2] };
			Eigen::Vector3d pTHold = { x[3],x[4],x(5) };
			Eigen::Vector3d pT3Chch1 = { x[6],x[7],x[8] };
			Eigen::Vector3d Tangent = calcTangentVector(x[10], x[9]);
			Eigen::Vector3d pT1Chch1 = pT3Chch1 - 110. * Tangent;
			Eigen::Vector3d pT1Chch2 = { x[11],x[12],x[13] };
			Eigen::Vector3d pT3Chch2 = { x[14] ,x[15],x[16] };
			well1.push_back(new Hold(pIHold, pTHold));
			well1.push_back(new CurveHoldCurveHold(pTHold, 0, 0, x[17], x[18], pT1Chch1, pT3Chch1));
			well1.push_back(new CurveHoldCurveHold(pT3Chch1, x[3], x[4], x[16], x[20], pT1Chch2, pT3Chch2));
			return well1;
		};
	}
}



