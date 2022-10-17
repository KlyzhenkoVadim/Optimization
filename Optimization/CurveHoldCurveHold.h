#ifndef CURVEHOLDCURVEHOLD_H_
#define CURVEHOLDCURVEHOLD_H_
#pragma once
#include "TrajectoryTemplate.h"
#include "CurveHold.h"

class CurveHoldCurveHold : public TrajectoryTemplate
{
private:
	Eigen::Vector3d p1; // координаты начальной точки траектории (координаты N,E,V)
	Eigen::Vector3d p4; // координаты конечной точки траектории (координаты N,E,V)
	Eigen::Vector3d t1; // единичный касательный вектор начальной точки
	Eigen::Vector3d t4;	// единичный касательный вектор конечной точки
	Eigen::Vector3d pInter;

	Eigen::Vector3d r1; // координаты центра окружности первого участка curve
	Eigen::Vector3d r4; // координаты центра окружности второго участка curve

	double tetta1; // значение зенитного угла (в градусах) касательного вектора t1 начальной точке
	double phi1; // значение азимутального угла (в градусах) касательного вектора t1 начальной точке
	double tetta4; // значение зенитного угла (в градусах) касательного вектора t1 конечной точке
	double phi4; // значение азимутального угла (в градусах) касательного вектора t1 конечной точке

	double R1; // радиус кривизны дуги к p1
	double R2; // радиус кривизны дуги к p4
	size_t nums;// число точек траектории


	double betta; // длина финального участка стабилизации hold, если betta = 0. то шаблон CurveHoldCurve
	double eps; 

	Eigen::Vector3d t;
	Eigen::Vector3d p1Inter;
	Eigen::Vector3d p4Inter;
	double alpha1;
	double alpha2;
	double holdLength;
	int condition;
	int fit();


public:
	CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1, double phi1, double R1, double R2, const Eigen::Vector3d& pT1,
		const Eigen::Vector3d& pT3, double eps = 10e-4, size_t nums = 100);

	CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1, double phi1, double R1, double R2, const Eigen::Vector3d& p4, 
		double tetta4, double phi4, double betta = 0.0, double eps = 10e-4, size_t nums = 100);

	int getCondition() override;
	void points(CoordinateSystem coordinateSystem) override;
	double length() override;

	void getInitPoint(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
	void getTarget1Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
	void getTarget3Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
};

#endif // CURVEHOLDCURVEHOLD_H_
