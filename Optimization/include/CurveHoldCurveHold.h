#ifndef CURVEHOLDCURVEHOLD_H_
#define CURVEHOLDCURVEHOLD_H_
#pragma once
#include "TrajectoryTemplate.h"
#include "CurveHold.h"

class CurveHoldCurveHold : public TrajectoryTemplate
{
private:
	Eigen::Vector3d p1; // ���������� ��������� ����� ���������� (���������� N,E,V)
	Eigen::Vector3d p4; // ���������� �������� ����� ���������� (���������� N,E,V)
	Eigen::Vector3d t1; // ��������� ����������� ������ ��������� �����
	Eigen::Vector3d t4;	// ��������� ����������� ������ �������� �����
	Eigen::Vector3d pInter;

	Eigen::Vector3d r1; // ���������� ������ ���������� ������� ������� curve
	Eigen::Vector3d r4; // ���������� ������ ���������� ������� ������� curve

	double tetta1; // �������� ��������� ���� (� ��������) ������������ ������� t1 ��������� �����
	double phi1; // �������� ������������� ���� (� ��������) ������������ ������� t1 ��������� �����
	double tetta4; // �������� ��������� ���� (� ��������) ������������ ������� t1 �������� �����
	double phi4; // �������� ������������� ���� (� ��������) ������������ ������� t1 �������� �����

	double R1; // ������ �������� ���� � p1
	double R2; // ������ �������� ���� � p4
	size_t nums;// ����� ����� ����������


	double betta; // ����� ���������� ������� ������������ hold, ���� betta = 0. �� ������ CurveHoldCurve
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
