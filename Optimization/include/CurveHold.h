#ifndef CURVEHOLD_H_
#define CURVEHOLD_H_
#pragma once
#include "TrajectoryTemplate.h"
class CurveHold : public TrajectoryTemplate
{
private:
	Eigen::Vector3d p1; // ���������� ��������� ����� ���������� (���������� N,E,V)
	Eigen::Vector3d p3; // ���������� �������� ����� ���������� (���������� N,E,V)
	Eigen::Vector3d t1; // ��������� ����������� ������ ��������� �����
	Eigen::Vector3d t2; // ��������� ����������� ������ �������� ����� ����

	// condition??
	

	double teta; // �������� ��������� ���� (� ��������) ������������ ������� t1 ��������� �����
	double phi; // �������� ������������� ���� (� ��������) ������������ ������� t1 ��������� �����
	double R; // ������ �������� ���������� (R>0)
	size_t nums;// ����� ����� ����������

	//std::vector<Eigen::Vector3d> pointsCartesian;
	//std::vector<Eigen::Vector4d> pointsMD;
	double alpha;
	double betta;
	int condition;
	int fit();

public:
	CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3, double teta, double phi, double R, size_t nums = 50);
	CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3, const Eigen::Vector3d& t1, double R, size_t nums = 50);
	//CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3, double teta, double phi, size_t nums = 100);
	int getCondition() override;
	void points(CoordinateSystem coordinateSystem) override;
	double length() override;
	double getAlpha();
	Eigen::Vector3d getPointInterpol();
	Eigen::Vector3d getTangent2();

	void getInitPoint(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
	void getTarget1Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
	void getTarget3Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
};

#endif // CURVEHOLD_H_
