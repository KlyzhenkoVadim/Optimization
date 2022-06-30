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

	std::vector<Eigen::Vector3d> pointsCartesian;
	std::vector<Eigen::Vector4d> pointsMD;
	double alpha;
	double betta;

public:
	CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3, double teta, double phi, double R, size_t nums = 100);
	void fit() override;
	void points(CoordinateSystem coordinateSystem) override;
	double length() override;
	double getAlpha();
	
	Eigen::Vector3d getInitPoint();
	Eigen::Vector3d getTargetPoint();
};

