#pragma once
#include "TrajectoryTemplate.h"
#include <cmath>
#include <iostream>
#include "Penalties.h"
/**
 * @brief ������� ��� ���������� ������� ���������� ���� ����������.
 * (separation factor). �� ���������� ������ travelling cyllinder(!!!)
 * @param pCartesianW1 - ����� ��������� (NS,EW,TVD) ������ ���������� (��������)
 * @param pMDW1 - ����� md � ��������� ������������ ������� (MD,NSt,EWt,TVDt)
 *  ������ ���������� (��������)
 * @param pCartesianW2 - ����� ��������� (NS,EW,TVD) ������ ���������� (������)
 * @param pMDW2 - ����� md � ��������� ������������ ������� (MD,NSt,EWt,TVDt)
 * ������ ���������� (������)
 * @param TVDstart - �������� �������, ������� � ������� ����������� ������ ����������
 * (������������� ��� ������ ������ ������������� �������.
 * @param actFunc - ����, ���� true - ����������� �������� �� ������� ����������,
 * ���� false - ����������� ��� ������ ����������.
 * @param penalty - �����, ���� ������ ���������� ������ 1.5.
 * @return ������ ����������, ��� �������� �������� �� ����.
*/
double sepFactor(const std::vector<Eigen::Vector3d>& pCartesianW1,
                 const std::vector<Eigen::Vector4d>& pMDW1,
                 const std::vector<Eigen::Vector3d>& pCartesianW2,
                 const std::vector<Eigen::Vector4d>& pMDW2, double TVDstart = 0,
                 bool actFunc = true, double penalty = 20);
/**
 * @brief ���������� ������� ��������� ������� ���������� (ddi)
 * @param well - ������ TrajectoryTemplate* ����������� ����������. 
 * @param pCartesian - ����������������� ����� �� ����������
 * @param actFunc - ���� true, ������������ �������� �� ddi. 
 * @param penalty - ����� �� ���������� ddi �������� 6.25
 * @return � ����������� �� actFunc (ddi ��� �������� �� ����)
*/
double DDI(std::vector<TrajectoryTemplate*>& well,
           const std::vector<Eigen::Vector3d>& pCartesian, bool actFunc = true,
           double penalty = 5);
/**
 * @brief ���������� ERD (Extended Reach Drilling) ��� �������������� �������.
 * @param points - ������ ����������������� ����� ����������
 * @return �������� �������.
*/
double ERD(const std::vector<Eigen::Vector3d>& points);
/**
 * @brief ���������� �������� ����������� ������� ���������� 
 * ��� ������� ��� �����������!
 * @param mainWell - ������ TrajectoryTemplate* �������������� ����������.
 * @param pCTrajectories - ����� �������� ����� ��� ����������� ����������.
 * @param pMDTrajectories - ����� ����������� �������� � md ����������� ����������.
 * @param SepFactorShift - �������� TVDStart ��� ������� SepFactor.
 * @param penalty - ����� �� ������������.
 * @return ������� �������: ������������� ����� + �������� sepFactor + �������� DDI (ERD)
*/
double orderScore1(
    std::vector<TrajectoryTemplate*>& mainWell,
    const std::vector<std::vector<Eigen::Vector3d>>& pCTrajectories,
    const std::vector<std::vector<Eigen::Vector4d>>& pMDTrajectories,
    double SepFactorShift = 0, double penalty = 1000);
/**
* @TODO �������, �� ������������
*/
double scoreSolver(std::vector<TrajectoryTemplate*>& tmp,
                   const WellTrajectoryConstraints& cs, double penalty = 1000);
/**
 * @brief ������� ��� ����������� dls
 * @param tangent1 ����������� ������ ������ ����
 * @param tangent2 ����������� ������ ����� ����
 * @return dog leg severity 
*/
double dls(Eigen::Vector3d& tangent1, Eigen::Vector3d& tangent2);
/**
 * @brief ���������� ���������� tortuosity ����������.
 * @param well ������ TrajectoryTemplate* ����������.
 * @return 
*/
double tortuosity(std::vector<TrajectoryTemplate*>& well);
