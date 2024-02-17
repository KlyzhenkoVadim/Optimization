#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include <TrajectoryTemplate.h>

double OneWellScore(std::vector<TrajectoryTemplate*>& mainWell,
                    double penalty = 1000);

/*
* \brief ������� ��� �������� ���������� �� ��������� ������ ����������
* \param x ������ ����������
* \param pinit - 3� ������ ������ ����������
* \param pT1 - 3� ������ ���� �1 ����� ����������
* \param pT3 - 3� ������ ���� �3 ����� ����������
* \param x[0] - ����� �� 0 �� 1, �����, ��� 900 * x[0] - ����� ������������� �������
* \param x[1] - dls ������� ������� �����������
* \param x[2] - dls ������� ������� �����������
* \param x[3] - dls �������� ������� �����������
* \param x[4] - dls ��������� ������� �����������
* \param x[5] - ����� �� 0 �� 1, �����, ��� 180 * x[5] - �������� ���� ����� ������� ChCh
* \param x[6] - ����� �� 0 �� 1, �����, ��� 360 * x[6] - ������������ ���� ����� ������� ChCh
* \param x[7] - ����� �� 0 �� 1 : pT1[0] + 1000 * (1-2*x[7]) - NS ���������� ����� ChCh
* \param x[8] - ����� �� 0 �� 1 : pT1[1] + 1000 * (1-2*x[8]) - EW ���������� ����� ChCh
* \return ������ TrajectoryTemplate* . ���������� ���������� ��� ������� (solve).
*/
std::vector<TrajectoryTemplate*> well2CHCH(const std::vector<double>& x,
                                           const Eigen::Vector3d& pinit,
                                           const Eigen::Vector3d& pT1,
                                           const Eigen::Vector3d& pT3);
/*
* \brief ������� ��� ����������� ����� ����������,
* ��������� �� �������� Hold-ChCh-ChCh.��������� ������������ � �������.
* \param pinit - 3� ������ ������ ���������� (����������)
* \param pT1 - 3� ������ ���� �1 ���������� (����������)
* \param pT3 - 3� ������ ���� �3 ���������� (����������)
*/
void OptimizeHorizontal(const Eigen::Vector3d& pinit,
                        const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3);

/*
* \brief ������� ��� ���������������� ����������� ������ ����������
* � ��������������, ��������� �� �������� Hold-ChCh-ChCh.
* ��������� ������������ � �������.
* \param pinit - ������ 3� �������� ������ ���������� (����������)
* \param pT1 - ������ 3� �������� ���� �1 ���������� (����������)
* \param pT3 - ������ 3� �������� ���� �3 ���������� (����������)
*/
void OptimizeHorizontals(const std::vector<Eigen::Vector3d>& pinits,
                         const std::vector<Eigen::Vector3d>& targets1,
                         const std::vector<Eigen::Vector3d>& targets3);
/*
* \brief ������� ��� ������������� ����������� ������ ����������,
* ��������� �� �������� Hold-ChCh-ChCh. ��������� ������������ � �������.
* \param pinit - ������ 3� �������� ������ ���������� (����������)
* \param pT1 - ������ 3� �������� ���� �1 ���������� (����������)
* \param pT3 - ������ 3� �������� ���� �3 ���������� (����������)
*/
void OptimizeTogether(const std::vector<Eigen::Vector3d>& pinits,
                      const std::vector<Eigen::Vector3d>& targets1,
                      const std::vector<Eigen::Vector3d>& targets3);

void testHorizontal(const Eigen::VectorXd& x,
                    std::vector<TrajectoryTemplate*>& tt);

void writeGG(const std::string& filename, const nlohmann::json& jresults);

void aloneOpt();