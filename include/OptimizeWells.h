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
* \brief ‘ункци€ дл€ создани€ траектории по заданному набору параметров
* \param x вектор параметров
* \param pinit - 3д вектор начала траектории
* \param pT1 - 3д вектор цели “1 конца траектории
* \param pT3 - 3д вектор цели “3 конца траектории
* \param x[0] - число от 0 до 1, такое, что 900 * x[0] - длина вертикального участка
* \param x[1] - dls первого участка искривлени€
* \param x[2] - dls второго участка искривлени€
* \param x[3] - dls третьего участка искривлени€
* \param x[4] - dls четвЄртого участка искривлени€
* \param x[5] - число от 0 до 1, такое, что 180 * x[5] - зенитный угол конца первого ChCh
* \param x[6] - число от 0 до 1, такое, что 360 * x[6] - азимутальный угол конца первого ChCh
* \param x[7] - число от 0 до 1 : pT1[0] + 1000 * (1-2*x[7]) - NS координата конца ChCh
* \param x[8] - число от 0 до 1 : pT1[1] + 1000 * (1-2*x[8]) - EW координата конца ChCh
* \return вектор TrajectoryTemplate* . ѕостроение происходит вне функции (solve).
*/
std::vector<TrajectoryTemplate*> well2CHCH(const std::vector<double>& x,
                                           const Eigen::Vector3d& pinit,
                                           const Eigen::Vector3d& pT1,
                                           const Eigen::Vector3d& pT3);
/*
* \brief ‘ункци€ дл€ оптимизации одной траектории,
* состо€щей из шаблонов Hold-ChCh-ChCh.–езультат записываетс€ в консоль.
* \param pinit - 3д вектор начала траектории (фиксирован)
* \param pT1 - 3д вектор цели “1 траектории (фиксирован)
* \param pT3 - 3д вектор цели “3 траектории (фиксирован)
*/
void OptimizeHorizontal(const Eigen::Vector3d& pinit,
                        const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3);

/*
* \brief ‘ункци€ дл€ последовательной оптимизации набора траекторий
* с регул€ризацией, состо€щих из шаблонов Hold-ChCh-ChCh.
* –езультат записываетс€ в консоль.
* \param pinit - вектор 3д векторов начала траектории (фиксирован)
* \param pT1 - вектор 3д векторов цели “1 траектории (фиксирован)
* \param pT3 - вектор 3д векторов цели “3 траектории (фиксирован)
*/
void OptimizeHorizontals(const std::vector<Eigen::Vector3d>& pinits,
                         const std::vector<Eigen::Vector3d>& targets1,
                         const std::vector<Eigen::Vector3d>& targets3);
/*
* \brief ‘ункци€ дл€ одновременной оптимизации набора траекторий,
* состо€щих из шаблонов Hold-ChCh-ChCh. –езультат записываетс€ в консоль.
* \param pinit - вектор 3д векторов начала траектории (фиксирован)
* \param pT1 - вектор 3д векторов цели “1 траектории (фиксирован)
* \param pT3 - вектор 3д векторов цели “3 траектории (фиксирован)
*/
void OptimizeTogether(const std::vector<Eigen::Vector3d>& pinits,
                      const std::vector<Eigen::Vector3d>& targets1,
                      const std::vector<Eigen::Vector3d>& targets3);

void testHorizontal(const Eigen::VectorXd& x,
                    std::vector<TrajectoryTemplate*>& tt);

void writeGG(const std::string& filename, const nlohmann::json& jresults);

void aloneOpt();