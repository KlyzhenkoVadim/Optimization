#pragma once
#include "TrajectoryTemplate.h"
#include <cmath>
#include <iostream>
#include "Penalties.h"
/**
 * @brief Функция для вычисления фактора разделения пары траекторий.
 * (separation factor). По простейшей модели travelling cyllinder(!!!)
 * @param pCartesianW1 - набор координат (NS,EW,TVD) первой траектории (референс)
 * @param pMDW1 - набор md и координат касательного вектора (MD,NSt,EWt,TVDt)
 *  первой траектории (референс)
 * @param pCartesianW2 - набор координат (NS,EW,TVD) второй траектории (оффсет)
 * @param pMDW2 - набор md и координат касательного вектора (MD,NSt,EWt,TVDt)
 * второй траектории (оффсет)
 * @param TVDstart - значение глубины, начиная с которой вычисляется фактор разделения
 * (предназначено для случая общего вертикального участка.
 * @param actFunc - флаг, если true - вычисляется сигмоида от фактора разделения,
 * если false - вычисляется сам фактор разделения.
 * @param penalty - штраф, если фактор разделения меньше 1.5.
 * @return фактор разделения, или значение сигмоиды от него.
*/
double sepFactor(const std::vector<Eigen::Vector3d>& pCartesianW1,
                 const std::vector<Eigen::Vector4d>& pMDW1,
                 const std::vector<Eigen::Vector3d>& pCartesianW2,
                 const std::vector<Eigen::Vector4d>& pMDW2, double TVDstart = 0,
                 bool actFunc = true, double penalty = 20);
/**
 * @brief Вычисление индекса сложности бурения траектории (ddi)
 * @param well - вектор TrajectoryTemplate* описывающий траекторию. 
 * @param pCartesian - интерполированные точки на траектории
 * @param actFunc - если true, возвращается сигмоида от ddi. 
 * @param penalty - штраф за превышение ddi значения 6.25
 * @return в зависимости от actFunc (ddi или сигмоида от него)
*/
double DDI(std::vector<TrajectoryTemplate*>& well,
           const std::vector<Eigen::Vector3d>& pCartesian, bool actFunc = true,
           double penalty = 5);
/**
 * @brief Вычисление ERD (Extended Reach Drilling) для горизонтальных скважин.
 * @param points - вектор интерполированных точек траектории
 * @return значение метрики.
*/
double ERD(const std::vector<Eigen::Vector3d>& points);
/**
 * @brief Вычисление целевого функционала текущей траектории 
 * при наличии уже построенных!
 * @param mainWell - вектор TrajectoryTemplate* оптимизируемой траектории.
 * @param pCTrajectories - набор векторов точек уже построенных траекторий.
 * @param pMDTrajectories - набор касательных векторов и md построенных траекторий.
 * @param SepFactorShift - параметр TVDStart для функции SepFactor.
 * @param penalty - штраф за непостроение.
 * @return целевая функция: нормированная длина + сигмоида sepFactor + сигмоида DDI (ERD)
*/
double orderScore1(
    std::vector<TrajectoryTemplate*>& mainWell,
    const std::vector<std::vector<Eigen::Vector3d>>& pCTrajectories,
    const std::vector<std::vector<Eigen::Vector4d>>& pMDTrajectories,
    double SepFactorShift = 0, double penalty = 1000);
/**
* @TODO Удалить, не используется
*/
double scoreSolver(std::vector<TrajectoryTemplate*>& tmp,
                   const WellTrajectoryConstraints& cs, double penalty = 1000);
/**
 * @brief Функция для определения dls
 * @param tangent1 касательный вектор начала дуги
 * @param tangent2 касательный вектор конца дуги
 * @return dog leg severity 
*/
double dls(Eigen::Vector3d& tangent1, Eigen::Vector3d& tangent2);
/**
 * @brief Вычисление суммарного tortuosity траектории.
 * @param well вектор TrajectoryTemplate* траектории.
 * @return 
*/
double tortuosity(std::vector<TrajectoryTemplate*>& well);
