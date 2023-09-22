#pragma once
#include "TrajectoryTemplate.h"
#include "CurveHold.h"

/**
 * Класс CurveHoldCurveHold представляет собо шаблон построения, состоящий из
 * двух последовательных шаблонов CurveHold. Используется для построения
 * профилей горизонтальных скважин.
 * Построение данного шаблона выполняется итерационным методом.
 */
class CurveHoldCurveHold : public TrajectoryTemplate
{
private:
    Eigen::Vector3d p1; // координаты начальной точки траектории (координаты N,E,V)
    Eigen::Vector3d p4; // координаты конечной точки траектории (координаты N,E,V)
    Eigen::Vector3d t1; // единичный касательный вектор начальной точки
    Eigen::Vector3d t4;	// единичный касательный вектор конечной точки
    Eigen::Vector3d pInter;  // координаты цели Т1. если последний Hold нулевой
                             // длины, то совпадает с p4

    Eigen::Vector3d r1; // координаты центра окружности первого участка curve
    Eigen::Vector3d r4; // координаты центра окружности второго участка curve

    double tetta1;  // значение зенитного угла (в градусах) касательного вектора
                    // t1 начальной точке
    double phi1;  // значение азимутального угла (в градусах) касательного
                  // вектора t1 начальной точке
    double tetta4;  // значение зенитного угла (в градусах) касательного вектора
                    // t4 конечной точке
    double phi4;  // значение азимутального угла (в градусах) касательного
                  // вектора t4 конечной точке

    double R1; // радиус кривизны дуги к p1
    double R2; // радиус кривизны дуги к p4
    size_t nums;// число точек траектории


    double betta;  // длина финального участка стабилизации hold, если betta = 0.
                   // то шаблон CurveHoldCurve
    double eps; // критерий остановки для построения шаблона.

    Eigen::Vector3d t;  // единичный касательный вектор участка hold между двумя
                        // кривыми (определяется итерационно)
    Eigen::Vector3d p1Inter; // координаты конца первого участка дуги.
    Eigen::Vector3d p4Inter; // координаты начала второго участка дуги.
    double alpha1; // угол первого участка дуги в радианах.
    double alpha2; // угол второго участка дуги в радианах.
    double holdLength; // длина промежуточного участка стабилизации.
    int condition;
    /**
      * @brief Функция построения шаблона (то есть поиска недостающих параметров
      * по входным параметрам). Реализацию данного метода см. в статье Sawaryn et. al. 2005.
     */
    int fit();


public:
    /**
     * @brief Конструктор шаблона построения по известным значениям координат
     * целей Т1 и Т3 (для горизонтальных скважин).
     * @param p1 - координаты начальной точки;
     * @param tetta1 - зенитный угол(°) в начальной точке;
     * @param phi1  - азимутальный угол(°) в начальной точке;
     * @param R1 - радиус кривизны первого участка искривления (в метрах);
     * @param R2 - радиус кривизны второго участка искривления (в метрах);
     * @param pT1 - координаты цели Т1;
     * @param pT3 - координаты цели T3;
     * @param eps - критерий остановки для построения;
     * @param nums - число точек для интерполяции.
     */
    CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1, double phi1,
                    double R1, double R2, const Eigen::Vector3d& pT1,
                    const Eigen::Vector3d& pT3, double eps = 10e-4,
                    size_t nums = 100);
     /**
      * @brief Конструктор шаблона построения по известному положению и
      * касательному вектору конечной точки и длине последнего участка стабилизации.
      * @param p1 - координаты начальной точки;
      * @param tetta1 - зенитный угол(°) в начальной точке;
      * @param phi1  - азимутальный угол(°) в начальной точке;
      * @param R1 - радиус кривизны первого участка искривления (в метрах);
      * @param R2 - радиус кривизны второго участка искривления (в метрах);
      * @param p4 - координаты конечной точки;
      * @param tetta4 - зенитный угол(°) в конечной точке;
      * @param phi4  - азимутальный угол(°) в конечной точке;
      * @param betta - длина последнего участка стабилизации;
      * @param eps - критерий остановки для построения;
      * @param nums - число точек для интерполяции.
      */
    CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1, double phi1,
                       double R1, double R2, const Eigen::Vector3d& p4,
                       double tetta4, double phi4, double betta = 0.0,
                       double eps = 10e-4, size_t nums = 100);
    // Перегруженные методы базового класса. (описание см. TrajectoryTemplate)
    int getCondition() override;
    void points(CoordinateSystem coordinateSystem) override;
    double length() override;
    double getTortuosity() override;
    Eigen::Vector3d FunctionPoint(double md) override;
    Eigen::Vector3d FunctionTangent(double md) override;

    void getInitPoint(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
    void getTarget1Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
    void getTarget3Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
};

