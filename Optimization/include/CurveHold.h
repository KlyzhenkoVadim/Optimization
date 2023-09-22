#pragma once
#include "TrajectoryTemplate.h"
/**
 * @brief Класс CurveHold представляет собой шаблон построения состоящий из дуги
 * и последующего участка стабилизации. (см. статью)
 */
class CurveHold : public TrajectoryTemplate
{
private:
    Eigen::Vector3d p1; // координата начальной точки траектории (координаты N,E,V)
    Eigen::Vector3d p3; // координата конечной точки траектории (координаты N,E,V)
    Eigen::Vector3d t1; // единичный касательный вектор начальной точки
    Eigen::Vector3d t2; // единичный касательный вектор конечной точки дуги

    double teta; // значения зенитного угла (в градусах) касательного вектора t1 начальной точке
    double phi; // значения азимутального угла (в градусах) касательного вектора t1 начальной точке
    double R; // радиус кривизны траектории (R>0)
    size_t nums;// число точек траектории

    double alpha; // угол дуги в радианах
    double betta; // длина участка стабилизации
    int condition;
    /**
     * @brief Функция построения шаблона (то есть поиска недостающих параметров
     * по входным параметрам).
    */
    int fit();

public:
    /**
    * @brief Конструктор шаблона построения по значениям зенитного и азимутального углов.
    * @param p1 - координаты начальной точки;
    * @param p3 - координаты конечной точки;
    * @param teta - зенитный угол(°) в начальной точке;
    * @param phi - азимутальный угол(°) в начальной точке;
    * @param R - радиус кривизны;
    * @param nums - параметр интерполяции, число точек.
    */
    CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3, double teta,
              double phi, double R, size_t nums = 50);
    /**
    * @brief Конструктор шаблона построения по касательному вектору.
    * @param p1 - координаты начальной точки;
    * @param p3 - координаты конечной точки;
    * @param t1 - касательный вектор в начальной точке
    * @param R - радиус кривизны;
    * @param nums - параметр интерполяции, число точек.
    */
    CurveHold(const Eigen::Vector3d& p1, const Eigen::Vector3d& p3,
              const Eigen::Vector3d& t1, double R, size_t nums = 50);
    // см. описание TrajectoryTemplate.
    int getCondition() override;
    void points(CoordinateSystem coordinateSystem) override;
    double length() override;
    double getTortuosity() override;
    void getInitPoint(CoordinateSystem coordinateSystem =
                          CoordinateSystem::CARTESIAN) override;
    void getTarget1Point(CoordinateSystem coordinateSystem =
                             CoordinateSystem::CARTESIAN) override;
    void getTarget3Point(CoordinateSystem coordinateSystem =
                             CoordinateSystem::CARTESIAN) override;
    Eigen::Vector3d FunctionPoint(double md) override;
    Eigen::Vector3d FunctionTangent(double md) override;
    /**
    * @brief Функция возвращает координаты конца дуги шаблона.
    */
    Eigen::Vector3d getPointInterpol();
    /**
    * @brief Функция возвращает касательный вектор конца дуги, а значит и конца шаблона.
    */
    Eigen::Vector3d getTangent2();

};

