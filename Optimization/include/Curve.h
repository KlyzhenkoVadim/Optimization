#pragma once
#include "TrajectoryTemplate.h"
// Тип построения дуги. Если известен конкретный радиус/dls, то DLS.
// Если необходимо построить дугу на конкретную глубину, то TVD (dls сам посчитается).
enum class TypeCurve { DLS,TVD };
/**
 * @brief Класс шаблона дуги кривой Curve. Может быть построен как по
 * фиксированному dls, так и на конкретный tvd().
 */
class Curve : public TrajectoryTemplate {
private:
    Eigen::Vector3d pi,pf; // декартовы координаты начали и конца дуг.
    double inc1, azi1; // зенитный и азимутальный углы(°) начала дуги.
    double inc2, azi2; // зенитный и азимутальный углы(°) конца дуги.
    double R; // радиус дуги.
    TypeCurve type; // тип построения
    double RTVD;  // величина зависящая от type. Если type = DLS, то радиус
                  // дуги, иначе TVD конца дуги.
    Eigen::Vector3d t1, t2; // касательные вектора начала и конца дуги.
    double alpha; // угол дуги в радианах.
    size_t nums; // параметр интерполяции - число точек.
    int condition; // состояние - построена(0) /не построена(-1).
    /**
     * @brief Функция построения шаблона (то есть поиска недостающих параметров
     * по входным параметрам).
    */
    int fit();

public:
    /**
    * @brief Конструктор.
    * @param pi - координаты начальной точки,
    * @param inc1,azi1 - зенитный и азимутальный углы(°) начала дуги.
    * @param inc2,azi2 - зенитный и азимутальный углы(°) конца дуги.
    * @param RTVD - в зависимости от типа построения радиус или tvd.
    * @param type - тип построения
    * @param nums - параметр интерполяции - число точек.
    */
    Curve(const Eigen::Vector3d& pi, double inc1, double azi1, double inc2,
          double azi2, double RTVD, TypeCurve type = TypeCurve::DLS,
          size_t nums = 20);
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
};
