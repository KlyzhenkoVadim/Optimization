#pragma once
#include "TrajectoryTemplate.h"
/**
* Выбор типа построения шаблона Hold
* @param md - известна длина отрезка;
* @param TVD - известна длина отрезка по вертикали;
* @param pEnd - известны координаты конечной точки.
*/
enum class typeHold {md,TVD,pEnd};
/**
 * @brief Класс построения шаблона стабилизации Hold (отрезок в 3D). Может быть
 * построен при помощи различной комбинации входных данных (см.typeHold)
 */
class Hold : public TrajectoryTemplate
{
private:
    Eigen::Vector3d pi; // Декартова координата начала Hold'a
    Eigen::Vector3d pf; // Декартова координата конца Hold'a
    Eigen::Vector3d direction; // направление Hold'a - единичный вектор 
    double Length;  // Длина отрезка (если известна изначально - определяется в
                    // конструкторе, если нет - вычисляется
    double MDTVD;  // Величина зависящая от значения type. Если type = md - то
                   // длина отрезка, если type = TVD, то TVD. Если type = pEnd,
                   // то не заполняется. (???)
    typeHold type = typeHold::pEnd; // тип построения
    size_t nums; // параметр интерполяции - число точек.
    int condition; // состояние - построена (0) / не построена (-1)
    /**
     * @brief Функция построения шаблона (то есть поиска недостающих параметров
     * по входным параметрам).
     */
    int fit();

public:
    /**
    * @brief Конструктор для случая type = pEnd (известны координаты конечной точки).
    * @param pi - координаты NS,EW,TVD начальной точки;
    * @param pf - координаты NS,EW,TVD конечной точки;
    * @param nums - число точек для интерполяции.
    */
    Hold(const Eigen::Vector3d& pi, const Eigen::Vector3d& pf, size_t nums = 5);
    /**
    * @brief Конструктор для случая type = md или type = TVD.
    * @param pi - координаты NS,EW,TVD начальной точки;
    * @param inc - зенитный угол (°);
    * @param azi - азимутальный угол (°);
    * @param L - значение md или TVD в зависимости от typeHold;
    * @param type - тип построения (должен принимать либо md, либо TVD ) (???);
    * @param nums - число точек для интерполяции.
    */
    Hold(const Eigen::Vector3d& pi, double inc, double azi, double L,
        typeHold type = typeHold::md, size_t nums = 5);
    // Перегруженные методы базового класса. (описание см. TrajectoryTemplate)
    int getCondition() override;
    void points(CoordinateSystem coordinateSystem) override;
    double length() override;
    double getTortuosity() override;
    void getInitPoint(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
    void getTarget1Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
    void getTarget3Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) override;
    Eigen::Vector3d FunctionPoint(double md) override;
    Eigen::Vector3d FunctionTangent(double md) override;
};
