#include "Curve.h"

Curve::Curve(const Eigen::Vector3d& pi, double inc1, double azi1, double inc2,
             double azi2, double RTVD, TypeCurve type, size_t nums) {
    this->pi = pi;
    this->azi1 = azi1;
    this->inc1 = inc1;
    this->azi2 = azi2;
    this->inc2 = inc2;
    this->t1 = calcTangentVector(azi1, inc1);
    this->t2 = calcTangentVector(azi2, inc2);
    this->nums = nums;
    this->RTVD = RTVD;
    this->type = type;
    this->condition =
        Curve::fit();  // построение происходит сразу в конструкторе. то есть
                       // при создании объекта можно сразу сказать, построен он
                       // или нет
}

int Curve::fit() 
{
    // Получение alpha при помощи скалярного произведения (t1,t2)
    double dotprod = t1.dot(t2) / t1.norm() / t2.norm();
    dotprod = fabs(dotprod - 1) < EPSILON ? 1 : dotprod;
    alpha = acos(dotprod);
    if (type == TypeCurve::DLS) {
        this->R = RTVD;
        pf = pi + R * tan(alpha / 2) * (t1 + t2);
    }
    else if (type == TypeCurve::TVD) {
        // не построится, если дуга находится в горизонтальной плоскости.
        if (abs(t1[2] + t2[2]) < EPSILON) {
            return -1;
        }
        //  если alpha == 0 то дуга нулевая, и финальная точка = начальной. вне
        //  зависимости от TVD
        if (abs(alpha) < EPSILON) {
            this->R = 1 / EPSILON;
            pf = pi;
        }
        else {
            // выражаем радиус и пользуемся привычной формулой.
            this->R = (RTVD - pi[2]) / (tan(alpha / 2) * (t1[2] + t2[2]));
            pf = pi + R * tan(alpha / 2) * (t1 + t2);
        }
    }
    return 0;
}

int Curve::getCondition()
{
    return condition;
}

void Curve::points(CoordinateSystem coordsys) {
    if (coordsys == CoordinateSystem::CARTESIAN) {
        // если угол нулевой, то возвращаем одну точку.
        if (abs(alpha) < EPSILON)
            pointsCartesian.push_back(pi);
        else {
            pointsCartesian = calcInterpolCartesianPoints(pi, t1, t2, R, alpha, nums);
        }
    }
    else {
        if (abs(alpha) < EPSILON)
            pointsMD.push_back({ 0,t1[0],t1[1],t1[2] });
        else {
            pointsMD = calcInterpolMDPoints(pi, t1, t2, R, alpha, nums);
        }
    }
}

double Curve::length() {
    return alpha * R;
}

double Curve::getTortuosity()
{
    return alpha;
}

Eigen::Vector3d Curve::FunctionPoint(double md) // md \in [0,1]
{
    Eigen::Vector3d t_star = FunctionTangent(md);
    return pi + R * tan(md*alpha/2) * (t1 + t_star);
}

Eigen::Vector3d Curve::FunctionTangent(double md)
{
    // формула взята из основной статьи (см.TrajectoryTemplate)
    if (alpha < EPSILON)
        return t1;
    return t1 * sin((1 - md) * alpha) / sin(alpha) + t2 * sin(md*alpha) / sin(alpha);
}

void Curve::getInitPoint(CoordinateSystem coordinateSystem ) {
    if (coordinateSystem == CoordinateSystem::CARTESIAN)
        pointInitial = this->pi;
    else
        pointInitialMD = { 0,t1[0],t1[1],t1[2] };
}

void Curve::getTarget1Point(CoordinateSystem coordinateSystem) {
    if (condition == 0)
    {
        if (coordinateSystem == CoordinateSystem::CARTESIAN)
            pointT1 = pf;
        else
            pointMDT1 = { length(),t2[0],t2[1],t2[2] };
    }
    else
    {
        Curve::getInitPoint(coordinateSystem);
        if (coordinateSystem == CoordinateSystem::CARTESIAN)
            pointT1 = pointInitial;
        else
            pointMDT1 = pointInitialMD;
    }
}

void Curve::getTarget3Point(CoordinateSystem coordinateSystem) {
    if (condition == 0)
    {
        if (coordinateSystem == CoordinateSystem::CARTESIAN)
            pointT3 = pf;
        else
            pointMDT3 = { length(),t2[0],t2[1],t2[2] };
    }
    else
    {
        Curve::getInitPoint(coordinateSystem);
        if (coordinateSystem == CoordinateSystem::CARTESIAN)
            pointT3 = pointInitial;
        else
            pointMDT3 = pointInitialMD;
    }
}