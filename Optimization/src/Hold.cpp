#include "Hold.h"

Hold::Hold(const Eigen::Vector3d& pi,const Eigen::Vector3d& pf, size_t nums) {
    this->pi = pi;
    this->pf = pf;
    this->nums = nums;
    this->direction = pf - pi;
    this->Length = direction.norm();
    direction.normalize();
    this->type = typeHold::pEnd;
    this->condition = Hold::fit();
}

Hold::Hold(const Eigen::Vector3d& pi, double inc, double azi, double L,
           typeHold type, size_t nums) {
    this->pi = pi;
    this->nums = nums;
    this->direction = calcTangentVector(azi, inc);
    this->type = type;
    this->MDTVD = L;
    this->condition = Hold::fit();
}

int Hold::fit()
{
    if (type == typeHold::pEnd)
    {
        return 0;
    }
    else
    {
        if (type == typeHold::md) {
            this->Length = MDTVD;
        }
        else if (type == typeHold::TVD) {
            if (abs(direction[2]) < EPSILON) {
                // если касательный вектор лежит в горизонтальной
                // плоскости, то по TVD не построится
                return -1;
            }
            this->Length = abs((MDTVD - pi[2]) / direction[2]); // md = (TVD2-TVD1)/cos(theta)
        }
        this->pf = pi + Length * direction;
        return 0;
    }
    
}

int Hold::getCondition()
{
    return condition;
}

double Hold::length() {
    return Length;
}

double Hold::getTortuosity()
{ // На прямом участке набор угла нулевой.
    return 0.;
}

Eigen::Vector3d Hold::FunctionPoint(double md) // md[0,1]
{
    return pi + md * length() * direction; // уравнение прямой в векторном виде.
}

Eigen::Vector3d Hold::FunctionTangent(double md) // md[0,1]
{
    return direction;
}

void Hold::getInitPoint(CoordinateSystem coordinateSystem) {
    if (coordinateSystem == CoordinateSystem::CARTESIAN)
        pointInitial = this->pi;
    else
        pointInitialMD = { 0,direction[0],direction[1],direction[2] };
}

void Hold::getTarget1Point(CoordinateSystem coordinateSystem) {
    if (condition == 0)
    {
        if (coordinateSystem == CoordinateSystem::CARTESIAN)
            pointT1 = this->pf;
        else
            pointMDT1 = { length(),direction[0],direction[1],direction[2] };
    }
    else
    {
        Hold::getInitPoint(coordinateSystem);
        if (coordinateSystem == CoordinateSystem::CARTESIAN)
        {
            pointT1 = pointInitial;
        }
        else
        {
            pointMDT1 = pointInitialMD;
        }
    }
}

void Hold::getTarget3Point(CoordinateSystem coordinateSystem) {
    if (condition == 0)
    {
        if (coordinateSystem == CoordinateSystem::CARTESIAN)
            pointT3 = this->pf;
        else
            pointMDT3 = { length(),direction[0],direction[1],direction[2] };
    }
    else
    {
        Hold::getInitPoint(coordinateSystem);
        if (coordinateSystem == CoordinateSystem::CARTESIAN)
        {
            pointT3 = pointInitial;
        }
        else
        {
            pointMDT3 = pointInitialMD;
        }
    }
}

void Hold::points(CoordinateSystem coordinateSystem) {
    if (coordinateSystem == CoordinateSystem::CARTESIAN) {
        if (Length < EPSILON)
            pointsCartesian.push_back(pi);
        else {
            for (size_t i = 0; i < nums; ++i) {
                pointsCartesian.push_back(pi + Length * direction * i / (nums - 1));
            }
        }
    }

    else {
        if (direction.norm() > EPSILON) {
            for (size_t id = 0; id < nums; ++id) {
                pointsMD.push_back({(Length * id / (nums - 1)), direction[0],
                                    direction[1], direction[2]});
            }
        }
    }
}