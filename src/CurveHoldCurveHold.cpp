#include <algorithm>
#include "CurveHoldCurveHold.h"

CurveHoldCurveHold::CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1,
                                       double phi1, double R1, double R2,
                                       const Eigen::Vector3d& p4, double tetta4,
                                       double phi4, double betta, double eps,
                                       size_t nums) {

    this->p1 = p1;
    this->p4 = p4;
    this->t1 = calcTangentVector(phi1, tetta1);
    this->phi1 = phi1;
    this->phi4 = phi4;
    this->tetta1 = tetta1;
    this->tetta4 = tetta4;
    this->t4 = calcTangentVector(phi4, tetta4);
    this->R1 = R1;
    this->R2 = R2;
    this->eps = eps;
    this->nums = nums;
    this->betta = betta;
    this->pInter = p4 - t4 * betta;
    this->condition = CurveHoldCurveHold::fit();
};

CurveHoldCurveHold::CurveHoldCurveHold(const Eigen::Vector3d& p1, double tetta1,
                                       double phi1, double R1, double R2,
                                       const Eigen::Vector3d& pT1,
                                       const Eigen::Vector3d& pT3, double eps,
                                       size_t nums) {

    this->p1 = p1;
    this->p4 = pT3;
    this->pInter = pT1;
    this->phi1 = phi1;
    this->tetta1 = tetta1;
    this->tetta4 = 0.;
    this->phi4 = 0.;
    this->t1 = calcTangentVector(phi1, tetta1);
    Eigen::Vector3d tmpVec = pT3 - pT1;
    this->betta = tmpVec.norm();
    tmpVec.normalize();
    this->t4 = tmpVec;
    this->R1 = R1;
    this->R2 = R2;
    this->nums = nums;
    this->eps = eps;
    this->condition = CurveHoldCurveHold::fit();
};

int CurveHoldCurveHold::fit() 
{
    //Проверка на пересечение окружностей
    Eigen::Vector3d b1 = (pInter - p1).cross(t1);
    Eigen::Vector3d n1 = t1.cross(b1);
    n1.normalize();

    Eigen::Vector3d b4 = (p1 - pInter).cross(-t4);
    Eigen::Vector3d n4 = -t4.cross(b4);
    n4.normalize();
    r1 = p1 + R1 * n1;
    r4 = pInter + R2 * n4;

    if (b1.norm() < EPSILON && b4.norm() < EPSILON) {
        t = t1;
        holdLength = (pInter - p1).norm();
        alpha1 = 0;
        alpha2 = 0;
        p1Inter = p1;
        p4Inter = p4;
    }
    
    else
    {

        if ((pInter - r1).norm() < R1 || (p1 - r4).norm() < R2) {
            return -1;
        }

        /**
         * @brief Данный метод в зависимости от флага true/false
         * строит шаблон curveHold
         * @param p - координаты конечной точки;
         * @param flag = false от точки p1 c направлением t1 до p с радиусом R1;
         * @param flag = true от точки p4 c направлением (-t4) до p с радиусом R2;
         * @return Если построение возможно, то возвращает угол участка кривизны
         * этого CurveHold, иначе -1.
         */
        auto getaAlphaCurrCurvate = [&](bool flag, const auto& p) { 
            // flag = 0 - first curveHold
            if (false == flag) {
                CurveHold firstCurveHold(p1, p, t1, R1);
                if (firstCurveHold.getCondition() == 0)
                    return firstCurveHold.getTortuosity();
                else
                    return -1.;
            }
            else {
                CurveHold secondCurveHold(pInter, p, -t4, R2);
                if (secondCurveHold.getCondition() == 0)
                    return secondCurveHold.getTortuosity();
                else
                    return -1.;
            }
        };

        std::vector<double> alphaFirst(2);
        std::vector<double> alphaSecond(2);
        alphaFirst[0] = getaAlphaCurrCurvate(0, pInter);
        if (alphaFirst[0] < -EPSILON)
            return -1;
        Eigen::Vector3d p1j = p1 + R1 * tan(alphaFirst[0] / 2) * t1;
        alphaSecond[0] = getaAlphaCurrCurvate(1, p1j);
        if (alphaSecond[0] < -EPSILON)
            return -1;
        Eigen::Vector3d pInterj = pInter - R2 * tan(alphaSecond[0] / 2) * t4;

        const size_t maxIter = 1000;
        size_t iter = 0;
        // схему расчёта см. Sawaryn et al 2005
        while (iter < maxIter) {
            alphaFirst[1] = getaAlphaCurrCurvate(0, pInterj);
            // Если хотя бы один расчёт CurveHold'a не построился - не построен
            // весь шаблон
            if (alphaFirst[1] < -EPSILON)
                return -1;
            p1j = p1 + R1 * tan(alphaFirst[1] / 2) * t1;
            alphaSecond[1] = getaAlphaCurrCurvate(1, p1j);
            if (alphaSecond[1] < -EPSILON)
                return -1;
            pInterj = pInter - R2 * tan(alphaSecond[1] / 2) * t4;

            if (sqrt((alphaFirst[1] - alphaFirst[0]) *
                         (alphaFirst[1] - alphaFirst[0]) +
                     (alphaSecond[1] - alphaSecond[0]) *
                         (alphaSecond[1] - alphaSecond[0])) < eps) {
                break;
            }

            std::reverse(alphaFirst.begin(), alphaFirst.end());
            std::reverse(alphaSecond.begin(), alphaSecond.end());
            iter += 1;  // !!!
        }
        // Если не было выхода из цикла, шаблон не построен.
        if (iter == 999) {
            return -1;
        }
        // в результате схемы есть две точки на касательных
        // из них мы получим направление промежуточного холда.
        t = (pInterj - p1j);
        t.normalize();
        alpha1 = alphaFirst[1];
        alpha2 = alphaSecond[1];
        // Отсюда получим неизвестные точки участков искривления
        // и длину холда.
        p1Inter = p1 + R1 * tan(alpha1 / 2) * (t1 + t);
        p4Inter = pInter - R2 * tan(alpha2 / 2) * (t4 + t);
        holdLength = (p4Inter - p1Inter).norm();
        // здесь проверка того, 
        // что касательные вектора промежуточных точек
        //  на двух участках искривления сонаправлены! 
        // (иначе шаблон будет с изломами, похож на гиперболу)
        Eigen::Vector3d tmp = (p4Inter - p1Inter) / holdLength;
        // А если холд = 0 ???
        for (size_t idx = 0; idx < t.size(); ++idx) {
            if (fabs(tmp[idx] - t[idx]) > 1e-4) {
                return -1;
            }
        }
    }
    return 0;
}

int CurveHoldCurveHold::getCondition()
{
    return condition;
}

void CurveHoldCurveHold::getInitPoint(CoordinateSystem coordinateSystem) {
    if (coordinateSystem == CoordinateSystem::CARTESIAN)
        pointInitial = this->p1;
    else
        pointInitialMD = { 0,t1[0],t1[1],t1[2] };
}

void CurveHoldCurveHold::getTarget3Point(CoordinateSystem coordinateSystem) {
    if (coordinateSystem == CoordinateSystem::CARTESIAN)
        pointT3 = this->p4;
    else
        pointMDT3 = { length(),t4[0],t4[1],t4[2] };
}

void CurveHoldCurveHold::getTarget1Point(CoordinateSystem coordinateSystem) {
    if (coordinateSystem == CoordinateSystem::CARTESIAN)
        pointT1 = this->pInter;
    else
        pointMDT1 = { length() - betta,t4[0],t4[1],t4[2] };
}

Eigen::Vector3d CurveHoldCurveHold::FunctionPoint(double md) // md [0,1]
{
    // Необходимо определить, какому сегменту (их всего может быть 4)
    // принадлежит точка
    // md*L [0,arc1]
    // md*L (arc1,arc1+HoldLength]
    // md*L (arc1+HoldLength,arc1+HoldLength+arc2]
    // md*L (length-betta,length]
    double L = length(); // полная длина
    double arc1 = alpha1 * R1; // длина первого curve'a
    double arc2 = alpha2 * R2; // длина второго curve'a
    if (md * L < arc1)
    {
        return p1 + R1 * tan(md * L / arc1 * alpha1 / 2) *
                        (t1 + FunctionTangent(md));
    }
    if (md * L - arc1 > 0 && md * L - arc1 < holdLength)
    {
        return p1Inter + (md * L - arc1) * t;
    }
    if (md * L - arc1 - holdLength > 0 && md * L - arc1 - holdLength < arc2)
    {
        double s = (md * L - arc1 - holdLength) / arc2;
        return p4Inter + R2 * tan(s * alpha2 / 2) * (t + FunctionTangent(md));
    }
    else
        return p4 - (1 - md) * L * t4;
}

Eigen::Vector3d CurveHoldCurveHold::FunctionTangent(double md) // md [0,1]
{
    double L = length();
    double arc1 = alpha1 * R1, arc2 = alpha2 * R2;
    if (md * L < arc1)
    {
        return t1 * sin((1 - md * L / arc1) * alpha1) / sin(alpha1) +
               t * sin(md * L / arc1 * alpha1) / sin(alpha1);
    }
    if (md * L - arc1 > 0 && md * L - arc1 < holdLength)
    {
        return t;
    }
    if (md * L - arc1 - holdLength > 0 && md * L - arc1 - holdLength < arc2)
    {
        double s = (md * L - arc1 - holdLength) / arc2;
        return t * sin((1 - s) * alpha2) / sin(alpha2) +
               t4 * sin(s * alpha2) / sin(alpha2);
    }
    else
        return t4;
}


void CurveHoldCurveHold::points(CoordinateSystem coordinateSystem) {
    // 
    double arc1 = R1 * alpha1; // длина первого curve'a
    double arc2 = R2 * alpha2; // длина второго curve'a
    double h = length() / nums; // шаг интерполяции
    // На каждый сегмент число точек пропорционально длине,
    // однако на дугу не меньше 5 точек, а на отрезок не меньше 2
    int arc1Nums = std::max(5, int(arc1 / h)); // число точек на первом curve
    int arc2Nums = std::max(5, int(arc2 / h)); // число точек на втором curve
    int nHold1 = std::max(2, int(holdLength / h)); // число точек на промежуточном hold'e
    int nHold2 = std::max(2, int(betta / h)); // число точек на последнем hold'e

    if (coordinateSystem == CoordinateSystem::CARTESIAN) {
        std::vector<Eigen::Vector3d> pointsArc1;
        // интерполяция первой дуги, если alpha1 нулевой, то начальная точка
        if (abs(alpha1) < EPSILON) {
            pointsArc1.push_back(p1);
        }
        else {
            pointsArc1 =
                calcInterpolCartesianPoints(p1, t1, t, R1, alpha1, arc1Nums);
        }
        // интерполяция промежуточного холда.
        std::vector<Eigen::Vector3d> pointsHold1(nHold1 + 1);
        for (size_t idx = 0; idx < nHold1 + 1; ++idx) {
            pointsHold1[idx] = p1Inter + idx * holdLength / nHold1 * t;
        }
        // интерполяция второй дуги, если alpha2 нулевой, то точка конца дуги.
        std::vector<Eigen::Vector3d> pointsArc2;
        if (abs(alpha2) < EPSILON) {
            pointsArc2.push_back(pInter);
        }
        else {
            pointsArc2 = calcInterpolCartesianPoints(p4Inter, t, t4, R2, alpha2,
                                                     arc2Nums);
        }
        // интерполяция последнего холда.
        std::vector<Eigen::Vector3d> pointsHold2(nHold2 + 1);
        for (size_t idx = 0; idx < nHold2 + 1; ++idx) {
            pointsHold2[idx] = pInter + idx * betta / nHold2 * t4;
        }
        // здесь мы объединяем вектора точек каждого из сегментов в один вектор.
        pointsCartesian = pointsArc1;
            if (holdLength > EPSILON) {
            std::copy(pointsHold1.begin(), pointsHold1.end(),
                      std::back_inserter(pointsCartesian));
            }
            std::copy(pointsArc2.begin(), pointsArc2.end(),
                      std::back_inserter(pointsCartesian));
            if (betta > EPSILON) {
            std::copy(pointsHold2.begin(), pointsHold2.end(),
                      std::back_inserter(pointsCartesian));
            }
    }
    else {
        // для интерполяции касательных векторов аналогично.
        // интерполяция первой дуги.
        std::vector<Eigen::Vector4d> pointsArc1;
        if (abs(alpha1) < EPSILON) {
            pointsArc1.push_back({ 0,t1[0],t1[1],t1[2] });
        }
        else {
            pointsArc1 = calcInterpolMDPoints(p1, t1, t, R1, alpha1,arc1Nums);
        }
        // интерполяция промежуточного холда.
        std::vector<Eigen::Vector4d> pointsHold1(nHold1 + 1);
        for (size_t idx = 0; idx < nHold1 + 1; ++idx) {
            pointsHold1[idx] = { arc1 + idx * holdLength / nHold1, t[0], t[1], t[2] };
        }
        // интерполяция второй дуги.
        std::vector<Eigen::Vector4d> pointsArc2;
        if (abs(alpha2) < EPSILON) {
            pointsArc2.push_back({ 0,t4[0],t4[1],t4[2] });
        }
        else {
            // md интерполируются с нуля, поэтому нужно сдвинуть на arc1.
            pointsArc2 = calcInterpolMDPoints(p4Inter, t, t4, R2, alpha2, arc2Nums);
        }
        // сдвигаем интерполированные md на константу arc1.
        for (size_t idx = 0; idx < pointsArc2.size(); ++idx) {
            pointsArc2[idx][0] += arc1 + holdLength;
        }
        // интерполяция последнего холда.
        std::vector<Eigen::Vector4d> pointsHold2(nHold2 + 1);
        for (size_t idx = 0; idx < nHold2 + 1; ++idx) {
            pointsHold2[idx] = {arc1 + arc2 + holdLength + idx * betta / nHold2,
                                t4[0], t4[1], t4[2]};
        }
        pointsMD = pointsArc1;
        // объединяем всё в один вектор.
        if (holdLength > EPSILON) {
            std::copy(pointsHold1.begin(), pointsHold1.end(),
                      std::back_inserter(pointsMD));
        }
        std::copy(pointsArc2.begin(), pointsArc2.end(),
                  std::back_inserter(pointsMD));
        if (betta > EPSILON) {
            std::copy(pointsHold2.begin(), pointsHold2.end(),
                      std::back_inserter(pointsMD));
        }
    }
}

double CurveHoldCurveHold::length() {
    double arc1 = alpha1 < EPSILON ? 0 : R1 * alpha1;
    double arc2 = alpha2 < EPSILON ? 0 : R2 * alpha2;
    return arc1 + arc2 + holdLength + betta;
}

double CurveHoldCurveHold::getTortuosity()
{
    return alpha1 + alpha2;
}