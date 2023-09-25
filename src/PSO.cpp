#include "PSO.h"
#include <random>

PsoValueType PSO(std::function<double(const std::vector<double>&)> func,
                 const std::vector<double>& minValues,
                 const std::vector<double>& maxValues, size_t numAgents,
                 size_t numIterations, const std::vector<double>& inertia,
                 double socCoef, double indCoef, bool history) 
{
    size_t dimension{ minValues.size() };
    // Проверка того, что количество ограничений сверху и снизу совпадает.
    if (dimension != maxValues.size())
    {
        return PsoValueType();
    }
    // Создадим вектор максимальных скоростей.
    std::vector<double> vMax(dimension);
    for (size_t idx = 0; idx < dimension; ++idx) 
    {
        vMax[idx] = (fabs(maxValues[idx] - minValues[idx]) / 2);
    }
    // Вектор положений частиц в пространстве параметров
    std::vector<std::vector<double>> positions(numAgents); 
    // Вектор скоростей частиц в пространстве параметров
    std::vector<std::vector<double>> velocities(numAgents);
    // Глобальный минимум среди всех частиц
    std::vector<double> gBestPos(dimension);
    // "Персональный" минимум для каждой частицы
    // минимум среди всех положений отдельно взятой частицы.
    std::vector<std::vector<double>> pBestPos(numAgents);
    // создание случайного распределения
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0,1);
    double minFunc{ 1e3 };
    size_t idxMin;
    // Задание начального положения частицам (случайно).
    // В начальный момент скорости нулевые.
    for (size_t i = 0; i < numAgents; ++i) 
    {
        for (size_t j = 0; j < dimension; ++j) 
        {
            positions[i].push_back(minValues[j] +
                                   (maxValues[j] - minValues[j]) * dist(gen));
            velocities[i].push_back(0.);
        }
    }
    // Инициализация pBestPos и gBestPos 
    pBestPos = positions;
    idxMin = 0;
    minFunc = func(positions[0]);
    for (size_t k = 1; k < numAgents; ++k) 
    {
        double tmp = func(positions[k]);
        if (tmp - minFunc < pso::EPSILON) 
        {
            minFunc = tmp;
            idxMin = k;
        }
    }
    gBestPos = positions[idxMin];
    PsoValueType result;
    // Если требуется история изменения целевой в процессе оптимиации
    if (history)
    {
        result.costHist.reserve(numIterations + 1);
        result.costHist.push_back(minFunc);
    }
    double r1, r2;
    double R{ 0.5 };
    // Шаг итерации стандратный, за исключением точек, которые находятся в
    // текущем минимуме.
    // см. DOI:10.3233/FI-2010-370
    for (size_t iteration = 0; iteration < numIterations; ++iteration)
    {
        // r1 r2 \in \mathbb{U}[0;1]
        for (size_t idAg = 0; idAg < numAgents; ++idAg)
        {
            r1 = dist(gen);
            r2 = dist(gen);
            for (size_t idim = 0; idim < dimension; ++idim) 
            {
                auto& v = velocities[idAg][idim];
                auto& x = positions[idAg][idim];
                // см. DOI:10.3233/FI-2010-370
                if (idAg == idxMin)
                {
                    v = inertia[iteration] * v + R * (1. - 2. * dist(gen));
                }
                else 
                {
                    v = (inertia[iteration] * v +
                         r1 * indCoef * (pBestPos[idAg][idim] - x) +
                        r2 * socCoef * (gBestPos[idim] - x));
                }
                // проверка того, что скорости не вышла за границы
                if (v - vMax[idim] > pso::EPSILON)
                {
                    v = vMax[idim];
                }
                else if (v + vMax[idim] < -pso::EPSILON)
                {
                    v = -vMax[idim];
                }
                // Если вышли за границу, то шаг не делается.
                if ((x + v - maxValues[idim] < pso::EPSILON) &&
                    (x + v - minValues[idim] > -pso::EPSILON))
                {
                    x += v;
                }
                else if ((x - v - maxValues[idim] < pso::EPSILON) &&
                    (x - v - minValues[idim] > -pso::EPSILON)) // TODO: check . . .
                    x -= v;
            }
        }
        // Перезапись "персональных" минимальных положений
        // и глобального.
        double tmpf{ minFunc };
        for (size_t idxAg = 0; idxAg < numAgents; ++idxAg) 
        {
            double tmpfuncMean = func(positions[idxAg]);
            if (tmpfuncMean - func(pBestPos[idxAg]) < pso::EPSILON)
                pBestPos[idxAg] = positions[idxAg];
            if (tmpfuncMean - minFunc < pso::EPSILON) 
            {
                idxMin = idxAg;
                gBestPos = positions[idxAg];
                minFunc = func(gBestPos);
            }
        }
        if (history)
        {
            result.costHist.push_back(minFunc);
        }		
    }
    result.cost = minFunc;
    result.optVec = gBestPos;
    
    return result;
}