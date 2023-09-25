#pragma once
#include <vector>
#include <functional>

namespace pso
{
    constexpr double EPSILON = 1e-10;
}

/* \brief Результат оптимизации PSO;
* \param optVec - оптимальный вектор размерности dimension;
* \param cost - минимальное значение целевой функции;
* \param costHist - значение минимальной целевой функции на каждой итерации,
вектор размерности numIterations+1.
*/
struct PsoValueType
{
    std::vector<double> optVec;
    double cost;
    std::vector<double> costHist;
};
/* \brief Алгоритм роя частиц.
* \param func - реализация метода оптимизируемой функции;
* \param minValues - ограничения на вектор параметров снизу;
* \param maxValues - ограничения сверху;
* \param numAgents - гиперпараметр числа частиц;
* \param numIterations - максимальное число итераций;
* \param inertia - гиперпараметр, коэффициент инерционности;
* \param socCoef - гиперпараметр, "социальный" коэффициент;
* \param indCoef - гиперпараметр, "индивидуальный" коэффициент;
* \param isHistoryNeeded - если true, то записывается costHist.
*/
PsoValueType PSO(std::function<double(const std::vector<double>&)> func,
    const std::vector<double>& minValues, const std::vector<double>& maxValues,
    size_t numAgents, size_t numIterations = 100,
    const std::vector<double>& inertia = std::vector<double>(500, 0.9),
    double socCoef = 0.3, double indCoef = 0.5, bool isHistoryNeeded = false);

