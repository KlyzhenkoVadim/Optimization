#pragma once

#include <PSO.h>
#include <TrajectoryTemplate.h>
#include <functional>
#include <optional>
#include <string>

namespace well_trajectory {
/**
 * @brief Структура для обозначения положения кустовой площадки (глубина
 * считается нулевой).
 */
struct Point2d {
    double north;
    double east;
};

struct Layer {
    double TVD, theta, phi;
};

struct Target {
    double northT1;
    double eastT1;
    double tvdT1;

    double northT3;
    double eastT3;
    double tvdT3;

    double H;

    double frGas;
    double frLiq;

    std::string name;
    std::string date;
};

/**
 * @brief Структура предназначена для того, чтобы задать ограничения следующего
 * характера: На глубине TVD = 900 траектория должна входить в соответствующий
 * пласт под углами в диапазоне [thetaMin,thetaMax] - для зенитного угла;
 * [phiMin,phiMax] - для азимута.
 * @todo Включить в солвер
 */
struct Constraint {
    Layer lMin, lMax;
};

/**
 * Структура для задания ограничений бурения.
 * @param minDepthFirstHold - минимальная глубина первого вертикального участка;
 * @param maxDLS - максимальный DLS(Dogleg severity, скорость набора угла),
 * величиная, обратно пропорциональная радиусу;
 * @param maxMD - максимальная длина по стволу скважины (measured depth);
 * @param maxDistEastWest, maxDistNorthSouth - максимальный значения отхода
 * траектории по горизонтали (отдельно для восток-запад и север-юг).
 */
struct WellTrajectoryConstraints {
    double minDepthFirstHold{std::numeric_limits<double>::quiet_NaN()};
    double maxDLS{std::numeric_limits<double>::quiet_NaN()};
    double maxMD{std::numeric_limits<double>::quiet_NaN()};
    double maxDistEastWest{std::numeric_limits<double>::quiet_NaN()};
    double maxDistNorthSouth{std::numeric_limits<double>::quiet_NaN()};
};

/**
 * @brief Функция, принимающая вектор параметров x (например, dls,hold)
 * . .)
 * @return вектор шаблонов построения.
 */
using wellType = std::function<std::vector<TrajectoryTemplate*>(
    const std::vector<double>& x)>;

/**
 * Класс для оптимизации траекторий по заданным ограничениям и координатам
 * кустовой площадки и целей. Пример использования:
 * WellTrajectorySolver solver;
 * solver.loadInputData(file_path);
 * solver.optimize();
 * solver.getTrajectoryLength();
 */
class WellTrajectorySolver {
   private:
    std::string name_{"well"};
    std::string date_{"01.01.2000"};
    Eigen::Vector3d pointInitial_;  // Координата устья траектории
    Eigen::Vector3d pointsT1_;  // Координаты цели траектории T1.
    // Координаты цели траектории T3 (может не быть).
    std::optional<Eigen::Vector3d> pointsT3_;
    // Траектория - вектор шаблонов построения.
    std::vector<TrajectoryTemplate*> trajectory_;
    // Характеризует количество НЕпостроенных шаблонов в траектории
    int condition_{0};
    double length_{0.0};
    // Вектор декартовых координат оптимальной траектории NS,EW,TVD
    std::vector<Eigen::Vector3d> pCtrajectory_;
    // Вектор инклинометрии оптимальной траектории MD,INC[°],AZI[°]
    std::vector<Eigen::Vector3d> pInclinometry_;
    // функция от вектора параметров. Для ННС и ГС
    // имеют разный вид.
    wellType mainWell_;
    WellTrajectoryConstraints constraints_{400., 1.5, 1 / EPSILON, 1 / EPSILON,
                                           1 / EPSILON};
    // Структура для записи результатов оптимизации
    PsoValueType optData_;
    // Функция для проверки наличия T3.
    bool isWellHorizontal() const;

   public:
    WellTrajectorySolver(){};
    /**
     * @brief Задаются гиперпараметры для оптимизатора.
     * Пока не используется.
     */
    void setPSOdata();
    /**
     * @brief Входные параметры для построения оптимальной траектории.
     * @param pInitial - 2D координаты кустовой площадки;
     * @param targets - структура для задания целей T1 и,возможно, Т3;
     * @param cs - структура для задания ограничений.
     */
    void setData(
        const Point2d& pInitial, const Target& targets,
        const WellTrajectoryConstraints& cs = WellTrajectoryConstraints());
    /**
     * @brief Парсер для входного json-файла.
     * @param filename - полный путь json-файла. 
     * Пример файла выглядит следующим образом:
     * {"Platform": {
     *  "name": "test_case_1","date": "01.01.1970",
     *  "coord": [0.0,0.0],"T1": [500.0,0.0,1000.0]},
     * "Constraints":{"MinHoldLength": 100,"MaxDLS": 1}}
     * @return false - если всё успешно,
     * true - если найдены ошибки
     * (логика такая для common_solver_wrapper'a)
     */
    bool loadInputData(const std::string& filename);
    /**
    * @brief Основная функция для вычислений.
    */
    void optimize();
    /**
    * 
    * @brief Функция, возвращающая приватное поле optData_.
    * 
    */
    PsoValueType getPSOdata() const;
    // north, east, depth
    /**
    * @brief Функция, возвращающая приватное поле pCtrajectory_.
    */
    std::vector<Eigen::Vector3d> getTrajectoryPoints() const;
    /**
    * @brief Функция, возвращающая приватное поле pInclinometry_.
    */
    std::vector<Eigen::Vector3d> getInclinometry() const;
    /**
    * @brief Функция записывает результаты в json-файл.
    * @param filename_json - название файла с результатами,
    * выглядит следующим образом:
    * {"length": 5000,   "points": "./inclinometry.csv"}
    * @param filename_inclinometry - название файла, 
    * в который будет записана лишь инклинометрия в порядке md,inc(°),azi(°)
    * 
    */
    void writeResults(const std::string& filename_json,
                      const std::string& filename_inclinometry) const;
    /**
    * @brief Функция, возвращающая приватное поле length_.
    */
    double getTrajectoryLength() const;
};
/**
* @brief Целевая функция для оптимизации траектории.
* @param well - вектор шаблонов построения;
* @param cs - ограничения на бурение;
* @param penalty - штраф за непостроение траектории.
*/
double scoreSolver(std::vector<TrajectoryTemplate*>& well,
                   const WellTrajectoryConstraints& cs, double penalty = 1000.);
}  // namespace well_trajectory
