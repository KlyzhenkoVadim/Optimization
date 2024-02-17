#pragma once
#include <Eigen/Dense>
#include <vector>

// все обозначения сохранены из статьи: A compendium of directional calculations
// based on minimum curvature method. Sawaryn, Thorogood.

/**
* @brief Функция вычисления касательного вектора
* @param azimuth - азимутальный угол(°);
* @param inclination - зенитный угол(°)
* @return вектор с координатами NS,EW,TVD.
*/
Eigen::Vector3d calcTangentVector(double azimuth, double inclination);
/**
* @brief Функция, обратная calcTangentVector
* @param t - касательный вектор. \|t\| = 1
* @return std::pair<double,double> - где first - зенитный угол(°), 
* second - азимутальный угол(°)
*/
std::pair<double, double> CartesianToSpherical(Eigen::Vector3d t);
/**
* @brief Выбор системы координат для интерполяции
* @param CARTESIAN - NS,EW,TVD - точка на траектории
* @param MD - MD,NS,EW,TVD - где NS,EW,TVD - касательный вектор точки.
*/
enum class CoordinateSystem { CARTESIAN, MD };

constexpr double PI = 3.14159265358979323846;
constexpr double EPSILON = 1e-9;

using interpolatedValueType = std::pair<std::vector<double>, std::vector<Eigen::Vector3d>>;

// TODO: Выделить, какие из виртуальных методов будут константными (!!!)
// (Почти все, кроме, кажется, solve() (!!!!)

/**
* @brief TrajectoryTemplate - абстрактный базовый класс шаблона траектории.
* Его наследуют конкретные реализации шаблонов Hold,Curve,CurveHold,CurveHoldCurveHold.
* Это сделано для того, чтобы представить траекторию в виде вектора из TrajectoryTemplate*
* Пример создания траектории:
* std::vector<TrajectoryTemplate*> well;
* well.push_back(new Hold(params...));
* well.push_back(new CurveHoldCurveHold(params...));
*/
class TrajectoryTemplate
{
private:
    /**
    * @brief Функция, интерполирующая дугу угла alpha 
    * с известными касательными векторами на концах.
    * @param t1,t2 - касательные вектора начала и конца дуги соответственно;
    * @param alpha - угол дуги в радианах
    * @param nums - число возвращаемых значений
    * @return first - вектор значений угла в радианах от 0 до alpha
    * @return second - вектор касательных векторов в соответствующих точках.
    */
 interpolatedValueType calcInterpolPoints(const Eigen::Vector3d& t1,
                                          const Eigen::Vector3d& t2,
                                          double alpha, size_t nums = 100);

public:
    virtual ~TrajectoryTemplate() = default;
    /**
    * @brief Виртуальный метод, проверка построился ли шаблон.
    */
    virtual int getCondition() = 0;
    /**
    * @brief Виртуальный метод для интерполяции
    * @param coordinateSystem - выбор системы координат
    */
    virtual void points(CoordinateSystem coordinateSystem) = 0;
    /**
    * @brief Виртуальный метод для длины шаблона.
    */
    virtual double length() = 0;
    /**
     * @brief Виртуальный метод для извилистости (в данном случае, это сумма
     * углов всех дуг в шаблоне измеренная в радианах.
     */
    virtual double getTortuosity() = 0;
    /**
     * @brief Виртуальный метод для получения координат точки NS,EW,TVD в
     * конкретном md.
     * @param md - натуральный (md\in[0;1]) параметр кривой.
     */
    virtual Eigen::Vector3d FunctionPoint(double md) = 0;
    /**
     * @brief Виртуальный метод для получения касательного вектора в точке на
     * кривой с конкретным параметром.
     * @param md - натуральный (md\in[0;1]) параметр кривой.
     */
    virtual Eigen::Vector3d FunctionTangent(double md) = 0;
    /**
    * @brief Виртуальные функции для получения начальной точки.
    * @param coordinateSystem - система координат.
    */
    virtual void getInitPoint(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) = 0;
    /**
    * @brief Виртуальные функции для получения цели Т1.
    * @param coordinateSystem - система координат.
    */
    virtual void getTarget1Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) = 0;
    /**
    * @brief Виртуальные функции для получения цели Т3.
    * @param coordinateSystem - система координат.
    */
    virtual void getTarget3Point(CoordinateSystem coordinateSystem = CoordinateSystem::CARTESIAN) = 0;
    std::vector<Eigen::Vector3d> pointsCartesian; // вектор  декартовых координат NS,EW,TVD
    // вектор значений MD точек вдоль шаблона и NS,EW,TVD координат касательных
    // векторов в этих точках. 4D
    std::vector<Eigen::Vector4d> pointsMD;
    // вектора в декартовых координатах начальной точки, цели Т1 и Т3
    // соответственно.
    Eigen::Vector3d pointInitial,pointT1, pointT3;
    // вектор из md, ns,ew,tvd касательного вектора в начальной точке, целях Т1
    // и Т3 соответственно.
    Eigen::Vector4d pointInitialMD, pointMDT1, pointMDT3;
    /**
     * @brief Метод для интерполяции произвольной дуги фиксированного радиуса в
     * декартовых координатах.
     * @param p1, t1 - декартовы координаты и касательный вектор начальной точки;
     * @param t2 - касательный вектор конечной точки;
     * @param R - радиус данной дуги;
     * @param alpha - угол дуги (в плоскости [t1xt2])
     * @param nums - число точек после интерполяции
     * @return вектор декартовых координат дуги размера nums.
     */
    std::vector<Eigen::Vector3d> calcInterpolCartesianPoints(const Eigen::Vector3d& p1,
        const Eigen::Vector3d& t1,
        const Eigen::Vector3d& t2,
        double R, double alpha, size_t nums = 100);
    /**
     * @brief Метод для интерполяции произвольной дуги фиксированного радиуса
     * для касательных векторов.
     * @param p1, t1 - декартовы координаты и касательный вектор начальной
     * точки;
     * @param t2 - касательный вектор конечной точки;
     * @param R - радиус данной дуги;
     * @param alpha - угол дуги (в плоскости [t1xt2])
     * @param nums - число точек после интерполяции
     * @return вектор координат md,nsTangent,ewTangent,tvdTangent дуги размера nums.
     */
    std::vector<Eigen::Vector4d> calcInterpolMDPoints(const Eigen::Vector3d& p1,
        const Eigen::Vector3d& t1,
        const Eigen::Vector3d& t2,
        double R, double alpha, size_t nums = 100);
};	

/**
 * @brief Метод, позволяющий получить длину траектории как вектора шаблонов.
 * @param Well - вектор шаблонов построения
 * @TODO: Well должен передаваться по константной ссылке(!!!!)
 */
double allLength(const std::vector<TrajectoryTemplate*>& Well);
/**
* @brief Метод, позволяющий получить декартовых координаты траектории.
* @param Well - вектор шаблонов построения
* @return вектор декартовых координат.
* @return Размер (???) если Well состоит из k шаблонов, каждый из которых имеет nums_i i = 1...k точек, то
* @return размер вектора nums_1 + . . . + nums_k.
* @TODO: Well должен передаваться по константной ссылке(!!!!)
*/
std::vector<Eigen::Vector3d> allPointsCartesian(std::vector<TrajectoryTemplate*>& Well);
/**
* @brief Метод, позволяющий получить касательные вектора вдоль траектории.
* @param Well - вектор шаблонов построения
* @return вектор касательных векторов.
* @return Размер (???) если Well состоит из k шаблонов, каждый из которых имеет nums_i i = 1...k точек, то
* @return размер вектора nums_1 + . . . + nums_k.
* @TODO: Well должен передаваться по константной ссылке(!!!!)
*/
std::vector<Eigen::Vector4d> allPointsMD(std::vector<TrajectoryTemplate* >& Well);
/**
* @brief Аналог функции FunctionPoint для конкретных шаблонов, только вдоль всей траектории.
* @param md - натуральный параметр ВСЕЙ траектории$
* @param well - траектория, вектор шаблонов построения.
* @return 3d координаты точки траектории при параметре md
* @TODO: Well должен передаваться по константной ссылке(!!!!)
*/
Eigen::Vector3d FunctionWellPoint(double md, std::vector<TrajectoryTemplate*>& well); // md[0,1]
/**
* @brief Аналог функции FunctionTangent для конкретных шаблонов, только вдоль всей траектории.
* @param md - натуральный параметр ВСЕЙ траектории
* @param well - траектория, вектор шаблонов построения.
* @return 3d координаты касательного вектора при параметре md (норма вектора = 1)
* @TODO: Well должен передаваться по константной ссылке(!!!!)
*/
Eigen::Vector3d FunctionWellTangent(double md, std::vector<TrajectoryTemplate*>& well); // md[0,1]
/**
* @brief Функция, выполняющая построение каждого шаблона траектории, из которого она состоит.
* @param Well - вектор шаблонов.
* @return 0, если траекторию по заданным параметрам (то есть все шаблоны) можно построить;
* @return -i, если i шаблонов траектории не построены.
*/
int solve(std::vector<TrajectoryTemplate*>& Well);
