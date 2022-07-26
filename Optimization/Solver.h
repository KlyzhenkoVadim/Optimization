#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>
/*
struct GeoPoint {
	double northT1;
	double eastT1;
	double tvdT1;

	double H;

	double northT3;
	double eastT3;
	double tvdT3;

	double frGas;
	double frLiq;

	std::string name;
	std::string date;
};

struct Point2d
{
	double east;
	double north;
};

struct WellPad
{
	bool isWellPad;

	Point2d coord; // координаты центра кустовой площадки

	std::vector<GeoPoint> geoAims; // координаты геологических целей
};

enum class LineMode { VERTICAL, HORIZONTAL };

using comparatorType = std::function<bool(const GeoPoint&, const GeoPoint&)>;

class Solver
{
private:
	size_t maxWellNum;
	double maxRadius;
	std::vector<GeoPoint> geoAims;
	std::vector<WellPad> wellPads;
	std::vector<std::pair<Point2d, Point2d>> lines;
	double getDistance2d(double x1, double y1, double x2, double y2);
	bool checkAdjGeoAimsSize();
	comparatorType getComparator(LineMode lineMode);
public:
	Solver(size_t maxWellNum, double maxRadius);
	void loadGeoAimFromCSV(const std::string& filename, char delim = ';');
	void randomGen(size_t count);
	void structedGen(size_t count, double step);
	void findWellPads(WellPad& wellPad, LineMode lineMode);
	void deleteSingleWellPads();
	void kMean();
	void solve();
	double getCost();
	void save2File(const std::string& filename);
};
*/