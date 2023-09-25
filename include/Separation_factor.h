#pragma once
#include <TrajectoryTemplate.h>

using point3d = Eigen::Vector3d;

struct line
{
	point3d pi, pf;
};

double SegmentDistance(point3d pi, point3d pf, point3d qi, point3d qf);

void testSegmentsDist();
