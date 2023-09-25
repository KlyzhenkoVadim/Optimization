#include <Separation_factor.h>

double brute_force_segdist(point3d pi, point3d pf, point3d qi, point3d qf)
{
	size_t N = 10000;
	double distance = std::numeric_limits<double>::infinity();
	double dd = (pi - qi).dot(pi-qi), da = (pi-qi).dot(pf-pi), 
		db = (pi-qi).dot(qf-qi), aa = (pf-pi).dot(pf-pi),bb = (qf-qi).dot(qf-qi), ab = (pf-pi).dot(qf-qi);
	double t = 0, s = 0;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			t = double(i) / (N - 1);
			s = double(j) / (N - 1);
			distance = std::min(distance,
				dd + 2*da*t - 2*db*s + aa*t*t - 2*ab*t*s + bb*s*s);
		}
	}
	return sqrt(distance);
}

double SegmentDistance(point3d pi, point3d pf, point3d qi, point3d qf)
{
	double distance = std::numeric_limits<double>::infinity();
	double t, s; // t,s
	point3d a = pf - pi;
	point3d b = qf - qi;
	point3d delta = pi - qi;
	double a_norm = a.norm(), b_norm = b.norm(), ab = a.dot(b), da = delta.dot(a);
	double cross2 = (a.cross(b)).norm(), db = delta.dot(b);
	cross2 *= cross2;
	if (cross2 < EPSILON) // t = b_norm / a_norm * s - delta.dot(a) / a_norm / a_norm;
	{
		double intercept = da / ab;
		if (intercept >= 0)
		{
			t = 0;
			s = std::min(1., intercept);
		}
		else
		{
			s = 0;
			t = std::min(1., -intercept);
		}
		distance = (delta + a * t - b * s).norm();
	}
	else 
	{
		t = -(b_norm * b_norm * da - db * ab) / cross2;
		s = (a_norm * a_norm * db - ab * da) / cross2;
		if (t * (t - 1) < EPSILON and s * (s - 1) < EPSILON)
		{
			distance = (delta + a * t - b * s).norm();
		}
		else 
		{
			t = std::max(0., std::min(1., t));
			s = std::max(0., std::min(1., s));
			double s_star = std::max(0., std::min(1., (db + ab * t) / b_norm / b_norm));
			double t_star = std::max(0., std::min(1., (-da + ab * s) / a_norm / a_norm));
			distance = std::min((delta + a * t - b * s_star).norm(), (delta + a * t_star - b * s).norm());
		}
	}
	return distance;
}

void testSegmentsDist()
{
	/*REGIONS:
	*     s
	*	1 | 2  | 3
	*  --------------
		8 | 0  | 4
	   --------------t
		7 | 6  | 5
	*/
	std::vector<line> lines1, lines2;
	//parallel
	lines1.push_back(line({ 1,5,0 }, { 6,3,0 }));
	lines2.push_back(line({ -6,-3,0 }, { -1,-5,5 })); 
	//skew lines
	lines1.push_back(line({ 1,2,0 }, { 5,10,0 }));
	lines2.push_back(line({ 19,2,5 }, { 9,8,5 }));
	// reg0 t,s[0,1]
	lines1.push_back(line({ 0,5,0 }, { 0,2,0 }));
	lines2.push_back(line({ 0,3,0 }, { 0,0,5 }));
	// reg1 s > 1 t < 0
	lines1.push_back(line({ 0,5,0 }, {0,2,0}));
	lines2.push_back(line({0,0,5}, {0,4,2}));
	// reg5 s < 0 t > 1
	lines1.push_back(line({ 0,2,0 }, { 0,5,0 }));
	lines2.push_back(line({ 0,4,2 }, { 0,0,5 }));
	// reg3 s > 1 t > 1
	lines1.push_back(line({ 0,2,0 }, { 0,5,0 }));
	lines2.push_back(line({ 0,0,5 }, { 0,4,2 }));
	// reg7 s < 0 t < 0
	lines1.push_back(line({ 0,5,0 }, { 0,2,0 }));
	lines2.push_back(line({ 0,4,2 }, { 0, 0, 5 }));
	// reg 2 s > 1 t [0,1]
	lines1.push_back(line({0,2,0}, {0,5,0}));
	lines2.push_back(line({0,0,4}, {0,1,3}));
	// reg 6 s < 0 t [0,1]
	lines1.push_back(line({0,2,0}, {0,5,0}));
	lines2.push_back(line({0,1,3}, {0,0,4}));
	// reg 5 s [0,1] t > 1
	lines1.push_back(line({0,-2,2}, {0,1,2}));
	lines2.push_back(line({0,0,5}, {0,5,0}));
	// reg 8 s [0,1] t < 0
	lines1.push_back(line({ 0,1,2 }, { 0,-2,2 }));
	lines2.push_back(line({ 0,0,5 }, { 0,5,0 }));

	for(size_t i = 0 ; i < lines1.size();++i)
		assert(abs(SegmentDistance(lines1[i].pi,lines1[i].pf, lines2[i].pi, lines2[i].pf) - 
			brute_force_segdist(lines1[i].pi, lines1[i].pf, lines2[i].pi, lines2[i].pf)) < 1e-3); 

}

double minDistPolylines(const std::vector<line>& polyline1, const std::vector<line>& polyline2)
{
	double minDist = std::numeric_limits<double>::max();
	for (size_t i = 0; i < polyline1.size(); ++i)
	{
		for (size_t j = 0; j < polyline2.size(); ++j)
		{
			minDist = std::min(minDist, SegmentDistance(polyline1[i].pi, polyline1[i].pf, polyline2[j].pi, polyline2[j].pf));
		}
	}
	return minDist;
}