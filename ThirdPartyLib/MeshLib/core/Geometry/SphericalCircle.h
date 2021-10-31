#ifndef _SPHERICALCIRCLE_H_
#define _SPHERICALCIRCLE_H_

#include "Point.h"
#include <Eigen/Dense>

namespace MeshLib {
	class CSphericalCircle
	{
	public:
		CSphericalCircle(CPoint center, double radius): m_center(center), m_radius(radius) {}
		CPoint& c() { return m_center; }
		double& r() { return m_radius; }
		CPoint normal() { return m_center; }
		void _spherical_to_euclidean();
	protected:
		CPoint m_center;
		double m_radius;
	};
	
	inline void CSphericalCircle::_spherical_to_euclidean() {
		double angle = m_radius;
		m_center = m_center * cos(angle);
		m_radius = sin(angle);
	}
	

	inline void _nearest_point(CPoint &p, CPoint &basepoint, CPoint &direction, CPoint &nearest, double &dist)
	{
		CPoint diff = basepoint - p;
		double t = -(diff * direction);
		nearest = basepoint + direction * t;
		dist = (nearest - p).norm(); 
	}

	inline int _spherical_circle_circle_intersection(CSphericalCircle c1, CSphericalCircle c2, CPoint &i1, CPoint &i2) 
	{
		c1._spherical_to_euclidean();
		c2._spherical_to_euclidean();
		
		// determine b
		Eigen::Vector2d b;
		b << c1.c() * c1.c(), c2.c() * c2.c();
		Eigen::MatrixXd A(2,3);
		A << c1.c()[0], c1.c()[1], c1.c()[2],
			c2.c()[0], c2.c()[1], c2.c()[2];
		Eigen::VectorXd x;
		x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
		//assert((A*x - b).norm() < 1e-7);
		CPoint basepoint(x(0), x(1), x(2));
		//assert(abs((basepoint - c1.c()) * c1.c()) < 1e-7);
		//assert(abs((basepoint - c2.c()) * c2.c()) < 1e-7);
		//basepoint = basepoint / basepoint.norm();
		CPoint intersection_vector = c1.c() ^ c2.c();
		if (intersection_vector.norm() < 1e-7) {
			return -1;
		}
		intersection_vector = intersection_vector / intersection_vector.norm();
		
		// find the nearest point to two center;
		CPoint p1;
		double dist1;
		CPoint nearest1;
		_nearest_point(c1.c(), basepoint, intersection_vector, nearest1, dist1);
		//assert(abs((nearest1 - c1.c()) * intersection_vector) < 1e-7);
		CPoint p2;	
		double dist2;
		CPoint nearest2;
		_nearest_point(c2.c(), basepoint, intersection_vector, nearest2, dist2);
		//assert(abs((nearest2 - c2.c()) * intersection_vector) < 1e-7);

		//assert((nearest1 - nearest2).norm() < 1e-7);

		if (dist1 > c1.r() || dist2 > c2.r()) {
			return -1;
		}

		double inter_dist1 = sqrt(c1.r() * c1.r() - dist1 * dist1);
		double inter_dist2 = sqrt(c2.r() * c2.r() - dist2 * dist2);
		//assert(abs(inter_dist1 - inter_dist2) < 1e-7);
		
		i1 = intersection_vector * inter_dist1 + nearest1;
		i2 = intersection_vector * (-inter_dist1) + nearest1;
		i1 = i1 / i1.norm();
		i2 = i2 / i2.norm();
		//assert(abs(i1.norm() - 1) < 1e-7);
		//assert(abs(i2.norm() - 1) < 1e-7);
		return 1;
	}
}


#endif