/* 
	Created by:			Julio Marco 
	Date of creation:	March 2013
 */

#ifndef __RECTANGLE_H
#define __RECTANGLE_H

#include "LinearAlgebra/Vector2.h"
#include "Geometry/Object.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"

/** 2D axis-aligned rectangle 
 */
class pRectangle : public Object2D
{
	Vector2 _min, _max; //max and min corners
	Vector2 normal;
	
	/* Normals */
	Vector2 n1, n2, n3, n4; // normals of 4
	
	void setup()
	{
		n1 = Vector2(-1., 0.);
		n2 = Vector2(0., 1.);
		n4 = -n2;
		n3 = -n1;
	}

public:
	Vector2 get_normal(const Vector2& surf_point) const 
	{
		Real eps = 1e-3;
		//Check which segment a point belongs to and return normal
		if (fabs(surf_point[0] - _min[0]) < eps) 
			return n1;
		else if (fabs(surf_point[0] -_max[0]) < eps) 
			return n3;
		else if (fabs(surf_point[1] - _min[1]) < eps) 
			return n4;
		else //surf_point[1] == _max[1]
			return n2;
	}

	pRectangle(const Vector2& p0, const Vector2& p1,
						Material<2>* _mat)
		: Object2D(_mat), 
		  _min(Vector2(std::min(p0[0], p1[0]), std::min(p0[1], p1[1]))), 
		  _max(Vector2(std::max(p0[0], p1[0]), std::max(p0[1], p1[1]))) 
	{
		setup();
	}

	void transform(const Matrix3x3& m)
	{
		// TODO
	}

	bool intersect(Ray2D& r, Intersection2D &it, float max_distance) const
	{
		Real t;
		AABB2D bb = get_bounding_box();

		if (!bb.intersect(r, t)) return false;
		
		if (t < max_distance && r.cond_set_parameter(t)){
			it.set(r, this, get_normal(r.get_position()), Vector2(0.,0.));
			return true;
		}
		return false;
	}
	
	bool intersect(const Ray2D& r, float max_distance) const
	{
		Real t;
		AABB2D bb = get_bounding_box();
		if (r.get_origin()[0] > _min[0] && r.get_origin()[0] < _max[0] && r.get_origin()[1] > _min[1] && r.get_origin()[1] < _max[1])
			return true;
		
		if (!bb.intersect(r, t) || t < SMALLEST_DIST || t > max_distance || r.get_parameter() < t) 
			return false;
		
		return true;
	}

	AABB2D get_bounding_box() const
	{
		return AABB2D(_min, _max);
	}
	
	Vector2 get_center() const
	{
		return 0.5*(_min+_max);
	}
	
	Real get_intersection_cost() const
	{
		return 0;
	}
};
#endif
