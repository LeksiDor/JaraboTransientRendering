#ifndef __CIRCLE_H
#define __CIRCLE_H

#include "LinearAlgebra/Vector2.h"
#include "Geometry/Object.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include <math.h>

/** 2D Circle. */
class Circle : public Object2D
{
	Vector2 center;
	Real radius;

	/*Intersection data*/
	Real cx, cy;

	void setup()
	{
		cx = center[0];
		cy = center[1];
	}

public:

	Vector2 get_normal(const Vector2& surf_point) const { return normalize(surf_point - center); }

	Circle(const Vector2& _center, const Real& _radius,
						Material<2>* _mat)
		: Object2D(_mat), center(_center), radius(_radius)
	{
		setup();
	}

	void transform(const Matrix3x3& m)
	{
		// TODO
	}

	bool intersect(Ray2D& r, Intersection2D &it, float max_distance) const
	{
		Real dx, dy, ox, oy;
		Real s1, s2, k1, k2, t1, t2, t;
		Real A, B, C, discr;
		
		//		printf("Circle::intersect(r, it)\n");
		
		dx = r.get_direction()[0];
		dy = r.get_direction()[1];
		ox = r.get_origin()[0];
		oy = r.get_origin()[1];
		
		k1 = cx*dy - ox*dy - cy*dx + oy*dx;
		k2 = radius*dy;
		
		B = -2.0f*radius*dx; //k3
		A = k1 - k2;
		C = k1 + k2;
		discr = B*B - 4.0f*A*C;
		
		if (discr < 0.0f) return false;
		
		discr = sqrt(discr);
		s1 = (-B + discr)/(2.0f*A);
		s2 = (-B - discr)/(2.0f*A);
		
		if (fabs(dx) > fabs(dy)) //choose the component of direction furthest from 0.0
		{ 
			t1 = (cx+radius*(1.0f-s1*s1)/(1.0f+s1*s1)-ox)/dx;
			t2 = (cx+radius*(1.0f-s2*s2)/(1.0f+s2*s2)-ox)/dx;
		}
		else
		{
			t1 = (cy+radius*(2.0f*s1)/(1.0f+s1*s1)-oy)/dy;
			t2 = (cy+radius*(2.0f*s2)/(1.0f+s2*s2)-oy)/dy;
		}
		
		if (t1 < 0.0f && t2 < 0.0f) return false; //circle is behind the ray
		else if (t1 < 0.0f && t2 > 0.0f) t = t2; //t1 behind, t2 in front
		else if (t1 > 0.0f && t2 < 0.0f) t = t1; //t1 in front, t2 behind
		else t = std::min(t1,t2); //both in front, choose closer
		
		if (t < max_distance && r.cond_set_parameter(t)) 
		{
			it.set(r, this, get_normal(r.get_position()), Vector2(0.,0.));
			return true;
		}
		return false;
	}
	
	bool intersect(const Ray2D& r, float max_distance) const
	{
		Real dx, dy, ox, oy;
		Real s1, s2, k1, k2, t1, t2, t;
		Real A, B, C, discr;
		
		dx = r.get_direction()[0];
		dy = r.get_direction()[1];
		ox = r.get_origin()[0];
		oy = r.get_origin()[1];
		
		k1 = cx*dy - ox*dy - cy*dx + oy*dx;
		k2 = radius*dy;
		
		B = -2.0f*radius*dx; //k3
		A = k1 - k2;
		C = k1 + k2;
		discr = B*B - 4.0f*A*C;
		
		if (discr < 0.0f) return false;
		
		discr = sqrt(discr);
		s1 = (-B + discr)/(2.0f*A);
		s2 = (-B - discr)/(2.0f*A);
		
		if (fabs(dx) > fabs(dy)) //choose the component of direction furthest from 0.0
		{ 
			t1 = (cx+radius*(1.0f-s1*s1)/(1.0f+s1*s1)-ox)/dx;
			t2 = (cx+radius*(1.0f-s2*s2)/(1.0f+s2*s2)-ox)/dx;
		}
		else
		{
			t1 = (cy+radius*(2.0f*s1)/(1.0f+s1*s1)-oy)/dy;
			t2 = (cy+radius*(2.0f*s2)/(1.0f+s2*s2)-oy)/dy;
		}
		
		if (t1 < 0.0f && t2 < 0.0f) return false; //circle is behind the ray
		else if (t1 < 0.0f && t2 > 0.0f) t = t2; //t1 behind, t2 in front
		else if (t1 > 0.0f && t2 < 0.0f) t = t1; //t1 in front, t2 behind
		else t = std::min(t1,t2); //both in front, choose closer
		
		if (t < SMALLEST_DIST || t > max_distance || r.get_parameter() < t) 
			return false;
		
		return true;
	}
	
	AABB2D get_bounding_box() const
	{
		return AABB2D(Vector2(cx-radius,cy-radius), Vector2(cx+radius,cy+radius));
	}

	Vector2 get_center() const
	{
		return center;
	}
	
	Real get_intersection_cost()const
	{
		return 0;
	}
	
	virtual void get_lambertiancoords(const VectorN<1>& c, Vector2& pos, Vector2& n)
	{ //c range is [0,1]
		Real angle;
		angle = c[0]*2*M_PI;
		pos = center + radius*Vector2(sin(angle), cos(angle));
		n = get_normal(pos);
	}
};
#endif
